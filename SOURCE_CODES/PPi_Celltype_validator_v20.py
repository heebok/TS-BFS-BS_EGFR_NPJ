#!/usr/bin/env python3
"""
PPI Cell-Type Expression Validator (v20)
==========================================
Validates whether each protein-protein interaction (PPI) pair in BFS Beam Search pathway data
can be realized in a specific cell type.

NEW IN v20 — REMOVE BUILTIN_ENSG + UNIFY CELL-LINE LOOKUP:
  Per the user's request: "이 코드들을 모두 없애버리는 방법이 앞으로 있을 문제를
  줄이지 않을까". Two structural cleanups that eliminate entire classes of bugs:

  1) BUILTIN_ENSG removed entirely. Every gene-symbol → ENSG resolution now goes
     through MyGene.info, with a small in-memory cache so repeated lookups in
     one session are still O(1). This permanently fixes the wrong-ENSG bugs
     (SOS1, and any others we hadn't discovered yet) — there are simply no
     hardcoded mappings that can drift from authoritative sources.

  2) Cell-line / cancer-group / single-cell lookup UNIFIED. Previously, the
     validation engine and the HTML "Expression matching diagnostics" table
     each had their own copy of matching code. The diagnostics table lacked
     family fallback (HEK293T → HEK293), so even when the engine resolved
     correctly, the diagnostics still reported "no data". Now both call the
     same _lookup_individual_cell_line / _lookup_cancer_category /
     _lookup_single_cell methods — one source of truth, no drift.

NEW IN v19 — Wrong-ENSG name-mismatch detection + true OR-best bulk merge.
NEW IN v18 — Cell-line family fallback (HEK293T → HEK293).
NEW IN v17 — SQLite cache thread-safe; path_grade WORST edge; KPI cards split.
NEW IN v16 — BULK TSV CACHE.
NEW IN v15 — Self-healing ENSG via MyGene.info.
NEW IN v14 — HPA JSON key fixes.
NEW IN v13 — Cancer cell line CATEGORY drop-down + OR-search.

APIs used: MyGene.info, Human Protein Atlas, STRING REST v12.x.
Output: interactive HTML report
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import json
import re
import os
import time
import urllib.request
import urllib.parse
import sqlite3
import zipfile
import io
import shutil
from datetime import datetime
from collections import defaultdict

# ─── Constants ────────────────────────────────────────────────────────────────────
APP_VERSION = "2.9.0"
STRING_API = "https://version-12-0.string-db.org/api/json"

# HPA API — correct URL and column structure
# rnascm = RNA single cell type specific nTPM (core column)
# rnasctm = RNA single cell type Tau score
# rnaclcss / rnaclcsm = RNA cell line cancer category specificity / nTPM
HPA_BASE = "https://www.proteinatlas.org"
HPA_JSON_ENSG = "https://www.proteinatlas.org/{ensg}.json"
HPA_JSON_SYM = "https://www.proteinatlas.org/{symbol}.json"
HPA_SEARCH = ("https://www.proteinatlas.org/api/search_download.php"
              "?search={symbol}&format=json"
              "&columns=g,eg,rnascm,rnasctm,subcell_location,rnaclcss,rnaclcsm,rnaclct"
              "&compress=no")

# v16 — Bulk TSV downloads (the full nTPM matrix that ENSG.json doesn't always include)
# Format: gene-major long-form rows with Ensembl gene, gene name, group/cell line, TPM, pTPM, nTPM.
HPA_BULK_CANCER_GROUP_URL = "https://www.proteinatlas.org/download/tsv/rna_cell_line_cancer.tsv.zip"
HPA_BULK_CELL_LINE_URL    = "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip"

# Cache lives in a hidden user-home directory (cross-platform safe)
def _default_cache_dir() -> str:
    return os.path.join(os.path.expanduser("~"), ".ppi_validator")

HPA_CACHE_DIR = _default_cache_dir()
HPA_CACHE_DB  = os.path.join(HPA_CACHE_DIR, "hpa_cache.db")

# MyGene.info (Gene Symbol → Ensembl ID)
MYGENE_API = "https://mygene.info/v3/query?q=symbol:{sym}&species=human&fields=ensembl.gene,name"

STRING_SCORE_THRESHOLD = 400  # medium confidence


# ─── v16 Bulk HPA TSV cache ──────────────────────────────────────────────────
# Downloads HPA's two cancer cell-line TSV files once, parses them into a local
# SQLite database, and serves expression lookups directly from disk afterwards.
# This eliminates the "Low cancer specificity → empty single-entry JSON" problem
# because the bulk files contain values for *every* gene in *every* group/cell line,
# regardless of how HPA classifies the gene's specificity.
class BulkHPACache:
    """Local SQLite cache of HPA bulk cell-line TSV data.

    Schema:
      cancer_group(ensg TEXT, gene TEXT, category TEXT, ntpm REAL,  PRIMARY KEY(ensg, category))
      cell_line   (ensg TEXT, gene TEXT, cell_line TEXT, ntpm REAL, PRIMARY KEY(ensg, cell_line))
      meta        (key TEXT PRIMARY KEY, value TEXT)
    """

    SCHEMA = """
        CREATE TABLE IF NOT EXISTS cancer_group(
            ensg TEXT NOT NULL, gene TEXT NOT NULL, category TEXT NOT NULL, ntpm REAL,
            PRIMARY KEY(ensg, category)
        );
        CREATE INDEX IF NOT EXISTS idx_cg_gene ON cancer_group(gene);
        CREATE TABLE IF NOT EXISTS cell_line(
            ensg TEXT NOT NULL, gene TEXT NOT NULL, cell_line TEXT NOT NULL, ntpm REAL,
            PRIMARY KEY(ensg, cell_line)
        );
        CREATE INDEX IF NOT EXISTS idx_cl_gene ON cell_line(gene);
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT);
    """

    def __init__(self, db_path: str = HPA_CACHE_DB, log_fn=None):
        self.db_path = db_path
        self.log = log_fn or (lambda msg: None)
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        # v17: Don't keep a long-lived connection. SQLite connections are
        # bound to the thread that created them, so reusing one across
        # GUI thread + worker threads silently fails. Each operation
        # opens a fresh connection — for read-only lookups this is fast.
        self._init_schema_once()

    def _open(self):
        """Open a fresh thread-local SQLite connection."""
        conn = sqlite3.connect(self.db_path, timeout=10.0)
        return conn

    def _init_schema_once(self):
        """Ensure the schema exists. Safe to call from any thread."""
        try:
            conn = self._open()
            try:
                conn.executescript(self.SCHEMA)
                conn.commit()
            finally:
                conn.close()
        except Exception as e:
            self.log(f"⚠️ HPA cache schema init warning: {e}")

    def close(self):
        """No-op — connections are now per-call. Kept for API compatibility."""
        pass

    # ── Status checks ────────────────────────────────────────────────
    def is_populated(self) -> bool:
        """Returns True iff both tables already have rows."""
        try:
            conn = self._open()
            try:
                cur = conn.cursor()
                cur.execute("SELECT COUNT(*) FROM cancer_group")
                cg_n = cur.fetchone()[0]
                cur.execute("SELECT COUNT(*) FROM cell_line")
                cl_n = cur.fetchone()[0]
                return cg_n > 0 and cl_n > 0
            finally:
                conn.close()
        except Exception as e:
            self.log(f"⚠️ HPA cache is_populated check failed: {e}")
            return False

    def get_meta(self, key: str) -> str | None:
        try:
            conn = self._open()
            try:
                cur = conn.cursor()
                cur.execute("SELECT value FROM meta WHERE key=?", (key,))
                row = cur.fetchone()
                return row[0] if row else None
            finally:
                conn.close()
        except Exception:
            return None

    def set_meta(self, key: str, value: str):
        try:
            conn = self._open()
            try:
                cur = conn.cursor()
                cur.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)", (key, value))
                conn.commit()
            finally:
                conn.close()
        except Exception as e:
            self.log(f"⚠️ HPA cache set_meta failed: {e}")

    def status_summary(self) -> str:
        try:
            conn = self._open()
            try:
                cur = conn.cursor()
                cur.execute("SELECT COUNT(DISTINCT ensg), COUNT(*) FROM cancer_group")
                cg_genes, cg_rows = cur.fetchone()
                cur.execute("SELECT COUNT(DISTINCT ensg), COUNT(*) FROM cell_line")
                cl_genes, cl_rows = cur.fetchone()
            finally:
                conn.close()
            built = self.get_meta("built_at") or "unknown"
            ver   = self.get_meta("hpa_version") or "unknown"
            return (f"HPA cache: {cg_genes} genes × {cg_rows//max(cg_genes,1)} cancer groups | "
                    f"{cl_genes} genes × {cl_rows//max(cl_genes,1)} individual cell lines | "
                    f"built {built} (HPA {ver})")
        except Exception as e:
            return f"HPA cache: not built yet ({e})"

    # ── Download + parse + ingest ───────────────────────────────────
    def build(self, progress_fn=None) -> bool:
        """Download both TSV zip files and populate SQLite. Returns True on success."""
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        try:
            self.log("📥 Downloading HPA cancer cell-line group TSV (~2-10 MB compressed)...")
            cg_rows = self._download_and_parse_tsv(
                HPA_BULK_CANCER_GROUP_URL,
                expected_columns=("Gene", "Gene name", "Cancer", "TPM", "pTPM", "nTPM"),
                key_col="Cancer",
                progress_fn=progress_fn, progress_label="cancer groups",
            )
            self.log(f"✅ Cancer-group rows: {len(cg_rows):,}")

            self.log("📥 Downloading HPA cell-line TSV (~10-30 MB compressed)...")
            cl_rows = self._download_and_parse_tsv(
                HPA_BULK_CELL_LINE_URL,
                expected_columns=("Gene", "Gene name", "Cell line", "TPM", "pTPM", "nTPM"),
                key_col="Cell line",
                progress_fn=progress_fn, progress_label="cell lines",
            )
            self.log(f"✅ Cell-line rows: {len(cl_rows):,}")

            self.log("💾 Writing SQLite cache...")
            conn = self._open()
            try:
                cur = conn.cursor()
                # Wipe old rows
                cur.execute("DELETE FROM cancer_group")
                cur.execute("DELETE FROM cell_line")
                cur.executemany(
                    "INSERT OR REPLACE INTO cancer_group(ensg, gene, category, ntpm) VALUES (?,?,?,?)",
                    cg_rows
                )
                cur.executemany(
                    "INSERT OR REPLACE INTO cell_line(ensg, gene, cell_line, ntpm) VALUES (?,?,?,?)",
                    cl_rows
                )
                conn.commit()
            finally:
                conn.close()
            self.set_meta("built_at", datetime.now().strftime("%Y-%m-%d %H:%M"))
            self.set_meta("hpa_version", "v25 or compatible")
            self.log(f"✅ Cache built: {self.status_summary()}")
            return True
        except Exception as e:
            self.log(f"❌ Cache build failed: {e}")
            import traceback; self.log(traceback.format_exc())
            return False

    @staticmethod
    def _download_and_parse_tsv(url: str, expected_columns: tuple,
                                  key_col: str, progress_fn=None,
                                  progress_label: str = "") -> list:
        """Download a .tsv.zip from HPA, parse it, and return a list of
        (ensg, gene, key_value, ntpm) tuples."""
        req = urllib.request.Request(url, headers={"User-Agent": "PPIValidator/1.0"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            zip_bytes = resp.read()

        rows = []
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
            tsv_names = [n for n in zf.namelist() if n.endswith(".tsv")]
            if not tsv_names:
                raise RuntimeError(f"No .tsv inside {url}")
            tsv_name = tsv_names[0]
            with zf.open(tsv_name) as f:
                # Decode line-by-line — these files can be 50+ MB uncompressed
                header_line = f.readline().decode("utf-8", errors="replace").rstrip("\r\n")
                headers = header_line.split("\t")
                # Resolve column indices robustly (column order can vary by HPA version)
                col_idx = {h: i for i, h in enumerate(headers)}
                missing = [c for c in expected_columns if c not in col_idx]
                if missing:
                    raise RuntimeError(
                        f"TSV header missing columns {missing} in {tsv_name}; got {headers}")
                idx_ensg = col_idx["Gene"]            # Ensembl ID column is 'Gene'
                idx_name = col_idx["Gene name"]
                idx_key  = col_idx[key_col]
                idx_ntpm = col_idx["nTPM"]

                count = 0
                for raw in f:
                    line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) <= max(idx_ensg, idx_name, idx_key, idx_ntpm):
                        continue
                    ensg = parts[idx_ensg].strip()
                    gene = parts[idx_name].strip()
                    key  = parts[idx_key].strip()
                    ntpm_raw = parts[idx_ntpm].strip()
                    try:
                        ntpm = float(ntpm_raw) if ntpm_raw else 0.0
                    except ValueError:
                        ntpm = 0.0
                    if ensg and gene and key:
                        rows.append((ensg, gene, key, ntpm))
                    count += 1
                    if progress_fn and (count & 0xFFFF) == 0:
                        progress_fn(count, label=progress_label)
        return rows

    # ── Lookup ────────────────────────────────────────────────────
    def get_cancer_groups(self, ensg: str = "", gene: str = "") -> dict:
        """Return {category_lower: nTPM_float} for a given gene (by ENSG or symbol).
        Empty dict if no rows."""
        try:
            conn = self._open()
            try:
                cur = conn.cursor()
                if ensg:
                    cur.execute("SELECT category, ntpm FROM cancer_group WHERE ensg=?", (ensg,))
                elif gene:
                    cur.execute("SELECT category, ntpm FROM cancer_group WHERE gene=?", (gene,))
                else:
                    return {}
                return {cat.lower(): float(n) for cat, n in cur.fetchall()}
            finally:
                conn.close()
        except Exception:
            return {}

    def get_cell_lines(self, ensg: str = "", gene: str = "") -> dict:
        try:
            conn = self._open()
            try:
                cur = conn.cursor()
                if ensg:
                    cur.execute("SELECT cell_line, ntpm FROM cell_line WHERE ensg=?", (ensg,))
                elif gene:
                    cur.execute("SELECT cell_line, ntpm FROM cell_line WHERE gene=?", (gene,))
                else:
                    return {}
                return {cl.lower(): float(n) for cl, n in cur.fetchall()}
            finally:
                conn.close()
        except Exception:
            return {}


# Module-level singleton — initialized on first use, populated lazily.
_BULK_CACHE: BulkHPACache | None = None


def get_bulk_cache(log_fn=None) -> BulkHPACache:
    """Get the shared BulkHPACache instance, creating it if needed."""
    global _BULK_CACHE
    if _BULK_CACHE is None:
        _BULK_CACHE = BulkHPACache(log_fn=log_fn)
    elif log_fn is not None:
        _BULK_CACHE.log = log_fn
    return _BULK_CACHE


# ─── v20: BUILTIN_ENSG REMOVED ────────────────────────────────────────────────
# Per the user's request, all hardcoded gene-symbol → ENSG mappings have been
# removed. They were a perpetual source of subtle bugs (e.g. SOS1 was wrong,
# pointing at a different gene; SHC1 was wrong before that). Any entry that
# silently drifts from the authoritative source produces foreign-gene data
# without obvious symptoms.
#
# All resolution now goes through MyGene.info — slower on first lookup
# (~200ms) but always authoritative. The session-level `ensembl_cache` dict
# means each gene is resolved once per validation run, so the cost is paid
# at most ~30 times for a typical pathway analysis.
#
# This block is intentionally left empty (rather than deleted entirely) to
# preserve the BUILTIN_ENSG name so any external test or import continues
# to work; it simply contains no entries.
BUILTIN_ENSG: dict = {}

# Major cell types (HPA reference)
CELL_TYPES = [
    "(Other — use manual field below)",  # v17: explicit opt-out for manual entry
    # ── Major cancer cell lines (most-requested in literature) ──
    "MDA-MB-468 (cancer cell line)", "MDA-MB-231 (cancer cell line)",
    "MCF7 (cancer cell line)", "T-47D (cancer cell line)",
    "HeLa (cancer cell line)", "A-549 (cancer cell line)",
    "HEK293T (cancer cell line)", "HEK293 (cancer cell line)",
    "HCT116 (cancer cell line)", "HT-29 (cancer cell line)",
    "K-562 (cancer cell line)", "Jurkat (cancer cell line)",
    "U-2 OS (cancer cell line)", "U-87 MG (cancer cell line)",
    "HepG2 (cancer cell line)", "Huh-7 (cancer cell line)",
    "PC-3 (cancer cell line)", "DU-145 (cancer cell line)",
    "SH-SY5Y (cancer cell line)", "Caco-2 (cancer cell line)",
    "A375 (cancer cell line)", "PANC-1 (cancer cell line)",
    "BT-474 (cancer cell line)", "SK-BR-3 (cancer cell line)",
    "Saos-2 (cancer cell line)",
    # ── Primary single-cell types (HPA single-cell type table) ──
    "cardiomyocytes", "skeletal muscle cells", "smooth muscle cells",
    "hepatocytes", "cholangiocytes", "neurons", "astrocytes", "oligodendrocytes",
    "microglia", "endothelial cells", "fibroblasts", "adipocytes",
    "T cells", "B cells", "NK cells", "monocytes", "macrophages",
    "dendritic cells", "neutrophils", "erythroid cells",
    "keratinocytes", "melanocytes",
    "proximal tubular cells", "podocytes",
    "alveolar cells", "club cells",
    "intestinal epithelial cells", "goblet cells",
    "pancreatic acinar cells", "beta cells",
    "breast glandular cells", "cervical cells",
    "spermatocytes", "oocytes",
]

# ─── HPA Cancer Cell Line Categories (v13 NEW) ──────────────────────────────
# 30 categories used by HPA's "RNA cancer cell line" classification.
# Reference: https://www.proteinatlas.org/humanproteome/cell+line
CANCER_CELL_LINE_CATEGORIES = [
    "(none — use cell line only)",
    "Adrenocortical cancer",
    "Bile duct cancer",
    "Bladder cancer",
    "Bone cancer",
    "Brain cancer",
    "Breast cancer",
    "Cervical cancer",
    "Colorectal cancer",
    "Esophageal cancer",
    "Gallbladder cancer",
    "Gastric cancer",
    "Head and neck cancer",
    "Kidney cancer",
    "Leukemia",
    "Liver cancer",
    "Lung cancer",
    "Lymphoma",
    "Myeloma",
    "Neuroblastoma",
    "Non-cancerous",
    "Ovarian cancer",
    "Pancreatic cancer",
    "Prostate cancer",
    "Rhabdoid",
    "Sarcoma",
    "Skin cancer",
    "Testis cancer",
    "Thyroid cancer",
    "Uncategorized",
    "Uterine cancer",
]

# Suggested default cancer category for each common HPA cell line.
# Used to auto-suggest a category when the user picks a known cell line.
CELL_LINE_DEFAULT_CATEGORY = {
    "hela":         "Cervical cancer",
    "siha":         "Cervical cancer",
    "mda-mb-468":   "Breast cancer",
    "mcf-7":        "Breast cancer",
    "mcf7":         "Breast cancer",
    "mda-mb-231":   "Breast cancer",
    "t-47d":        "Breast cancer",
    "t47d":         "Breast cancer",
    "bt-474":       "Breast cancer",
    "skbr3":        "Breast cancer",
    "a-549":        "Lung cancer",
    "a549":         "Lung cancer",
    "h1299":        "Lung cancer",
    "h460":         "Lung cancer",
    "calu-3":       "Lung cancer",
    "pc-9":         "Lung cancer",
    "hcc827":       "Lung cancer",
    "a-431":        "Skin cancer",
    "a375":         "Skin cancer",
    "a-375":        "Skin cancer",
    "sk-mel-30":    "Skin cancer",
    "sk-mel-28":    "Skin cancer",
    "g-361":        "Skin cancer",
    "hek293":       "Non-cancerous",
    "hek293t":      "Non-cancerous",
    "htert-rpe1":   "Non-cancerous",
    "huvec":        "Non-cancerous",
    "bj":           "Non-cancerous",
    "hap1":         "Leukemia",
    "k-562":        "Leukemia",
    "k562":         "Leukemia",
    "thp-1":        "Leukemia",
    "hl-60":        "Leukemia",
    "jurkat":       "Leukemia",
    "nb4":          "Leukemia",
    "rs4;11":       "Leukemia",
    "reh":          "Leukemia",
    "hdlm-2":       "Lymphoma",
    "raji":         "Lymphoma",
    "ramos":        "Lymphoma",
    "u-251mg":      "Brain cancer",
    "u-251 mg":     "Brain cancer",
    "u251":         "Brain cancer",
    "u-87mg":       "Brain cancer",
    "gamg":         "Brain cancer",
    "sh-sy5y":      "Neuroblastoma",
    "sk-n-be(2)":   "Neuroblastoma",
    "sk-n-sh":      "Neuroblastoma",
    "kelly":        "Neuroblastoma",
    "imr-32":       "Neuroblastoma",
    "caco-2":       "Colorectal cancer",
    "ht-29":        "Colorectal cancer",
    "sw480":        "Colorectal cancer",
    "hct-116":      "Colorectal cancer",
    "rko":          "Colorectal cancer",
    "lovo":         "Colorectal cancer",
    "hep-g2":       "Liver cancer",
    "hepg2":        "Liver cancer",
    "huh-7":        "Liver cancer",
    "snu-449":      "Liver cancer",
    "panc-1":       "Pancreatic cancer",
    "miapaca-2":    "Pancreatic cancer",
    "bxpc-3":       "Pancreatic cancer",
    "aspc-1":       "Pancreatic cancer",
    "rt-4":         "Bladder cancer",
    "t24":          "Bladder cancer",
    "5637":         "Bladder cancer",
    "pc-3":         "Prostate cancer",
    "lncap":        "Prostate cancer",
    "du145":        "Prostate cancer",
    "du-145":       "Prostate cancer",
    "22rv1":        "Prostate cancer",
    "ovcar-3":      "Ovarian cancer",
    "skov3":        "Ovarian cancer",
    "sk-ov-3":      "Ovarian cancer",
    "efo-21":       "Ovarian cancer",
    "a-204":        "Sarcoma",
    "u2os":         "Bone cancer",
    "u-2os":        "Bone cancer",
    "saos-2":       "Bone cancer",
    "rh30":         "Sarcoma",
    "susa":         "Testis cancer",
    "ntera-2":      "Testis cancer",
    "tcam-2":       "Testis cancer",
    "oe19":         "Esophageal cancer",
    "kyse-150":     "Esophageal cancer",
    "ags":          "Gastric cancer",
    "snu-1":        "Gastric cancer",
    "snu-16":       "Gastric cancer",
    "kato-iii":     "Gastric cancer",
    "tt":           "Thyroid cancer",
    "ftc-133":      "Thyroid cancer",
    "8505c":        "Thyroid cancer",
    "h295r":        "Adrenocortical cancer",
    "rpmi-8226":    "Myeloma",
    "u266":         "Myeloma",
    "mm.1s":        "Myeloma",
    "fadu":         "Head and neck cancer",
    "scc-25":       "Head and neck cancer",
    "detroit-562":  "Head and neck cancer",
}


# ─── v18 NEW: Cell-line family fallback ──────────────────────────────────────
# When the HPA bulk cell-line table doesn't contain the user's exact cell line,
# we look through "family" variants. The match with the highest nTPM wins
# (true OR-search across closely-related lines, per user request:
#  "HEK293T에서 없으면 HEK293에서도 찾도록 해주세요. 같은 세포에서 출발한 것").
#
# Each entry maps a normalized lookup key (lowercase, no dashes/spaces) to a
# list of normalized aliases to try IN ORDER. The first list item is the
# "preferred" form (closest to the user's query); subsequent items are looser
# fallbacks.
#
# Keys/values are normalized via _norm_cell_line() — both this map and the
# HPA cache lookup pass through that normalization, so "HEK293T", "hek 293t",
# "HEK-293T" all collapse to the same key.
CELL_LINE_FAMILIES = {
    # HEK293 family — the user's primary example
    "hek293t":      ["hek293t", "hek293", "293t", "293"],
    "hek293":       ["hek293", "hek293t", "293", "293t"],
    "hek293ad":     ["hek293ad", "hek293", "hek293t"],
    "293t":         ["293t", "hek293t", "hek293", "293"],
    "293":          ["293", "hek293", "hek293t", "293t"],

    # MDA-MB breast cancer family
    "mdamb468":     ["mdamb468", "mdamb231"],   # 468 first; 231 as backup
    "mdamb231":     ["mdamb231", "mdamb468"],
    "mdamb453":     ["mdamb453", "mdamb468", "mdamb231"],
    "mdamb436":     ["mdamb436", "mdamb231", "mdamb468"],
    "mdamb415":     ["mdamb415", "mdamb231", "mdamb468"],

    # MCF family
    "mcf7":         ["mcf7"],   # MCF7 is canonical; MCF-7 normalizes to same
    "mcf10a":       ["mcf10a"],

    # HCT colon cancer family
    "hct116":       ["hct116", "hct15", "hct8"],
    "hct15":        ["hct15", "hct116"],
    "ht29":         ["ht29", "hct116"],

    # K-562 / HL-60 leukemia
    "k562":         ["k562"],
    "hl60":         ["hl60"],

    # HepG2 liver
    "hepg2":        ["hepg2", "hep3b", "huh7"],
    "huh7":         ["huh7", "hepg2"],

    # Glioma U-87
    "u87mg":        ["u87mg", "u87"],
    "u87":          ["u87", "u87mg"],
    "u251mg":       ["u251mg", "u251"],
    "u251":         ["u251", "u251mg"],

    # SK-BR-3 / SKBR3 — same line, alternate spelling
    "skbr3":        ["skbr3"],

    # T-47D
    "t47d":         ["t47d"],

    # SH-SY5Y
    "shsy5y":       ["shsy5y", "shsy5", "sy5y"],

    # Caco-2
    "caco2":        ["caco2"],

    # PC-3 / DU-145 prostate
    "pc3":          ["pc3", "du145"],
    "du145":        ["du145", "pc3"],
}


def _norm_cell_line(name: str) -> str:
    """Normalize cell-line names for matching.
    Removes dashes, spaces, slashes, dots; lowercases; strips common suffixes."""
    if not name:
        return ""
    n = name.lower().strip()
    for suffix in [" (cancer cell line)", " cell line", " cells", " cell"]:
        if n.endswith(suffix):
            n = n[: -len(suffix)].strip()
    return n.replace("-", "").replace(" ", "").replace("/", "").replace(".", "")


def get_cell_line_family(cell_type: str) -> list:
    """Return ordered list of normalized aliases to try for a given cell type.
    First item is always the user's exact (normalized) input; further items
    are family fallbacks. Empty list = no aliases known beyond the input itself."""
    norm = _norm_cell_line(cell_type)
    if not norm:
        return []
    # Always start with the user's exact normalized form
    family = [norm]
    # Append known family aliases (skip duplicates)
    for alias in CELL_LINE_FAMILIES.get(norm, []):
        if alias not in family:
            family.append(alias)
    return family


def suggest_category_for_cell_line(cell_type: str) -> str:
    """Return the suggested default cancer cell line category for a given cell line.
    Returns empty string if not known."""
    if not cell_type:
        return ""
    key = cell_type.lower()
    # Strip common suffixes
    for suffix in [" (cancer cell line)", " cell line", " cells", " cell"]:
        if key.endswith(suffix):
            key = key[: -len(suffix)].strip()
    if key in CELL_LINE_DEFAULT_CATEGORY:
        return CELL_LINE_DEFAULT_CATEGORY[key]
    # Try without dashes / spaces
    key2 = key.replace("-", "").replace(" ", "")
    for k, v in CELL_LINE_DEFAULT_CATEGORY.items():
        if k.replace("-", "").replace(" ", "") == key2:
            return v
    return ""

# Subcellular location compatibility matrix (same group = interaction possible)
COMPARTMENT_GROUPS = {
    "cytosol": ["cytosol", "cytoplasm", "cytoskeleton"],
    "nucleus": ["nucleus", "nucleoplasm", "nucleoli", "nuclear membrane", "nuclear bodies", "nuclear speckles"],
    "membrane": ["plasma membrane", "cell membrane", "vesicles", "endosomes", "lysosomes"],
    "er": ["endoplasmic reticulum", "golgi apparatus"],
    "mitochondria": ["mitochondria", "mitochondrial matrix"],
    "secreted": ["extracellular", "secreted"],
}


def get_compartment_group(loc: str) -> str:
    loc_lower = loc.lower()
    for grp, terms in COMPARTMENT_GROUPS.items():
        if any(t in loc_lower for t in terms):
            return grp
    return "other"


# ─── API utilities ─────────────────────────────────────────────────────────────
def fetch_json(url: str, timeout: int = 15) -> dict | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PPIValidator/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8"))
    except Exception as e:
        return None


def gene_symbol_to_ensembl(symbol: str, cache: dict) -> str | None:
    """Convert Gene Symbol → Ensembl ID via MyGene.info (with session cache).

    v20: BUILTIN_ENSG no longer contains hardcoded entries (see comment by the
    BUILTIN_ENSG block). All resolution flows through MyGene.info, with the
    `cache` dict providing per-session memoization so each gene is resolved
    at most once per validation run. If MyGene.info is unreachable we fall
    back to a 'SYM:NAME' marker that lets HPA's symbol-direct API still try.
    """
    if symbol in cache:
        return cache[symbol]

    # Step 1: Vestigial BUILTIN_ENSG check (now always empty — see v20 notes).
    # Kept so that anyone wanting to add an emergency offline override can do so
    # by populating the BUILTIN_ENSG dict; otherwise this is a no-op.
    builtin = BUILTIN_ENSG.get(symbol.upper())
    if builtin:
        cache[symbol] = builtin
        return builtin

    # Step 2: MyGene.info — the authoritative resolver.
    try:
        url = MYGENE_API.format(sym=urllib.parse.quote(symbol))
        data = fetch_json(url, timeout=10)
        if data and data.get("hits"):
            hit = data["hits"][0]
            ens = hit.get("ensembl", {})
            ensg = None
            if isinstance(ens, dict):
                ensg = ens.get("gene")
            elif isinstance(ens, list) and ens:
                ensg = ens[0].get("gene")
            if ensg:
                cache[symbol] = ensg
                return ensg
    except Exception:
        pass

    # Step 3: Symbol-direct fallback for HPA APIs that accept a name instead of ENSG.
    cache[symbol] = f"SYM:{symbol}"
    return f"SYM:{symbol}"


def ntpm_to_level(ntpm) -> str:
    """Convert HPA nTPM value to an expression level."""
    try:
        v = float(ntpm)
        if v >= 100: return "High"
        if v >= 10:  return "Medium"
        if v >= 1:   return "Low"
        return "Not detected"
    except Exception:
        return "Not detected"


def _coerce_float(v):
    """Convert HPA value to float — handles strings like "176.8" and None safely."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        s = v.strip()
        if not s or s.lower() in ("null", "none", "n/a", "na"):
            return None
        try:
            return float(s)
        except ValueError:
            return None
    return None


def _parse_hpa_json(data, symbol: str = "") -> dict:
    """Parse HPA JSON response — based on the ACTUAL HPA v25.0 JSON structure.

    Verified top-level keys (from https://www.proteinatlas.org/ENSG{...}.json):
      "Gene"                                : "EGFR"
      "Ensembl"                             : "ENSG00000146648"

      # Single cell type (note: nCPM, NOT nTPM!)
      "RNA single cell type specificity"    : "Cell type enhanced"
      "RNA single cell type distribution"   : "Detected in many"
      "RNA single cell type specificity score": null | float
      "RNA single cell type specific nCPM"  : {"Cytotrophoblasts": "117.2", ...}
                                              # values are STRINGS

      # Cell line CANCER GROUP (this is the v13 fix — keys are under "RNA cell line"!)
      "RNA cell line specificity"           : "Low cancer specificity"
      "RNA cell line distribution"          : "Detected in all"
      "RNA cell line specificity score"     : null | float (Tau)
      "RNA cell line specific nTPM"         : {"Lung cancer": "176.8",
                                              "Breast cancer": "210.2", ...}
                                              # CG group dict, values STRINGS

      # Subcellular location is a list at top level:
      "Subcellular location"                : ["Plasma membrane", ...]
      "Subcellular main location"           : ["Plasma membrane"]

    Output dict keeps backwards-compatible field names used elsewhere in the program.
    """
    if isinstance(data, list):
        data = data[0] if data else {}
    if not isinstance(data, dict) or not data:
        return {}

    result = {
        "name": data.get("Gene", data.get("Gene name", symbol)),
        "ensg": data.get("Ensembl", ""),
        "tau_score": None,
        "cell_specificity": data.get("RNA single cell type specificity",
                                     data.get("RNA cell type specificity", "Unknown")),
        "subcellular_locations": [],
        "single_cell_expression": {},
        "single_cell_ntpm": {},  # actually nCPM but we keep the field name
        "cell_line_expression": {},
        "cell_line_ntpm": {},
        "cancer_group_expression": {},
        "cancer_group_ntpm": {},
        "cell_line_specificity": data.get("RNA cell line specificity", "Unknown"),
        "cell_line_distribution": data.get("RNA cell line distribution", ""),
        "cell_line_tau": None,
        "raw": data,
    }

    # ── Tau scores (single-cell, cell-line) ──────────────────────────
    sc_tau = (data.get("RNA single cell type specificity score") or
              data.get("RNA single cell type Tau score") or
              data.get("Tau score"))
    result["tau_score"] = _coerce_float(sc_tau)

    cl_tau = (data.get("RNA cell line specificity score") or
              data.get("RNA cell line Tau score"))
    result["cell_line_tau"] = _coerce_float(cl_tau)

    # ── Single-cell type nCPM (NOT nTPM — this is the v12/v13 bug) ─────
    sc = (data.get("RNA single cell type specific nCPM") or
          data.get("RNA single cell type specific nTPM") or
          data.get("RNA single cell type nTPM") or
          data.get("Single cell type") or {})
    if isinstance(sc, dict):
        for ct, val in sc.items():
            ct_lower = ct.lower()
            if isinstance(val, dict):
                if "nCPM" in val: numeric = _coerce_float(val["nCPM"])
                elif "nTPM" in val: numeric = _coerce_float(val["nTPM"])
                else: numeric = None
                if "Level" in val:
                    result["single_cell_expression"][ct_lower] = val["Level"]
                    result["single_cell_ntpm"][ct_lower] = numeric
                elif numeric is not None:
                    result["single_cell_expression"][ct_lower] = ntpm_to_level(numeric)
                    result["single_cell_ntpm"][ct_lower] = numeric
            else:
                numeric = _coerce_float(val)
                if numeric is not None:
                    result["single_cell_expression"][ct_lower] = ntpm_to_level(numeric)
                    result["single_cell_ntpm"][ct_lower] = numeric

    # ── Cancer cell line GROUP nTPM ───────────────────────────────────
    # *** This is the key v14 fix ***
    # HPA stores group-level (Brain cancer / Breast cancer / Lung cancer / ...)
    # under "RNA cell line specific nTPM". v13 was looking under
    # "RNA cancer cell line specific nTPM" which does NOT exist in HPA JSON.
    cg = (data.get("RNA cell line specific nTPM") or
          data.get("RNA cancer cell line specific nTPM") or
          data.get("RNA cell line cancer specific nTPM") or
          data.get("Cell line cancer") or
          data.get("Cancer cell line") or {})
    if isinstance(cg, dict):
        for cat, val in cg.items():
            cat_lower = cat.lower()
            if isinstance(val, dict):
                numeric = _coerce_float(val.get("nTPM") or val.get("nCPM"))
                if "Level" in val:
                    result["cancer_group_expression"][cat_lower] = val["Level"]
                if numeric is not None:
                    if cat_lower not in result["cancer_group_expression"]:
                        result["cancer_group_expression"][cat_lower] = ntpm_to_level(numeric)
                    result["cancer_group_ntpm"][cat_lower] = numeric
            else:
                numeric = _coerce_float(val)
                if numeric is not None:
                    result["cancer_group_expression"][cat_lower] = ntpm_to_level(numeric)
                    result["cancer_group_ntpm"][cat_lower] = numeric

    # ── Individual cell line nTPM (only present in some HPA datasets — leave empty if absent) ────
    cl = (data.get("RNA HPA cell line specific nTPM") or
          data.get("Cell line") or {})
    if isinstance(cl, dict):
        for cl_name, val in cl.items():
            cl_lower = cl_name.lower()
            if isinstance(val, dict):
                numeric = _coerce_float(val.get("nTPM"))
                if "Level" in val:
                    result["cell_line_expression"][cl_lower] = val["Level"]
                if numeric is not None:
                    if cl_lower not in result["cell_line_expression"]:
                        result["cell_line_expression"][cl_lower] = ntpm_to_level(numeric)
                    result["cell_line_ntpm"][cl_lower] = numeric
            else:
                numeric = _coerce_float(val)
                if numeric is not None:
                    result["cell_line_expression"][cl_lower] = ntpm_to_level(numeric)
                    result["cell_line_ntpm"][cl_lower] = numeric

    # ── Subcellular locations ───────────────────────────────────────────────
    # HPA v25 returns a flat list at top level
    sub = data.get("Subcellular location", data.get("Subcellular main location", []))
    if isinstance(sub, dict):
        locs = sub.get("Approved", []) or sub.get("Predicted", []) or sub.get("Main", []) or []
        result["subcellular_locations"] = locs if isinstance(locs, list) else [locs]
    elif isinstance(sub, list):
        result["subcellular_locations"] = [s for s in sub if s]
    elif isinstance(sub, str) and sub:
        result["subcellular_locations"] = [s.strip() for s in sub.split(";") if s.strip()]

    return result


def _hpa_data_returns_wrong_gene(parsed: dict, requested_symbol: str) -> bool:
    """v19: Detect the WORST class of ENSG bug — when an ENSG points at a *different*
    gene that has its own (non-empty) HPA data.

    v15's _hpa_data_is_empty() only catches ENSGs that return blank JSON. If the
    bad ENSG happens to point at some other gene, v15 silently returns that
    foreign gene's data. v19 cross-checks the parsed name against the requested
    symbol and rejects the result if they clearly disagree.

    We allow one-letter case differences and common HGNC alias forms (e.g.
    "MCF-7" vs "MCF7"); only reject when the parsed name is a recognizably
    different gene symbol.
    """
    if not parsed or not requested_symbol:
        return False
    parsed_name = (parsed.get("name") or "").strip()
    if not parsed_name:
        return False
    # Normalize for comparison: uppercase, strip dashes/spaces.
    norm_req = requested_symbol.upper().replace("-", "").replace(" ", "")
    norm_got = parsed_name.upper().replace("-", "").replace(" ", "")
    # Exact (post-normalization) match → fine.
    if norm_req == norm_got:
        return False
    # Allow short prefix overlap for known synonyms (e.g. "ERBB2" vs "ERBB2/HER2"):
    # if requested is fully contained in parsed name, accept.
    if norm_req in norm_got or norm_got in norm_req:
        return False
    # Otherwise, the ENSG returned data for a different gene.
    return True


def _hpa_data_is_empty(parsed: dict) -> bool:
    """Decide if a parsed HPA dict is effectively empty (all expression sources are blank).
    Used by fetch_hpa_data to detect a stale/wrong ENSG → trigger re-resolution via MyGene."""
    if not parsed:
        return True
    has_sc = bool(parsed.get("single_cell_expression"))
    has_cl = bool(parsed.get("cell_line_expression"))
    has_cg = bool(parsed.get("cancer_group_expression"))
    has_locs = bool(parsed.get("subcellular_locations"))
    # Treat as empty only when EVERY source is blank (rare for any real gene)
    return not (has_sc or has_cl or has_cg or has_locs)


def _resolve_ensg_via_mygene(symbol: str) -> str | None:
    """Force a fresh ENSG lookup via MyGene.info, ignoring BUILTIN_ENSG and any cache.
    Used as a self-healing step when the cached/builtin ENSG returns empty data."""
    if not symbol:
        return None
    try:
        url = MYGENE_API.format(sym=urllib.parse.quote(symbol))
        data = fetch_json(url, timeout=10)
        if not data or not data.get("hits"):
            return None
        hit = data["hits"][0]
        ens = hit.get("ensembl", {})
        if isinstance(ens, dict):
            return ens.get("gene")
        if isinstance(ens, list) and ens:
            return ens[0].get("gene")
    except Exception:
        return None
    return None


def fetch_hpa_data(ensg_or_marker: str, cache: dict, symbol: str = "",
                    log_fn=None) -> dict:
    """HPA JSON API — 3-step URL fallback strategy with self-healing ENSG resolution.

    ensg_or_marker: ENSG ID or "SYM:GENENAME" format
    log_fn:         optional logging callback to surface ENSG repair events to the user
    """
    cache_key = ensg_or_marker
    if cache_key in cache:
        return cache[cache_key]

    # Handle SYM: marker — try direct symbol query without ENSG
    if ensg_or_marker.startswith("SYM:"):
        sym = ensg_or_marker[4:]
    else:
        sym = symbol or ensg_or_marker

    empty_result = {
        "ensg": ensg_or_marker if not ensg_or_marker.startswith("SYM:") else "",
        "name": sym, "tau_score": None,
        "cell_specificity": "Unknown",
        "subcellular_locations": [],
        "single_cell_expression": {},
        "single_cell_ntpm": {},
        "cell_line_expression": {},
        "cell_line_ntpm": {},
        "cancer_group_expression": {},
        "cancer_group_ntpm": {},
        "cell_line_specificity": "Unknown",
        "cell_line_tau": None,
        "raw": None,
    }

    def _try_urls(ensg_id: str, sym_name: str) -> dict:
        """Try ENSG → SYM → SRCH URLs in order, return first non-empty parse."""
        urls = []
        if ensg_id and not ensg_id.startswith("SYM:"):
            urls.append(("ENSG", HPA_JSON_ENSG.format(ensg=ensg_id)))
        if sym_name:
            urls.append(("SYM", HPA_JSON_SYM.format(symbol=sym_name)))
            urls.append(("SRCH", HPA_SEARCH.format(symbol=urllib.parse.quote(sym_name))))
        last_parsed = {}
        for tag, url in urls:
            data = fetch_json(url, timeout=12)
            if data:
                p = _parse_hpa_json(data, sym_name)
                if not _hpa_data_is_empty(p):
                    return p
                if not last_parsed and p.get("name"):
                    last_parsed = p  # remember a metadata-only response as last resort
        return last_parsed

    # ── First attempt: with the provided ENSG (which may come from BUILTIN_ENSG cache) ──
    parsed = _try_urls(ensg_or_marker, sym)

    # ── v19: Detect *two* failure modes from the BUILTIN_ENSG cache:
    #     (a) ENSG returns empty data (v15's case)
    #     (b) ENSG returns data, but for a DIFFERENT gene than requested. This is
    #         the SOS1 bug — v15 happily accepted whatever the wrong ENSG resolved
    #         to. v19 catches this by name-cross-checking.
    # In either case, force-resolve via MyGene.info and retry. ─────────────────
    needs_repair = False
    repair_reason = ""
    if not ensg_or_marker.startswith("SYM:") and sym:
        if _hpa_data_is_empty(parsed):
            needs_repair = True
            repair_reason = "empty data"
        elif _hpa_data_returns_wrong_gene(parsed, sym):
            needs_repair = True
            got_name = (parsed.get("name") or "?")
            repair_reason = f"returned foreign gene '{got_name}' instead of '{sym}'"

    if needs_repair:
        if log_fn:
            log_fn(f"  ⚠️  [{sym}] ENSG '{ensg_or_marker}' {repair_reason} — re-resolving via MyGene.info")
        fresh_ensg = _resolve_ensg_via_mygene(sym)
        if fresh_ensg and fresh_ensg != ensg_or_marker:
            if log_fn:
                log_fn(f"  🔄 [{sym}] MyGene.info → {fresh_ensg}, retrying HPA")
            parsed_retry = _try_urls(fresh_ensg, sym)
            # Accept the retry only if it's non-empty AND the name matches.
            retry_ok = (not _hpa_data_is_empty(parsed_retry)
                        and not _hpa_data_returns_wrong_gene(parsed_retry, sym))
            if retry_ok:
                parsed = parsed_retry
                ensg_or_marker = fresh_ensg
                if log_fn:
                    log_fn(f"  ✅ [{sym}] Recovered with corrected ENSG {fresh_ensg}")
            else:
                if log_fn:
                    log_fn(f"  ⚠️  [{sym}] MyGene retry also failed; keeping prior result with caveat")
                # If the retry didn't work AND the original was wrong-gene, we
                # blank out the foreign data to avoid contaminating downstream
                # validation. The bulk cache (gene-name-keyed) will still
                # supply correct data via name-based supplement.
                if repair_reason.startswith("returned foreign gene"):
                    parsed = {}

    if not parsed:
        cache[cache_key] = empty_result
        return empty_result

    # Merge
    final = {**empty_result, **parsed}
    if not final["ensg"] and not ensg_or_marker.startswith("SYM:"):
        final["ensg"] = ensg_or_marker
    cache[cache_key] = final
    return final


def fetch_string_scores(gene_a: str, gene_b: str, species: int = 9606) -> dict | None:
    """
    Fetch STRING evidence scores for a specific protein pair.

    Returns a dictionary with the combined score and channel-specific scores
    on the standard STRING JSON scale (0.0–1.0):
      score   = combined score
      nscore  = gene neighborhood
      fscore  = gene fusion
      pscore  = phylogenetic profile
      ascore  = coexpression
      escore  = experimental
      dscore  = database
      tscore  = text mining
    """
    url = (f"{STRING_API}/network?"
           f"identifiers={urllib.parse.quote(gene_a)}%0D{urllib.parse.quote(gene_b)}"
           f"&species={species}&required_score=0&network_type=functional"
           f"&caller_identity=PPIValidator")
    data = fetch_json(url)
    if data and isinstance(data, list):
        best_payload = None
        best_score = None
        for interaction in data:
            a_name = (interaction.get("preferredName_A", "") or "").upper()
            b_name = (interaction.get("preferredName_B", "") or "").upper()
            a_id = (interaction.get("stringId_A", "") or "").upper()
            b_id = (interaction.get("stringId_B", "") or "").upper()
            combined = interaction.get("score", None)
            if combined is None:
                continue

            matched = False
            if {a_name, b_name} == {gene_a.upper(), gene_b.upper()}:
                matched = True
            elif gene_a.upper() in a_id and gene_b.upper() in b_id:
                matched = True
            elif gene_b.upper() in a_id and gene_a.upper() in b_id:
                matched = True

            if not matched:
                continue

            try:
                combined = float(combined)
            except Exception:
                continue

            payload = {
                "score": combined,
                "nscore": float(interaction.get("nscore", 0.0) or 0.0),
                "fscore": float(interaction.get("fscore", 0.0) or 0.0),
                "pscore": float(interaction.get("pscore", 0.0) or 0.0),
                "ascore": float(interaction.get("ascore", 0.0) or 0.0),
                "escore": float(interaction.get("escore", 0.0) or 0.0),
                "dscore": float(interaction.get("dscore", 0.0) or 0.0),
                "tscore": float(interaction.get("tscore", 0.0) or 0.0),
            }

            if best_score is None or combined > best_score:
                best_score = combined
                best_payload = payload

        return best_payload
    return None


# ─── Path file parser ──────────────────────────────────────────────────────────
def parse_pathway_file(filepath: str):
    """
    Parse pathway files in multiple formats:
    - "Path #N: A -> B -> C"
    - "Unique Path #N (original #M): A -> B -> C"
    """
    paths = []
    seen = set()
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            # Original pathway formats
            m = re.search(r'(?:Path|Unique Path)\s*#\d+.*?:\s*(.+)', line)
            if m:
                seq = m.group(1).strip()
                nodes = [n.strip() for n in seq.split("->")]
                nodes = [n for n in nodes if n]
                # Exclude single-node paths
                if len(nodes) >= 2:
                    key = tuple(nodes)
                    if key not in seen:
                        seen.add(key)
                        paths.append(nodes)
    return paths


# ─── Validation engine ────────────────────────────────────────────────────────────────
class ValidationEngine:
    def __init__(self, log_fn, cancer_category: str = "", use_bulk_cache: bool = True):
        """
        Args:
            log_fn: logging function
            cancer_category: HPA cancer cell line category (e.g. "Brain cancer").
                             Empty string disables category-based fallback.
            use_bulk_cache:  when True (default), the engine consults the local
                             SQLite mirror of HPA's bulk TSV files to fill in
                             cancer-group / cell-line nTPM for any gene whose
                             single-entry JSON returned an empty dict. This is
                             how v16 resolves the "Low cancer specificity" gap.
        """
        self.log = log_fn
        self.cancer_category = (cancer_category or "").strip()
        self.ensembl_cache = {}
        self.hpa_cache = {}
        self.string_cache = {}
        self.use_bulk_cache = use_bulk_cache
        self.bulk = get_bulk_cache(log_fn=log_fn) if use_bulk_cache else None

    @staticmethod
    def _level_rank(level) -> int:
        return {"High": 3, "Medium": 2, "Low": 1, "Not detected": 0}.get(level, 0)

    def _lookup_single_cell(self, hpa_data: dict, cell_type: str):
        """Try to resolve expression from the single-cell-type table.
        Returns (level, source_str, found_bool). source includes [SC] tag."""
        sc = hpa_data.get("single_cell_expression", {})
        ntpm = hpa_data.get("single_cell_ntpm", {})
        ct_lower = cell_type.lower()

        def ntpm_str(key):
            v = ntpm.get(key)
            return f" (nTPM={v:.1f})" if v is not None else ""

        # Step 1: exact key match
        if ct_lower in sc:
            return sc[ct_lower], f"[SC] exact '{ct_lower}'{ntpm_str(ct_lower)}", True

        # Step 2: partial word-set matching
        ct_words = set(ct_lower.replace("(", "").replace(")", "").split())
        best_key, best_ov = None, 0
        for k in sc:
            kw = set(k.replace("(", "").replace(")", "").split())
            ov = len(ct_words & kw)
            if ov > best_ov and ov >= 1:
                best_ov, best_key = ov, k
        if best_key:
            return sc[best_key], f"[SC] partial '{best_key}'{ntpm_str(best_key)}", True

        return None, "", False

    def _lookup_individual_cell_line(self, hpa_data: dict, cell_type: str):
        """Try to resolve expression from the individual-cell-line table.

        v18: searches not only the user's exact cell line but also the entire
        cell-line FAMILY (e.g. HEK293T → HEK293 → 293T → 293). Among ALL hits
        across the family, returns the one with the HIGHEST nTPM — true
        OR-search per the user's request.
        """
        cl = hpa_data.get("cell_line_expression", {})
        cl_ntpm = hpa_data.get("cell_line_ntpm", {})
        if not cl:
            return None, "", False

        family = get_cell_line_family(cell_type)
        if not family:
            return None, "", False

        # Build a normalized index of cache keys for fast lookup.
        # Multiple raw keys may normalize to the same alias (e.g. "HEK 293",
        # "HEK-293", "hek293" all → "hek293"); collect ALL of them and let
        # the OR-best pick the highest-expressing variant.
        cache_index: dict = {}  # alias_norm → [(raw_key, ntpm_or_None), ...]
        for raw_key in cl:
            alias = _norm_cell_line(raw_key)
            cache_index.setdefault(alias, []).append(
                (raw_key, cl_ntpm.get(raw_key))
            )

        # Walk the family in order, collecting every match across every alias.
        # The OR-best is the one with the highest numeric nTPM.
        candidates = []  # list of (raw_key, level, ntpm, family_position)
        for pos, alias in enumerate(family):
            for raw_key, ntpm in cache_index.get(alias, []):
                level = cl.get(raw_key)
                if level is None:
                    continue
                candidates.append((raw_key, level, ntpm, pos))

        if not candidates:
            return None, "", False

        # OR-best: pick max nTPM (None nTPM treated as -1 to lose to any real value).
        # On ties, prefer earlier family position (closer to user's exact request).
        def _sort_key(c):
            raw_key, level, ntpm, pos = c
            return (ntpm if isinstance(ntpm, (int, float)) else -1.0, -pos)
        best = max(candidates, key=_sort_key)
        best_key, best_level, best_ntpm, best_pos = best

        ntpm_str = f" (nTPM={best_ntpm:.1f})" if isinstance(best_ntpm, (int, float)) else ""
        # Build a transparent source string. If we used a family fallback (not
        # the user's exact cell line), say so — the user should know.
        user_norm = family[0]
        best_norm = _norm_cell_line(best_key)
        if best_norm == user_norm:
            label = f"[CL] '{best_key}'{ntpm_str}"
        else:
            label = (f"[CL] family '{best_key}'{ntpm_str} "
                     f"(fallback from {cell_type})")
        return best_level, label, True

    def _lookup_cancer_category(self, hpa_data: dict, category: str):
        """Try to resolve expression from the HPA cancer cell line group table."""
        if not category:
            return None, "", False
        cg = hpa_data.get("cancer_group_expression", {})
        cg_ntpm = hpa_data.get("cancer_group_ntpm", {})
        if not cg:
            return None, "", False
        cat_lower = category.lower().strip()

        def ntpm_str(key):
            v = cg_ntpm.get(key)
            return f" (nTPM={v:.1f})" if v is not None else ""

        if cat_lower in cg:
            return cg[cat_lower], f"[CG] '{category}'{ntpm_str(cat_lower)}", True

        # Loose match
        target = cat_lower.replace(" ", "").replace("-", "")
        for k in cg:
            if k.replace(" ", "").replace("-", "") == target:
                return cg[k], f"[CG] '{k}'{ntpm_str(k)}", True

        return None, "", False

    def get_expression_level(self, hpa_data: dict, cell_type: str) -> tuple:
        """Resolve expression level using OR-based fallback strategy (v13).

        Order of attempts:
          1. Single-cell type table (HPA single-cell nTPM)
          2. Individual cell line table (e.g. HeLa specific nTPM)
          3. Cancer cell line CATEGORY table (e.g. Cervical cancer group nTPM)
          4. HPA specificity field (legacy fallback)

        When MULTIPLE sources return a value, the HIGHEST level (most permissive,
        OR semantics) is taken — this is what the user requested: "or 개념으로
        cell line + cancer category 함께 검색".

        Returns: (level, source_description)
        """
        candidates = []

        # 1. Single-cell type
        lv, src, ok = self._lookup_single_cell(hpa_data, cell_type)
        if ok:
            candidates.append((lv, src))

        # 2. Individual cell line
        lv, src, ok = self._lookup_individual_cell_line(hpa_data, cell_type)
        if ok:
            candidates.append((lv, src))

        # 3. Cancer cell line category
        lv, src, ok = self._lookup_cancer_category(hpa_data, self.cancer_category)
        if ok:
            candidates.append((lv, src))

        if candidates:
            # OR semantics — pick the entry with the highest expression rank.
            best = max(candidates, key=lambda x: self._level_rank(x[0]))
            if len(candidates) == 1:
                return best[0], best[1]
            # Multiple sources — show the chosen one + how many alternatives agreed
            other_summaries = [f"{c[1].split(' ')[0]}={c[0]}" for c in candidates if c is not best]
            extra = f" | also: {', '.join(other_summaries)}" if other_summaries else ""
            return best[0], f"OR-best {best[1]}{extra}"

        # 4. Specificity fallback
        spec = hpa_data.get("cell_specificity", "Unknown")
        spec_map = {
            "Cell type enriched": "High",
            "Tissue enriched": "High",
            "Cell type enhanced": "Medium",
            "Tissue enhanced": "Medium",
            "Group enriched": "Medium",
            "Low cell type specificity": "Low",
            "Low tissue specificity": "Low",
            "Not detected": "Not detected",
        }
        fb_level = "Not detected"
        for sk, lv in spec_map.items():
            if sk.lower() in str(spec).lower():
                fb_level = lv
                break

        sc = hpa_data.get("single_cell_expression", {})
        avail_sample = ", ".join(list(sc.keys())[:3]) if sc else "none"
        if sc:
            return fb_level, f"⚠️ No SC/CL/CG match (e.g. {avail_sample}) → specificity fallback({spec})"
        else:
            return fb_level, f"⚠️ No HPA expression data → specificity fallback({spec})"

    def check_colocalization(self, hpa_a: dict, hpa_b: dict) -> tuple[bool, str]:
        """Check subcellular co-localization between two proteins."""
        locs_a = set(get_compartment_group(l) for l in hpa_a.get("subcellular_locations", []))
        locs_b = set(get_compartment_group(l) for l in hpa_b.get("subcellular_locations", []))
        if not locs_a or not locs_b:
            return None, "Unknown (no localization data)"
        overlap = locs_a & locs_b
        if overlap:
            return True, ", ".join(overlap)
        return False, f"A:{','.join(locs_a)} vs B:{','.join(locs_b)}"

    def get_string_score(self, gene_a: str, gene_b: str) -> dict | None:
        key = tuple(sorted([gene_a, gene_b]))
        if key in self.string_cache:
            return self.string_cache[key]
        score_info = fetch_string_scores(gene_a, gene_b)
        self.string_cache[key] = score_info
        return score_info

    def score_ppi(self, expr_a, expr_b, coloc: bool | None,
                  tau_a: float | None, tau_b: float | None,
                  string_score: dict | float | None, cell_type: str) -> tuple[str, int, list]:
        """
        Calculate the overall validation score
        Returns: (grade, score, reason list)
        """
        score = 0
        reasons = []

        # 1. Expression score (both sides must be expressed)
        # If a tuple is provided (level, source), extract only level
        if isinstance(expr_a, tuple): expr_a = expr_a[0]
        if isinstance(expr_b, tuple): expr_b = expr_b[0]
        level_map = {"High": 3, "Medium": 2, "Low": 1, "Not detected": 0}
        lv_a = level_map.get(expr_a, 0)
        lv_b = level_map.get(expr_b, 0)

        if lv_a == 0 or lv_b == 0:
            reasons.append(f"⚠️ Not expressed: A={expr_a}, B={expr_b}")
        elif lv_a >= 2 and lv_b >= 2:
            score += 3
            reasons.append(f"✅ Both highly expressed: A={expr_a}, B={expr_b}")
        else:
            score += 1
            reasons.append(f"🔶 Low expression: A={expr_a}, B={expr_b}")

        # 2. Co-localization score
        if coloc is True:
            score += 2
            reasons.append("✅ Subcellular co-localization confirmed")
        elif coloc is False:
            reasons.append("❌ Subcellular location mismatch")
        else:
            reasons.append("❓ No localization data")

        # 3. STRING confidence
        if string_score is not None:
            if isinstance(string_score, dict):
                ss = float(string_score.get("score", 0.0) or 0.0)
                esc = float(string_score.get("escore", 0.0) or 0.0)
                dsc = float(string_score.get("dscore", 0.0) or 0.0)
                tsc = float(string_score.get("tscore", 0.0) or 0.0)
                asc = float(string_score.get("ascore", 0.0) or 0.0)
            else:
                ss = float(string_score)
                esc = dsc = tsc = asc = 0.0

            if ss >= 0.9:
                score += 3
                reasons.append(f"✅ High STRING confidence: {ss:.3f}")
            elif ss >= 0.7:
                score += 2
                reasons.append(f"🔶 Medium STRING confidence: {ss:.3f}")
            elif ss >= 0.4:
                score += 1
                reasons.append(f"🔶 Low STRING confidence: {ss:.3f}")
            else:
                reasons.append(f"⚠️ Very low STRING confidence: {ss:.3f}")

            channel_pairs = [("database", dsc), ("experimental", esc), ("textmining", tsc), ("coexpression", asc)]
            best_label, best_value = max(channel_pairs, key=lambda x: x[1])
            if best_value > 0:
                reasons.append(f"ℹ️ Top STRING evidence: {best_label}={best_value:.3f}")
        else:
            reasons.append("❓ STRING No data")

        # 4. Tau score (cell-type specificity)
        if tau_a is not None and tau_b is not None:
            avg_tau = (tau_a + tau_b) / 2
            if avg_tau >= 0.7:
                score += 2
                reasons.append(f"✅ High cell-type specificity Tau={avg_tau:.2f}")
            elif avg_tau >= 0.4:
                score += 1
                reasons.append(f"🔶 Medium cell-type specificity Tau={avg_tau:.2f}")
            else:
                reasons.append(f"🔶 Low cell-type specificity Tau={avg_tau:.2f}")

        # Grade assignment (maximum 10 points)
        if lv_a == 0 or lv_b == 0:
            grade = "NOT_EXPRESSED"
        elif score >= 7:
            grade = "HIGH"
        elif score >= 4:
            grade = "MEDIUM"
        elif score >= 2:
            grade = "LOW"
        else:
            grade = "UNLIKELY"

        return grade, score, reasons

    def validate_paths(self, paths: list, cell_type: str,
                       progress_fn=None, stop_flag=None) -> list:
        """Validate the full list of pathways."""
        results = []
        total_edges = sum(len(p) - 1 for p in paths)
        done = 0

        for path_idx, path in enumerate(paths):
            if stop_flag and stop_flag():
                break

            path_result = {
                "path_idx": path_idx + 1,
                "path": path,
                "edges": [],
                "proteins": {},
                "path_grade": "UNKNOWN",
                "path_score": 0,
            }

            # Collect HPA data for each protein
            for gene in path:
                if gene not in path_result["proteins"]:
                    self.log(f"🔍 [{path_idx + 1}/{len(paths)}] HPA query: {gene}")
                    ensg = gene_symbol_to_ensembl(gene, self.ensembl_cache)
                    if ensg:
                        hpa = fetch_hpa_data(ensg, self.hpa_cache, symbol=gene,
                                              log_fn=self.log)

                        # ── v16/v19: supplement with bulk-TSV cache ────────────────
                        # v16 logic: single-entry JSON returns null for "Low cancer
                        # specificity" genes; bulk TSV always has values, so we
                        # consult bulk cache and merge in any gaps.
                        # v19 fix: ALSO look up by gene name and OR-best the two
                        # paths (ENSG-keyed and name-keyed). This protects against
                        # any remaining wrong-ENSG entries — even if our ENSG
                        # resolves to a foreign gene's data, name-keyed lookup
                        # finds the correct rows in the bulk cache and the
                        # higher-nTPM value wins.
                        if self.bulk is not None and self.bulk.is_populated():
                            ensg_clean = hpa.get("ensg") or (ensg if not ensg.startswith("SYM:") else "")
                            # Two independent lookups — by ENSG and by gene name
                            cg_by_ensg = self.bulk.get_cancer_groups(ensg=ensg_clean) if ensg_clean else {}
                            cg_by_name = self.bulk.get_cancer_groups(gene=gene)
                            cl_by_ensg = self.bulk.get_cell_lines(ensg=ensg_clean) if ensg_clean else {}
                            cl_by_name = self.bulk.get_cell_lines(gene=gene)

                            # OR-best: across both paths, pick the highest nTPM per key.
                            def _or_best_merge(by_ensg: dict, by_name: dict) -> dict:
                                out = {}
                                for k, v in by_ensg.items():
                                    out[k] = v
                                for k, v in by_name.items():
                                    if k not in out or v > out[k]:
                                        out[k] = v
                                return out

                            bulk_cg = _or_best_merge(cg_by_ensg, cg_by_name)
                            bulk_cl = _or_best_merge(cl_by_ensg, cl_by_name)

                            # If ENSG-keyed and name-keyed disagree substantially, the
                            # ENSG was probably wrong. Surface this as a warning so
                            # the user can spot data-integrity problems.
                            if ensg_clean and cg_by_ensg and cg_by_name:
                                ensg_keys = set(cg_by_ensg.keys())
                                name_keys = set(cg_by_name.keys())
                                # Big mismatch = the two ENSGs are likely different genes
                                if ensg_keys and name_keys and len(ensg_keys & name_keys) < min(len(ensg_keys), len(name_keys)) * 0.5:
                                    self.log(f"  ⚠️ [{gene}] ENSG-vs-gene-name lookups disagree strongly — "
                                              f"ENSG '{ensg_clean}' may point to a different gene. "
                                              f"Using OR-best across both paths.")

                            # Merge: bulk fills in keys missing from JSON; if both
                            # have a value for the same key, JSON wins UNLESS the
                            # bulk-OR-best value is markedly higher (>2× the JSON
                            # value), which indicates the JSON came from a wrong
                            # ENSG match. In that case bulk wins.
                            existing_cg_n = hpa.get("cancer_group_ntpm", {})
                            existing_cg_e = hpa.get("cancer_group_expression", {})
                            for cat, ntpm in bulk_cg.items():
                                prior = existing_cg_n.get(cat)
                                if prior is None:
                                    existing_cg_n[cat] = ntpm
                                    existing_cg_e[cat] = ntpm_to_level(ntpm)
                                elif isinstance(prior, (int, float)) and ntpm > 2.0 * prior + 1.0:
                                    # Bulk thinks much higher — JSON likely wrong-gene
                                    existing_cg_n[cat] = ntpm
                                    existing_cg_e[cat] = ntpm_to_level(ntpm)
                            hpa["cancer_group_ntpm"] = existing_cg_n
                            hpa["cancer_group_expression"] = existing_cg_e

                            existing_cl_n = hpa.get("cell_line_ntpm", {})
                            existing_cl_e = hpa.get("cell_line_expression", {})
                            for cl, ntpm in bulk_cl.items():
                                prior = existing_cl_n.get(cl)
                                if prior is None:
                                    existing_cl_n[cl] = ntpm
                                    existing_cl_e[cl] = ntpm_to_level(ntpm)
                                elif isinstance(prior, (int, float)) and ntpm > 2.0 * prior + 1.0:
                                    existing_cl_n[cl] = ntpm
                                    existing_cl_e[cl] = ntpm_to_level(ntpm)
                            hpa["cell_line_ntpm"] = existing_cl_n
                            hpa["cell_line_expression"] = existing_cl_e

                            if bulk_cg or bulk_cl:
                                self.log(f"  💾 [{gene}] bulk cache: "
                                          f"{len(bulk_cg)} cancer groups, {len(bulk_cl)} cell lines "
                                          f"(ENSG-path={len(cg_by_ensg)}, name-path={len(cg_by_name)})")

                        # ── v14: detailed parsed-data diagnostics
                        sc_n = len(hpa.get("single_cell_expression", {}))
                        cl_n = len(hpa.get("cell_line_expression", {}))
                        cg_n = len(hpa.get("cancer_group_expression", {}))
                        self.log(
                            f"  [HPA] {'✅' if hpa.get('raw') else '⚠️'} {gene}: "
                            f"SC-spec={hpa.get('cell_specificity', '?')}, "
                            f"CL-spec={hpa.get('cell_line_specificity', '?')} | "
                            f"keys parsed: SC={sc_n}, CL={cl_n}, CG={cg_n}")
                        # If a cancer category is selected, show its actual nTPM value:
                        if self.cancer_category and cg_n > 0:
                            cat_l = self.cancer_category.lower()
                            cg_lvl = hpa.get("cancer_group_expression", {}).get(cat_l, "—")
                            cg_v = hpa.get("cancer_group_ntpm", {}).get(cat_l)
                            v_str = f"{cg_v:.1f}" if isinstance(cg_v, (int, float)) else "?"
                            self.log(f"     ↳ {self.cancer_category}: level={cg_lvl}, nTPM={v_str}")
                    else:
                        hpa = {"ensg": None, "name": gene,
                               "tau_score": None, "cell_specificity": "Unknown",
                               "subcellular_locations": [], "single_cell_expression": {},
                               "single_cell_ntpm": {},
                               "cell_line_expression": {}, "cell_line_ntpm": {},
                               "cancer_group_expression": {}, "cancer_group_ntpm": {},
                               "cell_line_specificity": "Unknown", "cell_line_tau": None}
                    path_result["proteins"][gene] = hpa
                    time.sleep(0.3)  # API rate limit

            # Edge-by-edge validation
            edge_grades = []
            for i in range(len(path) - 1):
                gene_a, gene_b = path[i], path[i + 1]
                hpa_a = path_result["proteins"][gene_a]
                hpa_b = path_result["proteins"][gene_b]

                expr_a, expr_src_a = self.get_expression_level(hpa_a, cell_type)
                expr_b, expr_src_b = self.get_expression_level(hpa_b, cell_type)
                # ── v14: log the OR-resolved expression for both proteins
                self.log(f"  ↳ {gene_a}: {expr_a}  ←  {expr_src_a}")
                self.log(f"  ↳ {gene_b}: {expr_b}  ←  {expr_src_b}")
                coloc, coloc_detail = self.check_colocalization(hpa_a, hpa_b)

                self.log(f"  → STRING: {gene_a} — {gene_b}")
                string_score = self.get_string_score(gene_a, gene_b)
                if isinstance(string_score, dict):
                    self.log(
                        "    [STRING] "
                        f"combined={float(string_score.get('score', 0.0) or 0.0):.3f} | "
                        f"exp={float(string_score.get('escore', 0.0) or 0.0):.3f}, "
                        f"db={float(string_score.get('dscore', 0.0) or 0.0):.3f}, "
                        f"text={float(string_score.get('tscore', 0.0) or 0.0):.3f}, "
                        f"coexp={float(string_score.get('ascore', 0.0) or 0.0):.3f}"
                    )
                time.sleep(0.2)

                grade, score, reasons = self.score_ppi(
                    expr_a, expr_b, coloc,
                    hpa_a.get("tau_score"), hpa_b.get("tau_score"),
                    string_score, cell_type
                )

                # v13: pull cancer-group level for both proteins (for the diagnostics column)
                cg_a_lvl = cg_b_lvl = None
                if self.cancer_category:
                    cat_l = self.cancer_category.lower()
                    cg_a_lvl = hpa_a.get("cancer_group_expression", {}).get(cat_l)
                    cg_b_lvl = hpa_b.get("cancer_group_expression", {}).get(cat_l)

                edge_result = {
                    "gene_a": gene_a, "gene_b": gene_b,
                    "expr_a": expr_a, "expr_b": expr_b,
                    "expr_src_a": expr_src_a, "expr_src_b": expr_src_b,
                    "coloc": coloc, "coloc_detail": coloc_detail,
                    "string_score": (string_score.get("score") if isinstance(string_score, dict) else string_score),
                    "string_channels": (string_score if isinstance(string_score, dict) else None),
                    "grade": grade, "score": score, "reasons": reasons,
                    "locs_a": hpa_a.get("subcellular_locations", []),
                    "locs_b": hpa_b.get("subcellular_locations", []),
                    "tau_a": hpa_a.get("tau_score"),
                    "tau_b": hpa_b.get("tau_score"),
                    "sc_keys_a": list(hpa_a.get("single_cell_expression", {}).keys())[:5],
                    "sc_keys_b": list(hpa_b.get("single_cell_expression", {}).keys())[:5],
                    "cg_level_a": cg_a_lvl,
                    "cg_level_b": cg_b_lvl,
                }
                path_result["edges"].append(edge_result)
                edge_grades.append(grade)

                done += 1
                if progress_fn:
                    progress_fn(done, total_edges)

            # Overall pathway grade — take the WORST edge.
            # v17 fix: previous min() with grade_order.index() returned the BEST
            # edge grade, not the worst, because HIGH has index 0 (smallest).
            # We want the highest index (= worst) to dominate the path grade.
            grade_order = ["HIGH", "MEDIUM", "LOW", "UNLIKELY", "NOT_EXPRESSED", "UNKNOWN"]
            worst = max(edge_grades, key=lambda g: grade_order.index(g) if g in grade_order else -1,
                        default="UNKNOWN")
            path_result["path_grade"] = worst
            path_result["path_score"] = sum(e["score"] for e in path_result["edges"])
            results.append(path_result)

        return results


# ─── HTML generator ──────────────────────────────────────────────────────────────
def grade_color(grade: str) -> str:
    return {
        "HIGH": "#22c55e",
        "MEDIUM": "#f59e0b",
        "LOW": "#f97316",
        "UNLIKELY": "#ef4444",
        "NOT_EXPRESSED": "#991b1b",
        "UNKNOWN": "#94a3b8",
    }.get(grade, "#94a3b8")


def grade_bg(grade: str) -> str:
    return {
        "HIGH": "#dcfce7",
        "MEDIUM": "#fef3c7",
        "LOW": "#ffedd5",
        "UNLIKELY": "#fee2e2",
        "NOT_EXPRESSED": "#fecaca",
        "UNKNOWN": "#f1f5f9",
    }.get(grade, "#f1f5f9")


def grade_label_kr(grade: str) -> str:
    return {
        "HIGH": "✅ High",
        "MEDIUM": "🔶 Medium",
        "LOW": "⚠️ Low",
        "UNLIKELY": "❌ Unlikely",
        "NOT_EXPRESSED": "🚫 Not expressed",
        "UNKNOWN": "❓ Unknown",
    }.get(grade, grade)


def expr_badge(level) -> str:
    if isinstance(level, tuple): level = level[0]  # tuple guard
    colors = {"High": "#166534", "Medium": "#92400e", "Low": "#1e40af", "Not detected": "#6b7280"}
    bgs = {"High": "#dcfce7", "Medium": "#fef3c7", "Low": "#dbeafe", "Not detected": "#f3f4f6"}
    c = colors.get(level, "#6b7280")
    b = bgs.get(level, "#f3f4f6")
    return f'<span style="background:{b};color:{c};padding:2px 8px;border-radius:12px;font-size:12px;font-weight:600">{level}</span>'


def generate_html(results: list, cell_type: str, filepath: str,
                  source_file: str = "", cancer_category: str = "") -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Summary statistics
    grade_counts = defaultdict(int)
    for r in results:
        grade_counts[r["path_grade"]] += 1
    total = len(results)
    high          = grade_counts["HIGH"]
    medium        = grade_counts["MEDIUM"]
    low           = grade_counts["LOW"]
    unlikely      = grade_counts["UNLIKELY"]
    not_expressed = grade_counts["NOT_EXPRESSED"]
    unknown       = grade_counts["UNKNOWN"]

    # Collect all proteins
    all_proteins = {}
    for r in results:
        for g, hpa in r["proteins"].items():
            if g not in all_proteins:
                all_proteins[g] = hpa

    # ── Build pathway card HTML ─────────────────────────────────────────
    path_cards = []
    for r in results:
        grade = r["path_grade"]
        path_str = " → ".join(r["path"])
        edge_rows = []
        for e in r["edges"]:
            ss = f"{float(e['string_score']):.3f}" if e['string_score'] is not None else "N/A"
            channels = e.get("string_channels") or {}
            string_sub = (
                f"Exp {float(channels.get('escore', 0.0) or 0.0):.3f} | "
                f"DB {float(channels.get('dscore', 0.0) or 0.0):.3f} | "
                f"Text {float(channels.get('tscore', 0.0) or 0.0):.3f} | "
                f"Coexp {float(channels.get('ascore', 0.0) or 0.0):.3f}"
            ) if channels else "No channel data"
            locs_a = ", ".join(e["locs_a"][:2]) or "Unknown"
            locs_b = ", ".join(e["locs_b"][:2]) or "Unknown"
            tau_a = f"{e['tau_a']:.2f}" if e['tau_a'] else "N/A"
            tau_b = f"{e['tau_b']:.2f}" if e['tau_b'] else "N/A"
            coloc_icon = "✅" if e["coloc"] else ("❌" if e["coloc"] is False else "❓")
            reasons_html = "<br>".join(e["reasons"])
            src_a = e.get("expr_src_a", "")
            src_b = e.get("expr_src_b", "")
            sc_a = ", ".join(e.get("sc_keys_a", [])[:3]) or "none"
            sc_b = ", ".join(e.get("sc_keys_b", [])[:3]) or "none"
            edge_rows.append(f"""
            <tr style="border-bottom:1px solid #e2e8f0">
              <td style="padding:8px;font-weight:600">{e['gene_a']}</td>
              <td style="padding:8px;text-align:center">→</td>
              <td style="padding:8px;font-weight:600">{e['gene_b']}</td>
              <td style="padding:8px;text-align:center">{expr_badge(e['expr_a'])}<br><small style="color:#94a3b8;font-size:10px;display:block;max-width:130px;overflow:hidden">{src_a[:40]}</small></td>
              <td style="padding:8px;text-align:center">{expr_badge(e['expr_b'])}<br><small style="color:#94a3b8;font-size:10px;display:block;max-width:130px;overflow:hidden">{src_b[:40]}</small></td>
              <td style="padding:8px;text-align:center;font-size:12px">{locs_a}<br><span style="color:#64748b">{tau_a}</span></td>
              <td style="padding:8px;text-align:center;font-size:12px">{locs_b}<br><span style="color:#64748b">{tau_b}</span></td>
              <td style="padding:8px;text-align:center">{coloc_icon} {e.get('coloc_detail', '')[:30]}</td>
              <td style="padding:8px;text-align:center;font-weight:600">{ss}<br><span style="color:#64748b;font-size:10px">{string_sub}</span></td>
              <td style="padding:8px;background:{grade_bg(e['grade'])};color:{grade_color(e['grade'])};font-weight:700;text-align:center">{grade_label_kr(e['grade'])}</td>
              <td style="padding:8px;font-size:11px;color:#475569">{reasons_html}</td>
            </tr>""")

        path_cards.append(f"""
        <div class="path-card" style="margin:16px 0;border:2px solid {grade_color(grade)};border-radius:12px;overflow:hidden">
          <div style="background:{grade_color(grade)};padding:10px 16px;display:flex;justify-content:space-between;align-items:center">
            <span style="color:white;font-weight:700;font-size:14px">Path #{r['path_idx']}: {path_str}</span>
            <span style="background:white;color:{grade_color(grade)};padding:3px 12px;border-radius:16px;font-weight:800;font-size:13px">{grade_label_kr(grade)} (Score: {r['path_score']})</span>
          </div>
          <div style="overflow-x:auto">
            <table style="width:100%;border-collapse:collapse;font-size:13px">
              <thead>
                <tr style="background:#f8fafc;font-size:12px;color:#64748b">
                  <th style="padding:8px;text-align:left">Protein A</th>
                  <th style="padding:8px"></th>
                  <th style="padding:8px;text-align:left">Protein B</th>
                  <th style="padding:8px;text-align:center">Expression (A)</th>
                  <th style="padding:8px;text-align:center">Expression (B)</th>
                  <th style="padding:8px;text-align:center">Location (A) / Tau</th>
                  <th style="padding:8px;text-align:center">Location (B) / Tau</th>
                  <th style="padding:8px;text-align:center">Co-localization</th>
                  <th style="padding:8px;text-align:center">STRING combined / channel scores</th>
                  <th style="padding:8px;text-align:center">Grade</th>
                  <th style="padding:8px;text-align:left">Evidence</th>
                </tr>
              </thead>
              <tbody>{"".join(edge_rows)}</tbody>
            </table>
          </div>
        </div>""")

    # ── Protein summary table ───────────────────────────────────────────
    prot_rows = []
    for gene, hpa in sorted(all_proteins.items()):
        ensg = hpa.get("ensg", "N/A") or "N/A"
        tau = f"{hpa['tau_score']:.2f}" if hpa.get("tau_score") else "N/A"
        spec = hpa.get("cell_specificity", "Unknown")
        locs = ", ".join(hpa.get("subcellular_locations", [])[:3]) or "N/A"
        hpa_link = f'<a href="https://www.proteinatlas.org/{ensg}" target="_blank" style="color:#3b82f6">{ensg}</a>' if ensg != "N/A" else "N/A"
        prot_rows.append(f"""
        <tr style="border-bottom:1px solid #e2e8f0">
          <td style="padding:8px;font-weight:700">{gene}</td>
          <td style="padding:8px">{hpa_link}</td>
          <td style="padding:8px">{tau}</td>
          <td style="padding:8px;font-size:12px">{spec}</td>
          <td style="padding:8px;font-size:12px">{locs}</td>
        </tr>""")

    # ── Bar chart data
    # ─ Expression matching diagnostics table (v13 OR-search)
    seen_genes_debug = set()
    debug_rows_list = []
    cat_lower = (cancer_category or "").lower()
    # Strip trailing " (cancer cell line)" etc. to match HPA cell line keys
    ct_clean = cell_type.lower()
    for sfx in [" (cancer cell line)", " cell line"]:
        if ct_clean.endswith(sfx):
            ct_clean = ct_clean[: -len(sfx)].strip()

    bg_map = {"High": "#dcfce7", "Medium": "#fef3c7", "Low": "#dbeafe",
              "Not detected": "#fce4ec"}
    rank_map = {"High": 3, "Medium": 2, "Low": 1, "Not detected": 0}

    # v20: One unified lookup engine for the diagnostics table.
    # We disable bulk-cache here (data already merged into hpa dicts during
    # validation) and pass a no-op log function. The engine just provides
    # the matching methods — _lookup_single_cell, _lookup_individual_cell_line,
    # _lookup_cancer_category — so the HTML report uses identical logic to
    # what the validator used.
    _summary_engine = ValidationEngine(log_fn=lambda m: None,
                                        cancer_category=cancer_category,
                                        use_bulk_cache=False)

    for r in results:
        for gene, hpa in r["proteins"].items():
            if gene in seen_genes_debug: continue
            seen_genes_debug.add(gene)
            sc = hpa.get("single_cell_expression", {})
            sc_n = hpa.get("single_cell_ntpm", {})
            cl = hpa.get("cell_line_expression", {})
            cl_n = hpa.get("cell_line_ntpm", {})
            cg = hpa.get("cancer_group_expression", {})
            cg_n = hpa.get("cancer_group_ntpm", {})

            # ── v20: Unified lookup via ValidationEngine ──────────────────
            # Previously this block had its own inline lookup code that drifted
            # from ValidationEngine's logic. Now we delegate to the engine's
            # methods so the HTML diagnostics ALWAYS match what the validator
            # actually used (including HEK293T → HEK293 family fallback).
            sc_lvl_v, sc_src_v, sc_ok = _summary_engine._lookup_single_cell(hpa, cell_type)
            cl_lvl_v, cl_src_v, cl_ok = _summary_engine._lookup_individual_cell_line(hpa, cell_type)
            cg_lvl_v, cg_src_v, cg_ok = _summary_engine._lookup_cancer_category(hpa, cancer_category)

            sc_lvl  = sc_lvl_v if sc_ok else "—"
            sc_label = sc_src_v if sc_ok else "no data"
            cl_lvl  = cl_lvl_v if cl_ok else "—"
            cl_label = cl_src_v if cl_ok else "no data"
            cg_lvl  = cg_lvl_v if cg_ok else "—"
            cg_label = cg_src_v if cg_ok else (
                "no data" if cancer_category else "no category selected"
            )

            # — final OR-best level
            ranked = [(sc_lvl, "SC"), (cl_lvl, "CL"), (cg_lvl, "CG")]
            ranked = [(lv, src) for lv, src in ranked if lv != "—"]
            if ranked:
                best_lv, best_src = max(ranked, key=lambda x: rank_map.get(x[0], 0))
                if len(ranked) > 1:
                    final_label = f"OR-best ({best_src}): {best_lv}"
                else:
                    final_label = f"{best_src}: {best_lv}"
            else:
                spec = hpa.get("cell_specificity", "Unknown")
                spec_map_local = {"Cell type enriched": "High", "Tissue enriched": "High",
                                  "Cell type enhanced": "Medium", "Tissue enhanced": "Medium",
                                  "Group enriched": "Medium", "Low cell type specificity": "Low",
                                  "Low tissue specificity": "Low", "Not detected": "Not detected"}
                best_lv = "Not detected"
                for sk, lv in spec_map_local.items():
                    if sk.lower() in str(spec).lower():
                        best_lv = lv; break
                final_label = f"⚠️ specificity fallback ({spec})"

            bg = bg_map.get(best_lv, "#f1f5f9")

            def cell_html(lvl, label):
                bg_l = bg_map.get(lvl, "#f1f5f9") if lvl != "—" else "transparent"
                return (f'<div style="background:{bg_l};padding:2px 6px;border-radius:8px;'
                        f'font-size:11px;display:inline-block;font-weight:600">{lvl}</div>'
                        f'<div style="font-size:10px;color:#64748b;margin-top:2px">{label}</div>')

            avail_sc = ", ".join(list(sc.keys())[:4]) or "—"
            avail_cg = ", ".join(list(cg.keys())[:4]) or "—"

            debug_rows_list.append(
                f'<tr style="border-bottom:1px solid #e2e8f0">'
                f'<td style="padding:8px;font-weight:700">{gene}</td>'
                f'<td style="padding:8px"><span style="background:{bg};padding:3px 10px;border-radius:10px;font-size:12px;font-weight:700">{best_lv}</span>'
                f'<div style="font-size:10px;color:#64748b;margin-top:3px">{final_label}</div></td>'
                f'<td style="padding:8px">{cell_html(sc_lvl, sc_label)}</td>'
                f'<td style="padding:8px">{cell_html(cl_lvl, cl_label)}</td>'
                f'<td style="padding:8px">{cell_html(cg_lvl, cg_label)}</td>'
                f'<td style="padding:8px;font-size:10px;color:#94a3b8">'
                f'<div><b>SC:</b> {avail_sc}</div>'
                f'<div><b>CG:</b> {avail_cg}</div></td></tr>')
    debug_rows = "".join(
        debug_rows_list) if debug_rows_list else "<tr><td colspan=6 style='padding:12px;color:#94a3b8'>No diagnostic data</td></tr>"

    grade_labels = ["HIGH", "MEDIUM", "LOW", "UNLIKELY", "NOT_EXPRESSED", "UNKNOWN"]
    grade_data = [grade_counts[g] for g in grade_labels]
    grade_colors_js = [grade_color(g) for g in grade_labels]
    grade_labels_kr = [grade_label_kr(g) for g in grade_labels]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PPI Cell-Type Feasibility Validation Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4"></script>
<style>
  * {{ box-sizing:border-box; margin:0; padding:0; }}
  body {{ font-family:Arial, sans-serif;background:#f0f4f8;color:#1e293b }}
  .container {{ max-width:1400px;margin:0 auto;padding:24px }}
  .header {{ background:linear-gradient(135deg,#0d2137 0%,#1e40af 100%);color:white;padding:32px;border-radius:16px;margin-bottom:24px }}
  .header h1 {{ font-size:24px;font-weight:800;margin-bottom:8px }}
  .header .sub {{ font-size:14px;opacity:0.8 }}
  .stat-grid {{ display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:16px;margin-bottom:24px }}
  .stat-card {{ background:white;border-radius:12px;padding:20px;text-align:center;box-shadow:0 1px 4px rgba(0,0,0,.08) }}
  .stat-card .num {{ font-size:36px;font-weight:800 }}
  .stat-card .lbl {{ font-size:13px;color:#64748b;margin-top:4px }}
  .section {{ background:white;border-radius:12px;padding:24px;margin-bottom:24px;box-shadow:0 1px 4px rgba(0,0,0,.08) }}
  .section h2 {{ font-size:18px;font-weight:700;margin-bottom:16px;padding-bottom:8px;border-bottom:2px solid #e2e8f0 }}
  .path-card {{ transition:box-shadow .2s }}
  .path-card:hover {{ box-shadow:0 4px 16px rgba(0,0,0,.12) }}
  .filter-bar {{ display:flex;gap:12px;margin-bottom:16px;flex-wrap:wrap }}
  .filter-btn {{ padding:6px 16px;border:2px solid #e2e8f0;border-radius:20px;cursor:pointer;font-size:13px;background:white;transition:all .2s }}
  .filter-btn.active {{ border-color:#3b82f6;background:#dbeafe;color:#1e40af;font-weight:600 }}
  .legend {{ display:flex;gap:16px;flex-wrap:wrap;margin-bottom:12px }}
  .legend-item {{ display:flex;align-items:center;gap:6px;font-size:13px }}
  .legend-dot {{ width:12px;height:12px;border-radius:50% }}
  table {{ width:100%;border-collapse:collapse }}
  th {{ padding:10px 8px;text-align:left;font-size:12px;color:#64748b;font-weight:600 }}
  tr:hover td {{ background:#f8fafc }}
  @media print {{ .filter-bar {{ display:none }} }}
</style>
</head>
<body>
<div class="container">

<!-- HEADER -->
<div class="header">
  <h1>🔬 PPI Cell-Type Feasibility Validation Report</h1>
  <div class="sub">
    Target cell type: <strong>{cell_type}</strong> &nbsp;|&nbsp;
    Cancer cell line category (OR-search): <strong>{cancer_category or "— (none)"}</strong> &nbsp;|&nbsp;
    Number of analyzed pathways: <strong>{total}</strong> &nbsp;|&nbsp;
    Analysis time: {now}
    {f'<br>Source file: {os.path.basename(source_file)}' if source_file else ''}
  </div>
</div>

<!-- Summary statistics -->
<div class="stat-grid">
  <div class="stat-card"><div class="num" style="color:#22c55e">{high}</div><div class="lbl">✅ High feasibility (HIGH)</div></div>
  <div class="stat-card"><div class="num" style="color:#f59e0b">{medium}</div><div class="lbl">🔶 Medium feasibility (MEDIUM)</div></div>
  <div class="stat-card"><div class="num" style="color:#fb923c">{low}</div><div class="lbl">⚠️ Low feasibility (LOW)</div></div>
  {('<div class="stat-card"><div class="num" style="color:#ef4444">' + str(unlikely) + '</div><div class="lbl">❌ Unlikely (UNLIKELY)</div></div>') if unlikely else ''}
  <div class="stat-card"><div class="num" style="color:#dc2626">{not_expressed}</div><div class="lbl">🚫 Not expressed</div></div>
  {('<div class="stat-card"><div class="num" style="color:#6b7280">' + str(unknown) + '</div><div class="lbl">❓ Unknown</div></div>') if unknown else ''}
  <div class="stat-card"><div class="num">{total}</div><div class="lbl">Total validated pathways</div></div>
  <div class="stat-card"><div class="num">{len(all_proteins)}</div><div class="lbl">Number of unique proteins</div></div>
</div>

<!-- Charts -->
<div class="section">
  <h2>📊 Validation grade distribution</h2>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;max-width:780px;margin:0 auto">
    <div style="position:relative;height:200px"><canvas id="gradeChart"></canvas></div>
    <div style="position:relative;height:200px"><canvas id="pieChart"></canvas></div>
  </div>
</div>

<!-- Validation legend -->
<div class="section">
  <h2>🎯 Validation methodology</h2>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px;font-size:13px">
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>① Expression validation</strong><br>Based on HPA single-cell expression data<br>High/Medium/Low/Not detected</div>
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>② Tau score</strong><br>Cell-type specificity metric (0~1)<br>≥0.7 = high specificity, ≥0.4 = medium specificity</div>
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>③ Subcellular co-localization</strong><br>Comparison of HPA subcellular locations<br>Same organelle group = co-localization possible</div>
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>④ STRING confidence</strong><br>STRING v11.5 combined PPI score<br>≥0.9 = high confidence, ≥0.7 = medium confidence</div>
  </div>
  <div class="legend" style="margin-top:16px">
    {"".join(f'<div class="legend-item"><div class="legend-dot" style="background:{grade_color(g)}"></div><span>{grade_label_kr(g)}</span></div>' for g in grade_labels)}
  </div>
</div>

<!-- Path filters -->
<div class="section">
  <h2>🧬 PPI validation results by pathway</h2>
  <div class="filter-bar">
    <button class="filter-btn active" onclick="filterPaths('ALL')">Show all</button>
    {"".join(f'<button class="filter-btn" onclick="filterPaths(&quot;{g}&quot;)" style="border-color:{grade_color(g)}">{grade_label_kr(g)}</button>' for g in grade_labels if grade_counts[g] > 0)}
  </div>
  <div id="path-container">
    {"".join(f'<div class="path-card grade-{r["path_grade"]}">{card}</div>' for r, card in zip(results, path_cards))}
  </div>
</div>

<!-- Protein summary -->
<div class="section">
  <h2>🔍 Expression matching diagnostics — OR-search across SC + CL + CG</h2>
  <p style="font-size:13px;color:#64748b;margin-bottom:12px">
    <strong>Selected cell type:</strong> {cell_type}
    {f' &nbsp;|&nbsp; <strong>Cancer cell line category:</strong> {cancer_category}' if cancer_category else ''}
    <br>
    <strong>SC</strong> = HPA single-cell type table &nbsp;|&nbsp;
    <strong>CL</strong> = individual HPA cell line (e.g. HeLa, A-549) &nbsp;|&nbsp;
    <strong>CG</strong> = HPA cancer cell line group (e.g. Cervical cancer, Lung cancer)
    <br>
    The <strong>OR-best</strong> column shows the highest expression level found across all
    three sources — this is the value used by the validator. If no source has data,
    a fallback to HPA's specificity field is used.
  </p>
  <div style="overflow-x:auto">
    <table><thead><tr style="background:#f8fafc">
      <th style="padding:8px">Gene</th>
      <th style="padding:8px">OR-best (used)</th>
      <th style="padding:8px">SC — single cell type</th>
      <th style="padding:8px">CL — individual cell line</th>
      <th style="padding:8px">CG — cancer cell line group</th>
      <th style="padding:8px">Available HPA keys</th>
    </tr></thead>
    <tbody>{debug_rows}</tbody></table>
  </div>
</div>

<div class="section">
  <h2>🔬 Protein expression profile summary</h2>
  <div style="overflow-x:auto">
    <table>
      <thead>
        <tr style="background:#f8fafc">
          <th>Gene</th>
          <th>Ensembl ID</th>
          <th>Tau Score</th>
          <th>Cell specificity</th>
          <th>Subcellular location</th>
        </tr>
      </thead>
      <tbody>{"".join(prot_rows)}</tbody>
    </table>
  </div>
</div>

<!-- Footer -->
<div style="text-align:center;color:#94a3b8;font-size:12px;padding:16px 0">
  PPI Cell-Type Validator v{APP_VERSION} &nbsp;|&nbsp; AI4Emotion &nbsp;|&nbsp;
  HPA API · STRING v12 API · MyGene.info based on &nbsp;|&nbsp; Generated: {now}
</div>
</div>

<script>
// ─── Charts
const grades  = {json.dumps(grade_labels_kr)};
const counts  = {json.dumps(grade_data)};
const colors  = {json.dumps(grade_colors_js)};

new Chart(document.getElementById('gradeChart'), {{
  type:'bar',
  data:{{ labels:grades, datasets:[{{ label:'Path count', data:counts,
    backgroundColor:colors, borderRadius:8 }}]}},
  options:{{ responsive:true, maintainAspectRatio:false,
    plugins:{{legend:{{display:false}},
    title:{{display:true,text:'Path count by grade'}}}} }}
}});

new Chart(document.getElementById('pieChart'), {{
  type:'doughnut',
  data:{{ labels:grades, datasets:[{{ data:counts, backgroundColor:colors }}]}},
  options:{{ responsive:true, maintainAspectRatio:false,
    plugins:{{legend:{{position:'right',labels:{{boxWidth:12,font:{{size:11}}}}}},
    title:{{display:true,text:'Proportion by grade'}}}} }}
}});

// ─── Filters
function filterPaths(grade) {{
  document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
  document.querySelectorAll('#path-container .path-card').forEach(card => {{
    card.style.display = (grade === 'ALL' || card.classList.contains('grade-' + grade)) ? 'block' : 'none';
  }});
}}
</script>
</body>
</html>"""

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(html)
    return filepath


# ─── Tkinter GUI ──────────────────────────────────────────────────────────────
class PPIValidatorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"PPI Cell-Type Validator v{APP_VERSION}")
        self.geometry("900x830")
        self.resizable(True, True)
        self.configure(bg="#f0f4f8")

        self._stop = False
        self._running = False
        self._paths = []
        self._source_file = ""
        # v17 fix: track the last auto-suggested category so we can tell whether
        # the current value came from auto-suggest (safe to overwrite) or from
        # an explicit user pick (must be preserved).
        self._last_auto_category = ""

        self._build_ui()

    def _build_ui(self):
        # ── Top header
        hdr = tk.Frame(self, bg="#0d2137", pady=14)
        hdr.pack(fill="x")
        tk.Label(hdr, text="🔬  PPI Cell-Type Validator",
                 font=("Arial", 16, "bold"), bg="#0d2137", fg="white").pack()
        tk.Label(hdr, text="Validate cell-type feasibility of PPIs in BFS pathway data",
                 font=("Arial", 10), bg="#0d2137", fg="#90caf9").pack()

        body = tk.Frame(self, bg="#f0f4f8", padx=20, pady=16)
        body.pack(fill="both", expand=True)

        # ── File selection
        f1 = ttk.LabelFrame(body, text="① Select pathway file", padding=12)
        f1.pack(fill="x", pady=6)
        self.file_var = tk.StringVar(value="Please select a file...")
        tk.Label(f1, textvariable=self.file_var, width=60, anchor="w",
                 relief="sunken", bg="white", font=("Arial", 9)).pack(side="left", padx=4)
        ttk.Button(f1, text="Open", command=self._choose_file).pack(side="left", padx=4)
        self.path_count_lbl = tk.Label(f1, text="", fg="#16a34a", font=("Arial", 9))
        self.path_count_lbl.pack(side="left", padx=8)

        # ── Cell type selection
        f2 = ttk.LabelFrame(body, text="② Select target cell type & cancer cell line category", padding=12)
        f2.pack(fill="x", pady=6)

        # Row 1: cell type
        row1 = tk.Frame(f2)
        row1.pack(fill="x", pady=2)
        tk.Label(row1, text="Cell type:", font=("Arial", 9), width=14, anchor="w").pack(side="left")
        self.cell_var = tk.StringVar(value=CELL_TYPES[0])
        cb = ttk.Combobox(row1, textvariable=self.cell_var, values=CELL_TYPES,
                          state="readonly", width=38)
        cb.pack(side="left", padx=4)
        cb.bind("<<ComboboxSelected>>", self._on_cell_type_changed)
        tk.Label(row1, text="Or enter manually:", font=("Arial", 9)).pack(side="left", padx=(10, 2))
        self.cell_custom = ttk.Entry(row1, width=18)
        self.cell_custom.pack(side="left", padx=4)
        self.cell_custom.bind("<FocusOut>", self._on_cell_type_changed)
        self.cell_custom.bind("<KeyRelease>", self._on_cell_type_changed)

        # Row 2: cancer cell line category (NEW in v13)
        row2 = tk.Frame(f2)
        row2.pack(fill="x", pady=4)
        tk.Label(row2, text="Cancer category:", font=("Arial", 9, "bold"), width=14, anchor="w",
                 fg="#1e40af").pack(side="left")
        self.category_var = tk.StringVar(value=CANCER_CELL_LINE_CATEGORIES[0])
        self.category_cb = ttk.Combobox(row2, textvariable=self.category_var,
                                         values=CANCER_CELL_LINE_CATEGORIES,
                                         state="readonly", width=38)
        self.category_cb.pack(side="left", padx=4)
        # v17 fix: record explicit user picks so auto-suggest doesn't clobber them.
        self.category_cb.bind("<<ComboboxSelected>>", self._on_category_user_pick)
        tk.Label(row2,
                 text="(OR-search: cell line + category → fills ‘No data’ gaps)",
                 font=("Arial", 8), fg="#64748b").pack(side="left", padx=8)

        # Try to suggest a category for the default cell type
        self._on_cell_type_changed()

        # ── Parameters
        f3 = ttk.LabelFrame(body, text="③ Analysis parameters", padding=12)
        f3.pack(fill="x", pady=6)
        tk.Label(f3, text="Maximum path count:", font=("Arial", 9)).grid(row=0, column=0, sticky="w")
        self.max_paths = ttk.Spinbox(f3, from_=1, to=999, width=6)
        self.max_paths.set(50)
        self.max_paths.grid(row=0, column=1, padx=8, sticky="w")
        tk.Label(f3, text="Minimum STRING score:", font=("Arial", 9)).grid(row=0, column=2, padx=16, sticky="w")
        self.min_score = ttk.Spinbox(f3, from_=0.0, to=1.0, increment=0.1, width=6)
        self.min_score.set(0.4)
        self.min_score.grid(row=0, column=3, padx=8, sticky="w")

        # ── Output file
        f4 = ttk.LabelFrame(body, text="④ HTML output file", padding=12)
        f4.pack(fill="x", pady=6)
        self.out_var = tk.StringVar(value="ppi_validation_result.html")
        tk.Label(f4, textvariable=self.out_var, width=55, anchor="w",
                 relief="sunken", bg="white", font=("Arial", 9)).pack(side="left", padx=4)
        ttk.Button(f4, text="Choose location", command=self._choose_output).pack(side="left", padx=4)

        # ── Progress bar
        self.progress = ttk.Progressbar(body, length=800, mode="determinate")
        self.progress.pack(fill="x", pady=8)
        self.status_lbl = tk.Label(body, text="Ready", font=("Arial", 9),
                                   fg="#475569", bg="#f0f4f8")
        self.status_lbl.pack()

        # ── Log
        f5 = ttk.LabelFrame(body, text="Execution log", padding=6)
        f5.pack(fill="both", expand=True, pady=6)
        self.log_box = scrolledtext.ScrolledText(f5, height=10, font=("Courier", 9),
                                                 state="disabled", bg="#0d1b2a", fg="#90ee90")
        self.log_box.pack(fill="both", expand=True)

        # ── Buttons
        btn_frame = tk.Frame(body, bg="#f0f4f8")
        btn_frame.pack(pady=8)
        self.run_btn = ttk.Button(btn_frame, text="▶  Start validation",
                                  command=self._start, style="Accent.TButton")
        self.run_btn.pack(side="left", padx=8, ipadx=12, ipady=4)
        ttk.Button(btn_frame, text="⏹  Stop",
                   command=self._stop_run).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="🌐  Open result",
                   command=self._open_result).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="💾  Build/Rebuild HPA cache",
                   command=self._rebuild_cache).pack(side="left", padx=4)

    def _log(self, msg: str):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", msg + "\n")
        self.log_box.see("end")
        self.log_box.configure(state="disabled")
        self.update_idletasks()

    def _on_cell_type_changed(self, event=None):
        """When the user changes the cell-type selection, auto-suggest a matching
        cancer cell line category.

        v17 (restored from earlier behavior): the suggestion is applied whenever
        the *current* category is either the placeholder OR was itself put there
        by a previous auto-suggest. Only an explicit user pick (tracked via the
        category combobox's <<ComboboxSelected>> event) is preserved across
        cell-type changes.
        """
        try:
            manual = self.cell_custom.get().strip()
            dropdown = self.cell_var.get()
            if dropdown == CELL_TYPES[0]:
                current_cell = manual
            else:
                current_cell = manual or dropdown
            if current_cell.endswith(" (cancer cell line)"):
                current_cell = current_cell[: -len(" (cancer cell line)")]

            suggested = suggest_category_for_cell_line(current_cell)
            current_cat = self.category_var.get()

            # Decide whether we may overwrite:
            #   1) Currently the placeholder → always safe.
            #   2) Currently equals what we last auto-suggested → still auto-state.
            #   3) Otherwise (user explicitly picked something) → never overwrite.
            placeholder = CANCER_CELL_LINE_CATEGORIES[0]
            is_auto_state = (current_cat == placeholder
                              or current_cat == self._last_auto_category)

            if not is_auto_state:
                return  # respect user's explicit pick

            if suggested:
                self.category_var.set(suggested)
                self._last_auto_category = suggested
            else:
                # No match for this cell type — clear back to placeholder so the
                # user clearly sees that no auto-match was available.
                self.category_var.set(placeholder)
                self._last_auto_category = ""
        except Exception as e:
            # Surface the error in the log — silent failures hide real bugs.
            try:
                self._log(f"⚠️ Auto-category suggestion failed: {e}")
            except Exception:
                pass

    def _on_category_user_pick(self, event=None):
        """User explicitly chose a category from the dropdown.
        Forget the auto-suggest history so subsequent cell-type changes don't
        overwrite the user's pick."""
        # The new selection is by definition NOT auto-suggested any more.
        self._last_auto_category = ""

    def _choose_file(self):
        path = filedialog.askopenfilename(
            title="Select pathway file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if path:
            self._source_file = path
            self.file_var.set(path)
            self._paths = parse_pathway_file(path)
            self.path_count_lbl.config(text=f"→ {len(self._paths)}pathways loaded")
            self._log(f"✅ File loaded: {os.path.basename(path)}  ({len(self._paths)} Path)")

    def _choose_output(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".html",
            filetypes=[("HTML", "*.html"), ("All files", "*.*")],
            initialfile="ppi_validation_result.html"
        )
        if path:
            self.out_var.set(path)

    def _start(self):
        if self._running:
            return
        if not self._paths:
            messagebox.showwarning("Warning", "Please select a pathway file first.")
            return
        # v17: When the dropdown is set to "(Other — use manual field below)",
        # always use the manual entry, never the dropdown value. This avoids the
        # ambiguity where the dropdown still held a stale value from a previous run.
        manual = self.cell_custom.get().strip()
        dropdown = self.cell_var.get()
        if dropdown == CELL_TYPES[0]:  # "(Other — use manual field below)"
            cell_type = manual
        else:
            cell_type = manual or dropdown
        # Strip the parenthetical " (cancer cell line)" suffix on dropdown picks
        # so downstream lookup uses the actual HPA cell line name (e.g. "MDA-MB-468").
        if cell_type.endswith(" (cancer cell line)"):
            cell_type = cell_type[: -len(" (cancer cell line)")]
        if not cell_type:
            messagebox.showwarning("Warning",
                "Please select a cell type from the dropdown OR type one in the manual field.\n\n"
                "If the dropdown is set to '(Other — use manual field below)', the manual field is required.")
            return

        # v13: read cancer cell line category. Empty placeholder → no category search.
        cat_raw = self.category_var.get()
        cancer_category = "" if cat_raw == CANCER_CELL_LINE_CATEGORIES[0] else cat_raw

        self._running = True
        self._stop = False
        self.run_btn.config(state="disabled")
        self.progress["value"] = 0

        t = threading.Thread(target=self._run_validation,
                             args=(cell_type, cancer_category), daemon=True)
        t.start()

    def _stop_run(self):
        self._stop = True
        self._log("⏹  Stop requested...")

    def _open_result(self):
        out = self.out_var.get()
        if os.path.exists(out):
            import webbrowser
            webbrowser.open(f"file://{os.path.abspath(out)}")
        else:
            messagebox.showinfo("Notice", "Run the validation first to create the result file.")

    def _rebuild_cache(self):
        """Manually rebuild the HPA bulk TSV cache.
        Useful when HPA releases a new version or the cache becomes corrupted."""
        if self._running:
            messagebox.showwarning("Busy", "Validation is currently running. Please wait or stop it first.")
            return
        bulk = get_bulk_cache(log_fn=self._log)
        msg = (f"This will download HPA's bulk TSV files (~30 MB total) and rebuild\n"
               f"the local SQLite cache at:\n  {bulk.db_path}\n\n"
               f"Existing data will be replaced. Takes about 1-3 minutes depending\n"
               f"on your internet speed. Continue?")
        if not messagebox.askyesno("Rebuild HPA cache", msg):
            return

        def worker():
            self._log("\n" + "=" * 50)
            self._log("📦 Manual cache rebuild requested")
            self._log("=" * 50)
            self.status_lbl.config(text="Building HPA cache...")
            ok = bulk.build()
            if ok:
                self._log("✅ Cache rebuild complete.")
                messagebox.showinfo("Success", f"HPA cache rebuilt successfully.\n\n{bulk.status_summary()}")
            else:
                self._log("❌ Cache rebuild failed — see log for details.")
                messagebox.showerror("Failed",
                    "Cache rebuild failed.\n\nCheck the log for the error details.\n"
                    "Common causes: no internet, HPA server temporarily down, "
                    "firewall blocking proteinatlas.org.")
            self.status_lbl.config(text="Ready")

        threading.Thread(target=worker, daemon=True).start()

    def _update_progress(self, done: int, total: int):
        pct = (done / total * 100) if total else 0
        self.progress["value"] = pct
        self.status_lbl.config(text=f"Progress: {done}/{total} edges being validated...")
        self.update_idletasks()

    def _run_validation(self, cell_type: str, cancer_category: str = ""):
        try:
            max_p = int(self.max_paths.get())
            paths = self._paths[:max_p]
            self._log(f"\n{'=' * 50}")
            self._log(f"🚀 Validation started | Cell type: {cell_type}")
            self._log(f"   Cancer category (OR-search): {cancer_category or '— (none)'}")
            self._log(f"   Path count: {len(paths)} | API: HPA + STRING + MyGene.info")
            self._log(f"{'=' * 50}\n")

            # ── v16: ensure the bulk HPA cache is populated before validation ──
            bulk = get_bulk_cache(log_fn=self._log)
            if not bulk.is_populated():
                self._log("📦 First-time setup: building local HPA bulk cache.")
                self._log("   This downloads two TSV files (~30 MB total) and takes 1-3 minutes.")
                self._log("   It only happens once — future runs reuse the local SQLite database.")
                self.status_lbl.config(text="Building HPA cache (one-time setup)...")
                ok = bulk.build()
                if not ok:
                    self._log("⚠️ Bulk cache build failed — proceeding with single-entry JSON only.")
                    self._log("   You may still get partial results, but 'Low cancer specificity' genes")
                    self._log("   will report empty cancer-group data.")
                self.status_lbl.config(text="Validating...")
            else:
                self._log(f"💾 {bulk.status_summary()}")

            engine = ValidationEngine(self._log, cancer_category=cancer_category)
            results = engine.validate_paths(
                paths, cell_type,
                progress_fn=self._update_progress,
                stop_flag=lambda: self._stop
            )

            if results:
                out_path = self.out_var.get()
                generate_html(results, cell_type, out_path, self._source_file,
                              cancer_category=cancer_category)
                self._log(f"\n✅ HTML report saved: {out_path}")
                self._log(f"   Validated pathways: {len(results)}")

                grade_counts = defaultdict(int)
                for r in results:
                    grade_counts[r["path_grade"]] += 1
                for g, c in sorted(grade_counts.items()):
                    self._log(f"   {grade_label_kr(g)}: {c}")

                if messagebox.askyesno("Completed", f"✅ Validation completed!\n\nResult file: {out_path}\n\nOpen it in the browser?"):
                    import webbrowser
                    webbrowser.open(f"file://{os.path.abspath(out_path)}")
            else:
                self._log("⚠️ No results were generated.")

        except Exception as e:
            self._log(f"❌ Error: {e}")
            import traceback
            self._log(traceback.format_exc())
        finally:
            self._running = False
            self.run_btn.config(state="normal")
            self.status_lbl.config(text="Completed")
            self.progress["value"] = 100


# ─── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = PPIValidatorApp()
    app.mainloop()