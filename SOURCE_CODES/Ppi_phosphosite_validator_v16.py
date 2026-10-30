
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPI Phosphosite Validator  v1.6 (drop-in fix)
────────────────────────────────────────────────────────────────────────────────
v16 CHANGES (bug-fix release, no behavior/format changes otherwise):
  - OmniPath enzsub: added missing `fields=sources,references` param — this is
    why every enzsub residue hit previously came back with empty sources/
    references and n_references=0.
  - iPTMnet: added the missing trailing slash on `/ptmppi/` and `/substrate/`
    (the official API defines these routes WITH a trailing slash; the
    Actix-web backend 404s without it) — this is why iptmnet_evidence_count
    was 0 for every single edge.
  - iPTMnet: added an _api_errors list (there was none before), so failures
    are now visible instead of being swallowed by a bare `except: return []`.
  - OmniPath / SIGNOR / iPTMnet: exception logging now includes HTTP status
    code and a response-body snippet when available, instead of just the
    exception type/message.
  - SIGNOR: logs explicitly when a response has no data rows, or when none of
    the expected TSV columns (IDB/ENTITYB/TYPEB/MECHANISM/RESIDUE/EFFECT) are
    present in the header — this is the fastest way to tell "SIGNOR changed
    its schema" apart from "genuinely no interaction for this pair".
  - Added an iPTMnet connectivity check to the startup diagnostics (there was
    none before).
  - The full error list (not just the first 8 lines shown in the GUI log) is
    now written to <output_dir>/api_errors.log at the end of every run.
────────────────────────────────────────────────────────────────────────────────
경로(Pathway) 상의 모든 PPI 엣지에 대해 인산화 부위(Phosphosite)의
생물학적 타당성을 4개 전문 DB에서 자동 교차 검증합니다.

  DB 1. OmniPath    — PhosphoSitePlus / SIGNOR / dbPTM 등 30+ 통합
  DB 2. SIGNOR 4.0  — 방향성 인과 관계 + 기전 + 활성화/억제 효과
  DB 3. iPTMnet     — PTM 의존적 PPI (결합 유도 vs 해리)
  DB 4. UniProt     — 단백질별 검증된 PTM 부위 (ECO 코드 포함)

입력 파일 형식:
  - 섹션 1: Gene@Site 관측 인산화 부위 매핑 테이블
  - 섹션 2: "경로 #N: A -> B -> C ..." 경로 목록

출력 파일:
  - phosphosite_validation_report.html
  - phosphosite_validation_results.tsv
  - phosphosite_validation_summary.txt
  - pathway_confidence_scores.tsv

UI: Pathway_Validation_System_v4 호환 스타일
────────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import json
import math
import os
import re
import threading
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    import requests
    _HAS_REQUESTS = True
except ImportError:
    _HAS_REQUESTS = False

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 1: 입력 파일 파서
# ══════════════════════════════════════════════════════════════════════════════

class PathwayFileParser:
    """
    BFS 경로 분석 결과 파일 파서.
    ① observed_phosphosites : {gene: {site, ...}}  관측 인산화 부위
    ② pathways               : [[node, ...], ...]   경로 목록
    ③ unique_edges           : {(geneA, geneB), …} 고유 PPI 쌍
    """

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.observed_phosphosites: Dict[str, Set[str]] = defaultdict(set)
        self.pathways: List[List[str]] = []
        self.unique_edges: Set[Tuple[str, str]] = set()

    def parse(self) -> "PathwayFileParser":
        with open(self.filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        self._parse_phosphosites(content)
        self._parse_pathways(content)
        self._extract_edges()
        return self

    def _parse_phosphosites(self, content: str):
        # Accept S/T/Y prefixed sites (e.g. S396) OR number-only sites (e.g. 396).
        # Number-only entries are stored as-is; matching uses position number comparison.
        pattern = re.compile(r"\|\s*([A-Z][A-Z0-9_]+)@([STY]?\d+)\s*\|", re.MULTILINE)
        for m in pattern.finditer(content):
            gene, site = m.group(1), m.group(2)
            self.observed_phosphosites[gene].add(site)

    def _parse_pathways(self, content: str):
        pattern = re.compile(r"(?:경로|[Pp]ath)\s*#?\d+:\s*(.+)")
        for m in pattern.finditer(content):
            nodes = [n.strip() for n in m.group(1).split("->")]
            nodes = [n for n in nodes if n]
            if len(nodes) >= 2:
                self.pathways.append(nodes)
        # 경로 # 형식 없이 -> 로만 연결된 줄도 허용
        if not self.pathways:
            for line in content.splitlines():
                line = line.strip()
                if "->" in line and not line.startswith("|") and not line.startswith("#"):
                    nodes = [n.strip() for n in line.split("->")]
                    nodes = [n for n in nodes if n]
                    if len(nodes) >= 2:
                        self.pathways.append(nodes)

    def _extract_edges(self):
        for path in self.pathways:
            for i in range(len(path) - 1):
                a, b = path[i], path[i + 1]
                if a != b:
                    self.unique_edges.add((a, b))


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 2: UniProt Accession 캐시 및 PTM 조회
# ══════════════════════════════════════════════════════════════════════════════

class UniProtAccessionCache:
    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

    def __init__(self, cache_file: str = "uniprot_accession_cache.json"):
        self.cache_file = cache_file
        self.cache: Dict[str, Optional[str]] = {}
        self._ptm_cache: Dict[str, List[Dict]] = {}
        self._load()

    def _load(self):
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    self.cache = data.get("acc", {})
                    self._ptm_cache = data.get("ptm", {})
            except Exception:
                pass

    def _save(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump({"acc": self.cache, "ptm": self._ptm_cache}, f, indent=2)
        except Exception:
            pass

    def get_accession(self, gene: str) -> Optional[str]:
        if gene in self.cache:
            return self.cache[gene]
        if not _HAS_REQUESTS:
            self.cache[gene] = None
            return None
        params = {
            "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
            "format": "json",
            "size": 1,
            "fields": "accession,gene_names",
        }
        try:
            r = requests.get(self.BASE_URL, params=params, timeout=10)
            r.raise_for_status()
            results = r.json().get("results", [])
            acc = results[0]["primaryAccession"] if results else None
            self.cache[gene] = acc
            self._save()
            time.sleep(0.2)
            return acc
        except Exception:
            self.cache[gene] = None
            return None

    def get_uniprot_ptm_sites(self, gene: str) -> List[Dict]:
        if gene in self._ptm_cache:
            return self._ptm_cache[gene]
        acc = self.get_accession(gene)
        if not acc or not _HAS_REQUESTS:
            self._ptm_cache[gene] = []
            return []
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        try:
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            data = r.json()
            sites = []
            # ── PTM 코멘트 전체 수집 (기능적 맥락) ─────────────────────────
            ptm_comment_text = ""
            for c in data.get("comments", []):
                if c.get("commentType") == "PTM":
                    for t in c.get("texts", []):
                        ptm_comment_text += " " + t.get("value", "")
            # ── Modified residue features ────────────────────────────────
            for feat in data.get("features", []):
                if feat.get("type") != "Modified residue":
                    continue
                desc_full = feat.get("description", "")
                if "Phospho" not in desc_full:
                    continue
                pos = feat.get("location", {}).get("start", {}).get("value")
                if pos is None:
                    continue
                eco  = (feat.get("evidences") or [{}])[0].get("evidenceCode", "")
                # 기본 설명: 첫 번째 세미콜론 앞
                desc_base = desc_full.split(";")[0].strip()  # e.g. "Phosphotyrosine"
                # 나머지 세미콜론 항목들: 키나아제명, 효과 등
                parts = [p.strip() for p in desc_full.split(";")[1:] if p.strip()]
                # 키나아제 추출: "by EGFR" / "by SRC and LCK"
                kinases = []
                for p in parts:
                    m = re.match(r'^by\s+(.+)$', p, re.IGNORECASE)
                    if m:
                        kinases = [k.strip() for k in re.split(r'\s+and\s+|,', m.group(1))]
                # 기능 키워드 추출: activates / inhibits / required for / blocks etc.
                func_notes = [p for p in parts if not p.lower().startswith("by ")
                              and p.lower() not in ("alternate", "by phosphorylation")]
                sites.append({
                    "position":      pos,
                    "description":   desc_base,
                    "description_full": desc_full,
                    "kinases":       kinases,          # e.g. ["EGFR", "SRC"]
                    "func_notes":    func_notes,       # e.g. ["activates kinase activity"]
                    "evidence_code": eco,
                    "ptm_comment":   ptm_comment_text.strip()[:300],
                })
            self._ptm_cache[gene] = sites
            self._save()
            time.sleep(0.2)
            return sites
        except Exception:
            self._ptm_cache[gene] = []
            return []


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 3: OmniPath API
# ══════════════════════════════════════════════════════════════════════════════

class OmniPathClient:
    # enzsub: 키나아제-기질 인산화 (잔기 수준)
    ENZSUB_URL      = "https://omnipathdb.org/enzsub"
    # interactions: 방향성 신호전달 상호작용 (잔기 없지만 방향·기전 포함)
    INTERACTIONS_URL = "https://omnipathdb.org/interactions"
    _api_errors: List[str] = []  # collect errors for diagnostics

    def __init__(self, cache_file: str = "omnipath_cache.json"):
        self.cache_file = cache_file
        self.cache: Dict[str, List[Dict]] = {}
        self._load()

    def _load(self):
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except Exception:
                pass

    def _save(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)
        except Exception:
            pass

    @staticmethod
    def _parse_sources(raw) -> List[str]:
        """OmniPath sources 필드: 문자열(';' 구분) 또는 리스트 모두 처리"""
        if isinstance(raw, list):
            return [s.strip() for s in raw if s.strip()]
        if isinstance(raw, str) and raw.strip():
            # 세미콜론 또는 쉼표로 구분
            sep = ";" if ";" in raw else ","
            return [s.strip() for s in raw.split(sep) if s.strip()]
        return []

    def get_enzyme_substrate(self, enzyme_gene: str, substrate_gene: str,
                             api_delay: float = 0.3) -> List[Dict]:
        """
        OmniPath enzsub 조회: 키나아제(A)가 기질(B)을 인산화하는 잔기 수준 정보.
        반환 결과에 잔기 정보가 없으면 interactions 엔드포인트로 보완.
        """
        key = f"enzsub__{enzyme_gene}__{substrate_gene}"
        if key in self.cache:
            return self.cache[key]
        if not _HAS_REQUESTS:
            self.cache[key] = []
            return []

        results = []

        # ① enzsub endpoint (residue-level kinase-substrate)
        try:
            r = requests.get(self.ENZSUB_URL, params={
                "enzymes":      enzyme_gene,
                "substrates":   substrate_gene,
                "modification": "phosphorylation",
                "genesymbols":  "yes",
                "format":       "json",
                "organisms":    9606,
                # v16 FIX: sources/references are NOT included by default by the
                # OmniPath webservice — they must be requested explicitly via
                # `fields`, exactly like the interactions call below already does.
                # Without this, every enzsub hit came back with sources=[],
                # references=[], n_references=0 even when residues were found.
                "fields":       "sources,references",
            }, timeout=20)
            r.raise_for_status()
            raw_list = r.json()
            if not isinstance(raw_list, list):
                OmniPathClient._api_errors.append(
                    f"OmniPath enzsub [{enzyme_gene}→{substrate_gene}]: "
                    f"unexpected non-list JSON response (type={type(raw_list).__name__}); "
                    f"raw snippet={str(raw_list)[:200]!r}")
                raw_list = []
            for item in raw_list:
                # Client-side substrate/enzyme check: API may return more rows than expected
                enz_sym = item.get("enzyme_genesymbol", "")
                sub_sym = item.get("substrate_genesymbol", "")
                if enz_sym and enz_sym.upper() != enzyme_gene.upper():
                    continue
                if sub_sym and sub_sym.upper() != substrate_gene.upper():
                    continue

                # sources field: semicolon-separated string or list
                raw_src = item.get("sources") or item.get("sources_curated") or ""
                srcs = self._parse_sources(raw_src)

                raw_ref = item.get("references") or item.get("references_curated") or ""
                refs = self._parse_sources(raw_ref)

                results.append({
                    "residue_type":   item.get("residue_type", ""),
                    "residue_offset": item.get("residue_offset"),
                    "modification":   item.get("modification", "phosphorylation"),
                    "sources":        srcs,
                    "references":     refs,
                    "n_references":   len(refs) if refs else int(item.get("n_references", 0) or 0),
                    "endpoint":       "enzsub",
                })
        except Exception as _e:
            # v16 FIX: include HTTP status + a short body snippet when available,
            # so a wrong format/endpoint/redirect shows up clearly in the log
            # instead of just "JSONDecodeError: Expecting value" with no context.
            _detail = ""
            _resp = getattr(_e, "response", None)
            if _resp is not None:
                _detail = f" | status={_resp.status_code} | body={_resp.text[:200]!r}"
            OmniPathClient._api_errors.append(
                f"OmniPath enzsub [{enzyme_gene}→{substrate_gene}]: {type(_e).__name__}: {_e}{_detail}")

        # ② interactions endpoint — ALWAYS called regardless of enzsub results
        # Provides directionality, stimulation/inhibition and additional source coverage
        try:
            r2 = requests.get(self.INTERACTIONS_URL, params={
                "sources":     enzyme_gene,
                "targets":     substrate_gene,
                "genesymbols": "1",
                "organisms":   "9606",
                "fields":      "sources,references,consensus_stimulation,consensus_inhibition",
                "format":      "json",
            }, timeout=20)
            r2.raise_for_status()
            raw2 = r2.json()
            if isinstance(raw2, list):
                for item in raw2:
                    if not isinstance(item, dict):   # list-of-lists 응답 방어
                        continue
                    srcs = self._parse_sources(item.get("sources") or "")
                    refs = self._parse_sources(item.get("references") or "")
                    results.append({
                        "residue_type":          "",
                        "residue_offset":        None,
                        "modification":          "interaction",
                        "sources":               srcs,
                        "references":            refs,
                        "n_references":          len(refs),
                        "consensus_stimulation": item.get("consensus_stimulation", ""),
                        "consensus_inhibition":  item.get("consensus_inhibition", ""),
                        "endpoint":              "interactions",
                    })
            else:
                # v16 FIX: previously silently dropped; now logged so a
                # persistently-empty interactions endpoint is diagnosable.
                OmniPathClient._api_errors.append(
                    f"OmniPath interactions [{enzyme_gene}→{substrate_gene}]: "
                    f"unexpected non-list JSON response (type={type(raw2).__name__}); "
                    f"raw snippet={str(raw2)[:200]!r}")
        except Exception as _e2:
            _detail = ""
            _resp2 = getattr(_e2, "response", None)
            if _resp2 is not None:
                _detail = f" | status={_resp2.status_code} | body={_resp2.text[:200]!r}"
            OmniPathClient._api_errors.append(
                f"OmniPath interactions [{enzyme_gene}→{substrate_gene}]: {type(_e2).__name__}: {_e2}{_detail}")

        self.cache[key] = results
        self._save()
        time.sleep(api_delay)
        return results


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 4: SIGNOR 4.0 API
# ══════════════════════════════════════════════════════════════════════════════

class SIGNORClient:
    """
    SIGNOR 4.0 API client.
    Official REST API: /getData.php?id=<UniProtKB_ID>&organism=9606  (TSV)
    Query strategy:
      1. Convert gene_a -> UniProt accession (via UniProtAccessionCache)
      2. GET /getData.php?id=<acc_a>&organism=9606 -> all interactions of protein A
      3. Filter rows where IDB == acc_b (gene_b UniProt) or ENTITYB == gene_b
      4. Further filter for phosphorylation mechanism
    Note: getInteractionByEntities.php is NOT an official SIGNOR endpoint (returns 404).
    """
    BASE_URL    = "https://signor.uniroma2.it/getData.php"
    _api_errors: List[str] = []

    def __init__(self, cache_file: str = "signor_cache.json",
                 uniprot_cache=None):
        self.cache_file    = cache_file
        self.uniprot_cache = uniprot_cache
        self.cache: Dict[str, List[Dict]] = {}
        self._load()

    def _load(self):
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except Exception:
                pass

    def _save(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)
        except Exception:
            pass

    def _gene_to_uniprot(self, gene: str) -> str:
        """Return primary reviewed UniProt accession for a human gene symbol."""
        if self.uniprot_cache:
            return self.uniprot_cache.get_accession(gene)
        if not _HAS_REQUESTS:
            return ""
        try:
            r = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                        "fields": "accession", "format": "json", "size": 1},
                timeout=10)
            r.raise_for_status()
            res = r.json().get("results", [])
            if res:
                return res[0].get("primaryAccession", "")
        except Exception:
            pass
        return ""

    def get_interaction(self, gene_a: str, gene_b: str,
                        api_delay: float = 0.3) -> List[Dict]:
        """
        Return phosphorylation interactions gene_a -> gene_b from SIGNOR.
        Uses /getData.php?id=<UniProt_A>&organism=9606 then filters for gene_b.
        """
        key = f"signor4__{gene_a}__{gene_b}"
        if key in self.cache:
            return self.cache[key]
        if not _HAS_REQUESTS:
            self.cache[key] = []
            return []

        results = []

        # Step 1: resolve UniProt accessions
        acc_a = self._gene_to_uniprot(gene_a)
        acc_b = self._gene_to_uniprot(gene_b)

        if not acc_a:
            SIGNORClient._api_errors.append(
                f"SIGNOR: could not resolve UniProt ID for {gene_a}")
            self.cache[key] = []
            self._save()
            return []

        # Step 2: fetch all interactions of protein A
        try:
            r = requests.get(self.BASE_URL,
                             params={"id": acc_a, "organism": "9606"},
                             timeout=20)
            r.raise_for_status()
            text = r.text.strip()
            if not text:
                SIGNORClient._api_errors.append(
                    f"SIGNOR getData [{gene_a}({acc_a})->{gene_b}]: empty response body")
                self.cache[key] = []
                self._save()
                return []

            lines = text.split("\n")
            if len(lines) < 2:
                SIGNORClient._api_errors.append(
                    f"SIGNOR getData [{gene_a}({acc_a})->{gene_b}]: response had no data rows "
                    f"(only {len(lines)} line(s)); first line={lines[0][:200]!r}")
                self.cache[key] = []
                self._save()
                return []

            header = [h.strip() for h in lines[0].split("\t")]
            # v16 FIX: if SIGNOR changed its column names, IDB/ENTITYB/TYPEB below
            # would silently return "" for every row and matched_b would always be
            # False -- log that case instead of failing silently.
            _expected_cols = {"IDB", "ENTITYB", "TYPEB", "MECHANISM", "RESIDUE", "EFFECT"}
            if not _expected_cols.intersection(header):
                SIGNORClient._api_errors.append(
                    f"SIGNOR getData [{gene_a}({acc_a})->{gene_b}]: none of the expected "
                    f"columns {sorted(_expected_cols)} found in response header={header} "
                    f"— SIGNOR may have changed its TSV schema or returned an error page")

            for line in lines[1:]:
                if not line.strip():
                    continue
                cols = line.split("\t")
                row  = dict(zip(header, cols))

                # Filter: target must be gene_b
                # Try IDB (UniProt accession) and ENTITYB (gene/protein name)
                idb   = row.get("IDB",     "").strip()
                entb  = row.get("ENTITYB", "").strip()
                typeb = row.get("TYPEB",   "").strip()

                # Only consider protein targets (skip complexes, phenotypes, stimuli)
                if typeb and typeb.upper() not in ("PROTEIN", "PROTEINFAMILY", ""):
                    continue

                matched_b = ((acc_b and idb == acc_b) or
                             entb.upper() == gene_b.upper())
                if not matched_b:
                    continue

                # Keep phosphorylation AND rows with empty mechanism
                # (SIGNOR uses "phosphorylation", "phospho-X", or leaves it blank)
                mech = row.get("MECHANISM", "").strip()
                is_phospho = (not mech) or ("phospho" in mech.lower())
                if not is_phospho:
                    continue

                residue = row.get("RESIDUE", "").strip()
                # SIGNOR residue format: e.g. "Y705", "T308", "S473"
                # Normalize: strip leading "pS"/"pT"/"pY" prefix if present
                if residue.startswith("p") and len(residue) > 1 and residue[1].isalpha():
                    residue = residue[1:]

                results.append({
                    "mechanism": mech or "phosphorylation",
                    "residue":   residue,
                    "effect":    row.get("EFFECT",  "").strip(),
                    "pubmed_id": row.get("PMID",    "").strip(),
                    "direct":    row.get("DIRECT",  "").strip(),
                    "score":     row.get("SCORE",   "").strip(),
                })

        except Exception as _e:
            _detail = ""
            _resp = getattr(_e, "response", None)
            if _resp is not None:
                _detail = f" | status={_resp.status_code} | body={_resp.text[:200]!r}"
            SIGNORClient._api_errors.append(
                f"SIGNOR getData [{gene_a}({acc_a})->{gene_b}]: {type(_e).__name__}: {_e}{_detail}")

        self.cache[key] = results
        self._save()
        time.sleep(api_delay)
        return results

# ══════════════════════════════════════════════════════════════════════════════
# 섹션 5: iPTMnet API
# ══════════════════════════════════════════════════════════════════════════════

class iPTMnetClient:
    BASE_URL = "https://research.bioinformatics.udel.edu/iptmnet/api"
    # v16 FIX: no diagnostics existed for this client at all — every failure was
    # swallowed by a bare `except Exception: return []`, so 0 results anywhere
    # looked identical to "no data" and to "server unreachable". Added the same
    # _api_errors mechanism the OmniPath/SIGNOR clients already have.
    _api_errors: List[str] = []

    def __init__(self, uniprot_cache: UniProtAccessionCache,
                 cache_file: str = "iptmnet_cache.json"):
        self.uniprot_cache = uniprot_cache
        self.cache_file = cache_file
        self.cache: Dict[str, List[Dict]] = {}
        self._load()

    def _load(self):
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except Exception:
                pass

    def _save(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)
        except Exception:
            pass

    def get_ptm_ppi(self, gene: str, api_delay: float = 0.3) -> List[Dict]:
        if gene in self.cache:
            return self.cache[gene]
        acc = self.uniprot_cache.get_accession(gene)
        if not acc or not _HAS_REQUESTS:
            self.cache[gene] = []
            return []
        try:
            # v16 FIX: the official iPTMnet API endpoint table lists this route
            # WITH a trailing slash ("/{id}/ptmppi/"). The iPTMnet backend is a
            # Rust/Actix-web service, which by default does NOT redirect a
            # missing-trailing-slash request to the slashed route — it 404s.
            # Without the slash, every single call here was failing and being
            # silently swallowed below, which is why iptmnet_evidence_count
            # was 0 for all 367 edges.
            r = requests.get(f"{self.BASE_URL}/{acc}/ptmppi/", timeout=15)
            r.raise_for_status()
            raw = r.json()
            if not isinstance(raw, list):
                iPTMnetClient._api_errors.append(
                    f"iPTMnet ptmppi [{gene}({acc})]: unexpected non-list JSON "
                    f"(type={type(raw).__name__}); raw snippet={str(raw)[:200]!r}")
                raw = []
            results = [
                {
                    "ptm_type": item.get("ptm_type", ""),
                    "site": item.get("site", ""),
                    "interactant_id": item.get("interactant_id", ""),
                    "interactant_name": item.get("interactant_name", ""),
                    "association": item.get("association", ""),
                    "pubmed_id": item.get("pubmed_id", ""),
                }
                for item in raw
                if "phospho" in item.get("ptm_type", "").lower()
            ]
            self.cache[gene] = results
            self._save()
            time.sleep(api_delay)
            return results
        except Exception as _e:
            _detail = ""
            _resp = getattr(_e, "response", None)
            if _resp is not None:
                _detail = f" | status={_resp.status_code} | body={_resp.text[:200]!r}"
            iPTMnetClient._api_errors.append(
                f"iPTMnet ptmppi [{gene}({acc})]: {type(_e).__name__}: {_e}{_detail}")
            self.cache[gene] = []
            return []

    def get_substrates(self, gene: str, api_delay: float = 0.3) -> List[Dict]:
        key = f"__sub__{gene}"
        if key in self.cache:
            return self.cache[key]
        acc = self.uniprot_cache.get_accession(gene)
        if not acc or not _HAS_REQUESTS:
            self.cache[key] = []
            return []
        try:
            # v16 FIX: trailing slash, see get_ptm_ppi() above for the reasoning.
            r = requests.get(f"{self.BASE_URL}/{acc}/substrate/", timeout=15)
            r.raise_for_status()
            raw = r.json()
            if not isinstance(raw, list):
                iPTMnetClient._api_errors.append(
                    f"iPTMnet substrate [{gene}({acc})]: unexpected non-list JSON "
                    f"(type={type(raw).__name__}); raw snippet={str(raw)[:200]!r}")
                raw = []
            results = [
                {
                    "ptm_type": item.get("ptm_type", ""),
                    "substrate_id": item.get("substrate_id", ""),
                    "substrate_name": item.get("substrate_name", ""),
                    "site": item.get("site", ""),
                    "source": item.get("source", ""),
                }
                for item in raw
                if "phospho" in item.get("ptm_type", "").lower()
            ]
            self.cache[key] = results
            self._save()
            time.sleep(api_delay)
            return results
        except Exception as _e:
            _detail = ""
            _resp = getattr(_e, "response", None)
            if _resp is not None:
                _detail = f" | status={_resp.status_code} | body={_resp.text[:200]!r}"
            iPTMnetClient._api_errors.append(
                f"iPTMnet substrate [{gene}({acc})]: {type(_e).__name__}: {_e}{_detail}")
            self.cache[key] = []
            return []


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 5.5: 인산화 부위 기능 주석 클라이언트  (v1.5 — UniProt 기반)
# ══════════════════════════════════════════════════════════════════════════════
#
# 각 인산화 부위의 생물학적 역할(활성화/억제/결합/위치이동 등)을
# 이미 작동 중인 UniProt REST API를 주 소스로 조회합니다.
# PSP AJAX / PhosphoELM은 보조 시도 (접근 가능할 때만 사용).
#
# 반환 필드 (사이트별):
#   effect       : activation / inhibition / binding / localization / unknown
#   on_function  : 구체적 기능 설명 (예: "activates kinase activity")
#   kinases      : 해당 부위를 인산화하는 키나아제 목록
#   ptm_comment  : UniProt PTM 코멘트 전체 (맥락 파악용)
#   n_refs       : 참조 문헌 수 (PSP 접근 가능할 때만)
#   source       : UniProt / PhosphoSitePlus / PhosphoELM
# ══════════════════════════════════════════════════════════════════════════════

class PhosphoSitePlusClient:
    """
    인산화 부위 기능 주석 클라이언트 (UniProt 기반, PSP/PhosphoELM 보조).

    조회 전략:
      1. UniProt REST API (주) — 이미 정상 작동, 기능 주석 포함
      2. PSP AJAX endpoint (보조) — 접근 가능 시 추가 정보
      3. PhosphoELM REST (보조) — PSP 실패 시
    """

    PSP_RESIDUE_URL = "https://www.phosphosite.org/ajax/residueList.action"
    PHELM_BASE_URL  = "http://phospho.elm.eu.org/byAccession"

    # 기능 키워드 → canonical effect 레이블
    _EFFECT_MAP = {
        "activat":   "activation",
        "stimulat":  "activation",
        "upregulat": "activation",
        "required for": "activation",
        "promotes":  "activation",
        "enhances":  "activation",
        "inhibit":   "inhibition",
        "downregul": "inhibition",
        "suppress":  "inhibition",
        "blocks":    "inhibition",
        "prevent":   "inhibition",
        "reduces":   "inhibition",
        "bind":      "binding",
        "interact":  "binding",
        "recruit":   "binding",
        "associat":  "binding",
        "docking":   "binding",
        "local":     "localization",
        "transloc":  "localization",
        "nuclear":   "localization",
        "cytoplasm": "localization",
        "export":    "localization",
        "import":    "localization",
        "degradat":  "degradation",
        "ubiquitin": "degradation",
        "stability": "stability",
        "stabiliz":  "stability",
        "turnover":  "stability",
    }

    EFFECT_COLORS = {
        "activation":   "#27ae60",
        "inhibition":   "#c0392b",
        "binding":      "#2980b9",
        "localization": "#8e44ad",
        "degradation":  "#e67e22",
        "stability":    "#16a085",
        "unknown":      "#95a5a6",
    }

    _api_errors: List[str] = []

    def __init__(self, uniprot_cache: "UniProtAccessionCache",
                 cache_file: str = "psp_function_cache.json"):
        self.uniprot_cache = uniprot_cache
        self.cache_file = cache_file
        self.cache: Dict[str, List[Dict]] = {}
        self._load()

    def _load(self):
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except Exception:
                pass

    def _save(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)
        except Exception:
            pass

    @staticmethod
    def _norm_site(raw: str) -> str:
        """'pY705' / 'Y705-p' / '705' → 'Y705'"""
        m = re.search(r'([STY])(\d+)', str(raw).upper())
        if m:
            return m.group(1) + m.group(2)
        m2 = re.search(r'(\d+)', str(raw))
        return m2.group(1) if m2 else ""

    @classmethod
    def _classify_effect(cls, text: str) -> str:
        t = text.lower()
        for kw, label in cls._EFFECT_MAP.items():
            if kw in t:
                return label
        return "unknown"

    # ── 1. UniProt 기반 주석 (주 소스) ──────────────────────────────────────
    def _from_uniprot(self, gene: str) -> List[Dict]:
        """UniProt의 확장 PTM features에서 기능 주석 추출 (이미 캐시됨)."""
        raw_sites = self.uniprot_cache.get_uniprot_ptm_sites(gene)
        results = []
        for s in raw_sites:
            pos  = s.get("position")
            if pos is None:
                continue
            desc_base  = s.get("description", "")   # e.g. "Phosphotyrosine"
            desc_full  = s.get("description_full", desc_base)
            func_notes = s.get("func_notes", [])
            kinases    = s.get("kinases", [])
            ptm_cmt    = s.get("ptm_comment", "")

            # 아미노산 코드 추출: Phosphotyrosine→Y, Phosphoserine→S, Phosphothreonine→T
            aa = ""
            if "tyrosine"    in desc_base.lower(): aa = "Y"
            elif "serine"    in desc_base.lower(): aa = "S"
            elif "threonine" in desc_base.lower(): aa = "T"
            site_code = f"{aa}{pos}" if aa else str(pos)

            # 기능 텍스트 조합: func_notes + ptm_comment 관련 문장
            func_text = " ".join(func_notes)
            # PTM 코멘트에서 이 부위 관련 문장 추출
            if ptm_cmt and site_code:
                for sent in re.split(r'[.;]', ptm_cmt):
                    if site_code in sent or (aa and str(pos) in sent):
                        func_text += " " + sent.strip()

            effect = self._classify_effect(func_text + " " + desc_full) if func_text else "unknown"

            results.append({
                "site":        site_code,
                "effect":      effect,
                "on_function": func_text[:150] or desc_full[:100],
                "on_process":  "",
                "domain":      "",
                "n_refs":      0,
                "diseases":    "",
                "kinases":     kinases,
                "source":      "UniProt",
            })
        return results

    # ── 2. PSP AJAX 보조 시도 ────────────────────────────────────────────────
    def _try_psp_ajax(self, acc: str) -> Optional[List[Dict]]:
        """PSP AJAX endpoint 시도. 실패 시 None 반환 (오류 로그만 남김)."""
        try:
            r = requests.get(
                self.PSP_RESIDUE_URL,
                params={"uniprotAccs": acc, "types": "Phosphorylation"},
                headers={"X-Requested-With": "XMLHttpRequest",
                         "Referer": "https://www.phosphosite.org/"},
                timeout=10,
            )
            r.raise_for_status()
            data = r.json()
            items = data if isinstance(data, list) else data.get("sites", data.get("data", []))
            results = []
            for item in items:
                if not isinstance(item, dict):
                    continue
                site_raw = (item.get("residue") or item.get("site") or
                            item.get("modification_site") or "")
                site = self._norm_site(site_raw)
                if not site:
                    continue
                on_func = str(item.get("on_function", "") or "")
                on_proc = str(item.get("on_process", "") or "")
                effect  = self._classify_effect(on_func + " " + on_proc)
                n_refs  = int(item.get("ms_cst", 0) or item.get("cst_refs", 0) or 0)
                results.append({
                    "site":        site,
                    "effect":      effect,
                    "on_function": on_func[:150],
                    "on_process":  on_proc[:80],
                    "domain":      str(item.get("domain", ""))[:60],
                    "n_refs":      n_refs,
                    "diseases":    str(item.get("on_other", ""))[:100],
                    "kinases":     [],
                    "source":      "PhosphoSitePlus",
                })
            return results if results else None
        except Exception as e:
            self._api_errors.append(f"PSP AJAX [{acc}]: {e}")
            return None

    # ── 3. PhosphoELM 보조 시도 ─────────────────────────────────────────────
    def _try_phelm(self, acc: str) -> Optional[List[Dict]]:
        """PhosphoELM /byAccession/<acc>/ 시도."""
        for url in [
            f"{self.PHELM_BASE_URL}/{acc}/",
            f"http://phospho.elm.eu.org/phosphosite/byAccession/{acc}/"
        ]:
            try:
                r = requests.get(url, headers={"Accept": "application/json"}, timeout=10)
                if r.status_code == 404:
                    continue
                r.raise_for_status()
                data = r.json()
                items = data if isinstance(data, list) else data.get("objects", data.get("entries", []))
                results = []
                for item in items:
                    if not isinstance(item, dict):
                        continue
                    pos_raw = item.get("position", "") or item.get("residue_number", "")
                    aa_raw  = item.get("code", "") or item.get("amino_acid", "")
                    site    = self._norm_site(f"{aa_raw}{pos_raw}" if aa_raw else str(pos_raw))
                    if not site:
                        continue
                    note   = str(item.get("notes", "") or "")
                    kinase = str(item.get("kinase", "") or "")
                    effect = self._classify_effect(note + " " + kinase)
                    n_refs = len(item.get("references", []))
                    results.append({
                        "site":        site,
                        "effect":      effect,
                        "on_function": (note or kinase)[:150],
                        "on_process":  "",
                        "domain":      str(item.get("domain", ""))[:60],
                        "n_refs":      n_refs,
                        "diseases":    "",
                        "kinases":     [kinase] if kinase else [],
                        "source":      "PhosphoELM",
                    })
                if results:
                    return results
            except Exception as e:
                self._api_errors.append(f"PhosphoELM [{acc}]: {e}")
        return None

    # ── 공개 메서드: 통합 조회 ───────────────────────────────────────────────
    def get_site_functions(self, gene: str, api_delay: float = 0.3) -> List[Dict]:
        """
        gene의 인산화 부위별 기능 주석 반환.
        캐시 → UniProt(주) → PSP AJAX(보조) → PhosphoELM(보조) 순으로 통합.
        PSP/PhosphoELM 데이터는 UniProt 항목의 n_refs/diseases 필드를 보완.
        """
        if gene in self.cache:
            return self.cache[gene]

        # 1. UniProt 기반 (항상 시도)
        results = self._from_uniprot(gene)

        if _HAS_REQUESTS:
            acc = self.uniprot_cache.get_accession(gene)
            if acc:
                # 2. PSP AJAX 보조
                psp_data = self._try_psp_ajax(acc)
                # 3. PhosphoELM 보조 (PSP 실패 시)
                if psp_data is None:
                    time.sleep(api_delay * 0.5)
                    psp_data = self._try_phelm(acc)

                if psp_data:
                    # PSP/PhosphoELM 데이터로 UniProt 항목 보완
                    psp_map = {self._norm_site(p["site"]): p for p in psp_data}
                    # UniProt 결과 업데이트
                    updated = set()
                    for r in results:
                        ns = self._norm_site(r["site"])
                        if ns in psp_map:
                            p = psp_map[ns]
                            if r["effect"] == "unknown" and p["effect"] != "unknown":
                                r["effect"] = p["effect"]
                            if not r["on_function"] and p["on_function"]:
                                r["on_function"] = p["on_function"]
                            if p["n_refs"] > 0:
                                r["n_refs"]   = p["n_refs"]
                            if p["diseases"]:
                                r["diseases"] = p["diseases"]
                            if p["source"] != "UniProt":
                                r["source"] = f"UniProt+{p['source']}"
                            updated.add(ns)
                    # PSP에만 있는 부위 추가 (UniProt에 없는 경우)
                    existing = {self._norm_site(r["site"]) for r in results}
                    for ns, p in psp_map.items():
                        if ns not in existing:
                            results.append(p)

        self.cache[gene] = results
        self._save()
        time.sleep(api_delay)
        return results

    def lookup_site(self, gene: str, site_code: str,
                    api_delay: float = 0.3) -> Optional[Dict]:
        """특정 site의 기능 주석 딕셔너리 반환."""
        all_sites = self.get_site_functions(gene, api_delay)
        norm = self._norm_site(site_code)
        for entry in all_sites:
            if self._norm_site(entry.get("site", "")) == norm:
                return entry
        # 숫자만 비교
        pos = re.search(r'\d+', norm)
        if pos:
            p = pos.group()
            for entry in all_sites:
                ep = re.search(r'\d+', entry.get("site", ""))
                if ep and ep.group() == p:
                    return entry
        return None

    def annotate_matched_sites(self, gene: str, matched_sites: List[Dict],
                               api_delay: float = 0.3) -> List[Dict]:
        """matched_sites 각 항목에 기능 주석 필드 추가하여 반환."""
        annotated = []
        for m in matched_sites:
            entry = dict(m)
            site_code = m.get("site", "") or m.get("observed", "")
            psp = self.lookup_site(gene, site_code, api_delay)
            if psp:
                entry["psp_effect"]   = psp.get("effect", "unknown")
                entry["psp_function"] = psp.get("on_function", "")
                entry["psp_process"]  = psp.get("on_process", "")
                entry["psp_domain"]   = psp.get("domain", "")
                entry["psp_n_refs"]   = psp.get("n_refs", 0)
                entry["psp_diseases"] = psp.get("diseases", "")
                entry["psp_kinases"]  = ", ".join(psp.get("kinases", []))
                entry["psp_source"]   = psp.get("source", "")
            else:
                entry["psp_effect"]   = ""
                entry["psp_function"] = ""
                entry["psp_process"]  = ""
                entry["psp_domain"]   = ""
                entry["psp_n_refs"]   = 0
                entry["psp_diseases"] = ""
                entry["psp_kinases"]  = ""
                entry["psp_source"]   = ""
            annotated.append(entry)
        return annotated


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 6: 검증 엔진
# ══════════════════════════════════════════════════════════════════════════════

class EdgeVerdict:
    GRADE_A = "A: Exact phosphosite match"
    GRADE_B = "B: Theoretical site exists (observed mismatch)"
    GRADE_C = "C: Non-phosphorylation PPI (e.g., scaffolding)"
    GRADE_D = "D: No evidence (potential false positive)"

    def __init__(self, gene_a: str, gene_b: str):
        self.gene_a = gene_a
        self.gene_b = gene_b
        self.grade = self.GRADE_D
        self.observed_sites_a: List[str] = []
        self.observed_sites_b: List[str] = []
        self.omnipath_sites: List[Dict] = []
        self.signor_data: List[Dict] = []
        self.iptmnet_data: List[Dict] = []
        self.uniprot_ptm_a: List[Dict] = []
        self.uniprot_ptm_b: List[Dict] = []
        self.matched_sites: List[Dict] = []      # Tier1: B observed ↔ DB substrate residues
        self.matched_sites_a: List[Dict] = []    # Tier2: A observed ↔ UniProt PTM-A (kinase activation)
        self.theoretical_sites: List[str] = []
        self.reason: str = ""
        self.confidence_score: float = 0.0
        self.pathway_appearances: int = 0
        # v1.5: PSP 기능 주석 (matched_sites에 psp_* 필드가 추가된 버전)
        self.matched_sites_annotated: List[Dict] = []  # matched_sites + PSP annotations

    def to_dict(self) -> Dict:
        # OmniPath sources: 리스트이므로 join 처리
        def _src_str(sites):
            raw = set()
            for o in sites:
                for s in o.get("sources", []):
                    if s:
                        raw.add(s)
            return "; ".join(sorted(raw))[:120]

        enzsub_with_res = [o for o in self.omnipath_sites if o.get("residue_offset") and o.get("endpoint") == "enzsub"]
        interactions    = [o for o in self.omnipath_sites if o.get("endpoint") == "interactions"]

        return {
            "edge":                     f"{self.gene_a} -> {self.gene_b}",
            "gene_a":                   self.gene_a,
            "gene_b":                   self.gene_b,
            "grade":                    self.grade,
            "confidence_score":         round(self.confidence_score, 3),
            "pathway_appearances":      self.pathway_appearances,
            "observed_sites_a":         "; ".join(self.observed_sites_a),
            "observed_sites_b":         "; ".join(self.observed_sites_b),
            "matched_sites":            "; ".join(f"{m['site']}({m['source']})" for m in self.matched_sites),
            "matched_sites_a":          "; ".join(f"{m['site']}({m['source']})" for m in self.matched_sites_a),
            "theoretical_sites":        "; ".join(s for s in self.theoretical_sites if not s.startswith(("UniProt:","iPTM:"))),
            "theoretical_sites_uniprot":"; ".join(s.replace("UniProt:","") for s in self.theoretical_sites if s.startswith("UniProt:")),
            "enzsub_residue_count":     len(enzsub_with_res),
            "interaction_record_count": len(interactions),
            "signor_evidence_count":    len(self.signor_data),
            "iptmnet_evidence_count":   len(self.iptmnet_data),
            "signor_residues":          "; ".join(s.get("residue", "") for s in self.signor_data if s.get("residue")),
            "signor_effects":           "; ".join(set(s.get("effect", "") for s in self.signor_data if s.get("effect"))),
            "omnipath_enzsub_sources":  _src_str(enzsub_with_res),
            "omnipath_interact_sources":_src_str(interactions),
            "uniprot_ptm_b_count":      len(self.uniprot_ptm_b),
            "reason":                   self.reason,
        }


class PhosphositeValidator:
    def __init__(self, omnipath: OmniPathClient, signor: SIGNORClient,
                 iptmnet: iPTMnetClient, uniprot: UniProtAccessionCache,
                 observed: Dict[str, Set[str]],
                 api_delay: float = 0.3,
                 skip_uniprot: bool = False,
                 skip_iptmnet: bool = False,
                 psp_client: Optional["PhosphoSitePlusClient"] = None):
        self.omnipath = omnipath
        self.signor = signor
        self.iptmnet = iptmnet
        self.uniprot = uniprot
        self.observed = observed
        self.api_delay = api_delay
        self.skip_uniprot = skip_uniprot
        self.skip_iptmnet = skip_iptmnet
        self.psp_client = psp_client  # v1.5: PSP 기능 주석 클라이언트

    def validate_edge(self, gene_a: str, gene_b: str) -> EdgeVerdict:
        v = EdgeVerdict(gene_a, gene_b)
        v.observed_sites_a = sorted(self.observed.get(gene_a, []))
        v.observed_sites_b = sorted(self.observed.get(gene_b, []))

        # ① OmniPath (enzsub + interactions 이중 조회)
        v.omnipath_sites = self.omnipath.get_enzyme_substrate(gene_a, gene_b, self.api_delay)

        # ② SIGNOR 4.0
        signor_all = self.signor.get_interaction(gene_a, gene_b, self.api_delay)
        v.signor_data = [r for r in signor_all if "phospho" in r.get("mechanism", "").lower()]

        # ③ iPTMnet
        if not self.skip_iptmnet:
            iptm_a = self.iptmnet.get_ptm_ppi(gene_a, self.api_delay)
            related = [r for r in iptm_a if gene_b.upper() in r.get("interactant_name", "").upper()
                       or gene_b.upper() in r.get("interactant_id", "").upper()]
            subs = self.iptmnet.get_substrates(gene_a, self.api_delay)
            b_subs = [r for r in subs if gene_b.upper() in r.get("substrate_name", "").upper()
                      or gene_b.upper() in r.get("substrate_id", "").upper()]
            v.iptmnet_data = related + b_subs

        # ④ UniProt PTM 부위 (기질 단백질 B의 검증된 인산화 잔기를 이론 부위로 사용)
        if not self.skip_uniprot:
            v.uniprot_ptm_a = self.uniprot.get_uniprot_ptm_sites(gene_a)
            v.uniprot_ptm_b = self.uniprot.get_uniprot_ptm_sites(gene_b)

        # ── 이론적 부위 수집 (우선순위 순) ──────────────────────────────────
        seen_sites: Set[str] = set()

        # OmniPath enzsub — 잔기 수준 (가장 직접적)
        for op in v.omnipath_sites:
            rt = op.get("residue_type", "")
            rp = op.get("residue_offset")
            if rt and rp:
                s = f"{rt}{rp}"
                if s not in seen_sites:
                    seen_sites.add(s)
                    v.theoretical_sites.append(s)

        # SIGNOR 잔기 필드
        for sig in v.signor_data:
            res = sig.get("residue", "").strip()
            if res and res not in seen_sites:
                seen_sites.add(res)
                v.theoretical_sites.append(res)

        # iPTMnet 기질 부위
        for it in v.iptmnet_data:
            site = it.get("site", "").strip()
            if site and site not in seen_sites:
                seen_sites.add(site)
                v.theoretical_sites.append(f"iPTM:{site}")

        # UniProt 기질(B)의 검증된 인산화 잔기 — OmniPath/SIGNOR에 잔기 없을 때 보완
        # 잔기 표기 변환: "Phosphoserine" @ pos → S{pos}
        _aa_map = {"Phosphoserine": "S", "Phosphothreonine": "T", "Phosphotyrosine": "Y"}
        for ptm in v.uniprot_ptm_b:
            pos = ptm.get("position")
            desc = ptm.get("description", "")
            prefix = _aa_map.get(desc, "")
            if pos and prefix:
                s = f"{prefix}{pos}"
                if s not in seen_sites:
                    seen_sites.add(s)
                    v.theoretical_sites.append(f"UniProt:{s}")

        self._match_sites(v)
        self._assign_grade(v)

        # ⑤ v1.5: PSP 기능 주석 — Grade A/B 매칭 부위에 효과/기능/다이나믹스 추가
        if self.psp_client and v.matched_sites:
            v.matched_sites_annotated = self.psp_client.annotate_matched_sites(
                v.gene_b, v.matched_sites, self.api_delay
            )
        else:
            v.matched_sites_annotated = list(v.matched_sites)

        return v

    def _match_sites(self, v: EdgeVerdict):
        """
        Two-tier matching for A→B phosphorylation edge:

        Tier 1 — Substrate match (primary, determines Grade A):
          Protein B observed phosphosites  vs  DB residues where A phosphorylates B
          Sources: OmniPath enzsub, SIGNOR residue, UniProt PTM of B
          → confirms "A phosphorylates B at this specific position"

        Tier 2 — Kinase activation match (supplementary):
          Protein A observed phosphosites  vs  UniProt validated PTM sites of A
          → confirms "A itself is phosphorylated (likely active)"
          Stored separately in v.matched_sites_a; does NOT trigger Grade A alone.

        Position comparison is always number-only regardless of S/T/Y prefix.
        """
        _aa_map = {"Phosphoserine": "S", "Phosphothreonine": "T", "Phosphotyrosine": "Y"}

        def _pos_map(site_list):
            """Build {position_number_str: original_label} map, prefixed label wins."""
            pm = {}
            for obs in site_list:
                pos = re.sub(r"[^0-9]", "", obs)
                if pos:
                    if pos not in pm or len(obs) > len(pm[pos]):
                        pm[pos] = obs
            return pm

        def _append_match(target_list, db_site, db_pos_str, obs_pos_map, source, n_refs):
            if db_pos_str in obs_pos_map:
                obs = obs_pos_map[db_pos_str]
                entry = {"site": db_site, "observed": obs,
                         "source": source, "n_references": n_refs}
                if not any(m["site"] == db_site and m["observed"] == obs
                           for m in target_list):
                    target_list.append(entry)

        # ── Tier 1: B observed sites vs DB substrate residues ─────────────
        obs_b_map = _pos_map(v.observed_sites_b)
        if obs_b_map:
            # OmniPath enzsub — residue on substrate B
            for op in v.omnipath_sites:
                if op.get("endpoint") != "enzsub":
                    continue
                rp = op.get("residue_offset")
                rt = op.get("residue_type", "")
                if not rp:
                    continue
                db_site = f"{rt}{rp}".strip() if rt else str(rp)
                srcs = op.get("sources", [])
                src_str = f"OmniPath({', '.join(srcs[:2])})" if srcs else "OmniPath"
                _append_match(v.matched_sites, db_site, str(rp),
                              obs_b_map, src_str, op.get("n_references", 0))

            # SIGNOR — residue on substrate B
            for sig in v.signor_data:
                res = sig.get("residue", "").strip()
                if not res:
                    continue
                db_pos = re.sub(r"[^0-9]", "", res)
                eff = sig.get("effect", "")
                src_str = f"SIGNOR({eff})" if eff else "SIGNOR"
                _append_match(v.matched_sites, res, db_pos,
                              obs_b_map, src_str, 1)

            # UniProt PTM of B
            for ptm in v.uniprot_ptm_b:
                pos    = ptm.get("position")
                desc   = ptm.get("description", "")
                eco    = ptm.get("evidence_code", "")
                prefix = _aa_map.get(desc, "")
                if pos and prefix:
                    db_site = f"{prefix}{pos}"
                    src_str = f"UniProt-B({eco})" if eco else "UniProt-B"
                    _append_match(v.matched_sites, db_site, str(pos),
                                  obs_b_map, src_str, 1)

        # ── Tier 2: A observed sites vs UniProt PTM sites of A ────────────
        # These confirm kinase A is itself phosphorylated (active state).
        # Stored in v.matched_sites_a; does NOT alone determine Grade A.
        obs_a_map = _pos_map(v.observed_sites_a)
        if obs_a_map:
            for ptm in v.uniprot_ptm_a:
                pos    = ptm.get("position")
                desc   = ptm.get("description", "")
                eco    = ptm.get("evidence_code", "")
                prefix = _aa_map.get(desc, "")
                if pos and prefix:
                    db_site = f"{prefix}{pos}"
                    src_str = f"UniProt-A({eco})" if eco else "UniProt-A"
                    _append_match(v.matched_sites_a, db_site, str(pos),
                                  obs_a_map, src_str, 1)

    def _assign_grade(self, v: EdgeVerdict):
        # enzsub 전용 레코드 (잔기 있는 것)
        enzsub_with_residue = [o for o in v.omnipath_sites
                               if o.get("residue_offset") and o.get("endpoint") == "enzsub"]
        # interactions 레코드 (잔기 없지만 방향 있음)
        interaction_records = [o for o in v.omnipath_sites
                               if o.get("endpoint") == "interactions"]
        # UniProt 이론 부위
        uniprot_theo = [s for s in v.theoretical_sites if s.startswith("UniProt:")]

        has_op_residue   = len(enzsub_with_residue) > 0
        has_op_interact  = len(interaction_records) > 0
        has_op_any       = len(v.omnipath_sites) > 0
        has_sig          = len(v.signor_data) > 0
        has_ipt          = len(v.iptmnet_data) > 0
        has_uniprot_theo = len(uniprot_theo) > 0
        has_match        = len(v.matched_sites) > 0       # Tier1: B substrate match
        has_match_a      = len(v.matched_sites_a) > 0    # Tier2: A kinase activation match

        score = 0.0
        reasons = []

        if has_match:
            # ── Grade A: observed site = DB residue (Tier 1: B substrate match) ─
            v.grade = EdgeVerdict.GRADE_A
            refs  = sum(m.get("n_references", 0) for m in v.matched_sites)
            # A kinase-activation match slightly boosts confidence
            a_bonus = 0.03 if has_match_a else 0.0
            score = min(1.0, 0.85 + len(v.matched_sites) * 0.05 + a_bonus)
            sites_str = ", ".join(m["site"] for m in v.matched_sites[:4])
            src_str   = ", ".join(set(m["source"].split("(")[0] for m in v.matched_sites))
            reasons.append(f"Observed site(s) ({sites_str}) match {src_str} validated residue(s) "
                           f"({refs} reference(s))")
            # Tier 2: report kinase A activation status
            if has_match_a:
                a_sites = ", ".join(f"{m['site']}(obs:{m['observed']})" for m in v.matched_sites_a[:3])
                reasons.append(f"Kinase {v.gene_a} activation confirmed: observed phosphosites "
                               f"match UniProt-A validated PTM residue(s) — {a_sites}")

        elif has_op_residue or has_sig:
            # ── Grade B-1: DB has residue-level phospho evidence, observed mismatch ──
            v.grade = EdgeVerdict.GRADE_B
            score = 0.55
            if has_op_residue:
                score += 0.10
                op_src = set(s for o in enzsub_with_residue for s in o.get("sources", []))
                theo_clean = [s for s in v.theoretical_sites if not s.startswith(("UniProt:", "iPTM:"))]
                reasons.append(f"OmniPath enzsub ({', '.join(list(op_src)[:3])}) records "
                                f"phosphorylation at residue(s) {', '.join(theo_clean[:4])}")
            if has_sig:
                score += 0.10
                sig_res = [s.get("residue") for s in v.signor_data if s.get("residue")]
                sig_eff = list(set(s.get("effect", "") for s in v.signor_data if s.get("effect")))
                reasons.append(f"SIGNOR: residue(s)={', '.join(sig_res[:3])}, effect(s)={', '.join(sig_eff)}")
            obs_all = v.observed_sites_a + v.observed_sites_b
            if obs_all:
                reasons.append(f"⚠ Observed site(s) ({', '.join(obs_all[:3])}) do not match "
                                f"DB functional residue position(s) — further validation recommended")
            else:
                reasons.append("⚠ No phosphosites observed in source data for either protein — "
                                "theoretical residues recorded only")

        elif has_op_interact or has_op_any:
            # ── Grade B-2: OmniPath directed interaction only (no residue) ──
            v.grade = EdgeVerdict.GRADE_B
            score = 0.55
            if has_op_interact:
                op_src = set(s for o in interaction_records for s in o.get("sources", []))
                stim_list = [o for o in interaction_records if o.get("consensus_stimulation") in ("1","True","true")]
                score += 0.05
                reasons.append(f"OmniPath interactions ({', '.join(list(op_src)[:3])}) records "
                                f"directed interaction (stimulatory={len(stim_list)}) — no residue-level data")
            else:
                op_src = set(s for o in v.omnipath_sites for s in o.get("sources", []))
                reasons.append(f"OmniPath ({', '.join(list(op_src)[:3])}) records interaction — residue not provided")
            if has_uniprot_theo:
                up_sites = [s.replace("UniProt:", "") for s in uniprot_theo[:4]]
                reasons.append(f"UniProt validated phosphosites on substrate {v.gene_b}: "
                                f"{', '.join(up_sites)} (direct mechanistic link unconfirmed)")

        elif has_ipt:
            # ── Grade B-3: iPTMnet PTM-dependent PPI only ────────────────
            v.grade = EdgeVerdict.GRADE_B
            score = 0.60
            for it in v.iptmnet_data[:2]:
                reasons.append(f"iPTMnet: {it.get('ptm_type')} @ {it.get('site')} → "
                                f"{it.get('association', 'N/A')}")

        elif has_uniprot_theo:
            # ── Grade C: UniProt phosphosites only (mechanism unconfirmed) ─
            v.grade = EdgeVerdict.GRADE_C
            score = 0.38
            up_sites = [s.replace("UniProt:", "") for s in uniprot_theo[:4]]
            reasons.append(f"UniProt records phosphorylated residue(s) {', '.join(up_sites)} "
                           f"on substrate {v.gene_b} — no direct evidence that {v.gene_a} "
                           f"phosphorylates these residue(s); non-phospho PPI or alternative kinase possible")

        else:
            # ── Grade D: no evidence in any DB ──────────────────────────
            rev = self.omnipath.get_enzyme_substrate(v.gene_b, v.gene_a, self.api_delay)
            if rev:
                v.grade = EdgeVerdict.GRADE_C
                score = 0.40
                reasons.append(f"Reverse-direction ({v.gene_b}→{v.gene_a}) phospho evidence found in OmniPath — "
                                f"forward phosphorylation unconfirmed; scaffold/complex association possible")
            else:
                v.grade = EdgeVerdict.GRADE_D
                score = 0.12
                reasons.append(f"No phosphorylation evidence for {v.gene_a}→{v.gene_b} "
                                f"in OmniPath, SIGNOR, iPTMnet, or UniProt — potential false positive")
                if not (v.observed_sites_a or v.observed_sites_b):
                    reasons.append("No phosphosites observed in source data for either protein — "
                                   "false-positive risk elevated")

        v.confidence_score = min(1.0, score)
        v.reason = " | ".join(reasons) or "No evidence found"


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 7: 경로 수준 신뢰도 계산
# ══════════════════════════════════════════════════════════════════════════════

def compute_pathway_scores(pathways, verdicts_map):
    results = []
    for i, path in enumerate(pathways, 1):
        scores, edge_details = [], []
        grade_counts = {"A": 0, "B": 0, "C": 0, "D": 0}
        for j in range(len(path) - 1):
            a, b = path[j], path[j + 1]
            if a == b:
                continue
            v = verdicts_map.get((a, b))
            if v:
                scores.append(v.confidence_score)
                edge_details.append(f"{a}→{b}[{v.grade[0]}:{v.confidence_score:.2f}]")
                grade_counts[v.grade[0]] += 1
        if not scores:
            continue
        geomean = math.exp(sum(math.log(max(s, 0.01)) for s in scores) / len(scores))
        pathway_score = geomean * 0.7 + min(scores) * 0.3
        results.append({
            "pathway_id": i,
            "pathway": " → ".join(path),
            "pathway_score": round(pathway_score, 3),
            "geomean": round(geomean, 3),
            "min_edge_score": round(min(scores), 3),
            "n_edges": len(scores),
            "grade_A": grade_counts["A"],
            "grade_B": grade_counts["B"],
            "grade_C": grade_counts["C"],
            "grade_D": grade_counts["D"],
            "edge_details": "; ".join(edge_details),
        })
    results.sort(key=lambda x: -x["pathway_score"])
    return results


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 8: 보고서 생성
# ══════════════════════════════════════════════════════════════════════════════

GRADE_COLORS = {"A": "#27ae60", "B": "#f39c12", "C": "#3498db", "D": "#e74c3c"}
GRADE_ICONS  = {"A": "✅", "B": "⚠️", "C": "🔵", "D": "❌"}


def _build_psp_function_html(v: "EdgeVerdict") -> str:
    """
    v1.5: Grade A/B 매칭 부위의 PSP 기능 주석 HTML 셀 내용 생성.

    표시 형식:
      [Y1068] 🟢 activation  →  GRB2 binding, MAPK cascade (23 refs)
      [T308]  🔴 inhibition  →  blocks nuclear translocation (5 refs)

    PSP 데이터 없으면 '—' 반환.
    """
    ann_sites = v.matched_sites_annotated if v.matched_sites_annotated else v.matched_sites
    if not ann_sites:
        return "—"

    parts = []
    for m in ann_sites[:4]:  # 최대 4개 표시
        site  = m.get("site", "?")
        eff   = m.get("psp_effect", "")
        fn    = m.get("psp_function", "")
        proc  = m.get("psp_process", "")
        nrefs = m.get("psp_n_refs", 0)
        dis   = m.get("psp_diseases", "")
        dom   = m.get("psp_domain", "")
        src   = m.get("psp_source", "")

        # PSP 데이터가 없으면 부위명만 표시
        if not eff:
            parts.append(f'<span style="color:#27ae60;font-weight:bold">{site}</span>')
            continue

        ec = PhosphoSitePlusClient.EFFECT_COLORS.get(eff, "#95a5a6")
        emoji = {"activation":"🟢","inhibition":"🔴","binding":"🔵",
                 "localization":"🟣","degradation":"🟠","stability":"🟤",
                 "unknown":"⚪"}.get(eff, "⚪")

        desc_parts = []
        if fn:   desc_parts.append(fn[:60])
        if proc and proc not in fn: desc_parts.append(proc[:40])
        desc = "; ".join(desc_parts) or eff
        ref_str = f' <span style="color:#aaa;font-size:.8em">({nrefs}refs)</span>' if nrefs else ""
        dis_str = (f' <span style="color:#e74c3c;font-size:.78em" title="Disease: {dis}">⚕</span>'
                   if dis else "")
        dom_str = (f' <span style="color:#9b59b6;font-size:.78em" title="Domain: {dom}">◈</span>'
                   if dom else "")
        tooltip = f"Source: {src} | Effect: {eff}" + (f" | Domain: {dom}" if dom else "") + \
                  (f" | Disease: {dis}" if dis else "")

        parts.append(
            f'<div style="margin-bottom:3px" title="{tooltip}">'
            f'<span style="color:#27ae60;font-weight:bold">[{site}]</span> '
            f'{emoji} <span style="color:{ec}">{eff}</span>'
            f'{dom_str}{dis_str}{ref_str}'
            f'<br><span style="color:#555;font-size:.82em">{desc}</span>'
            f'</div>'
        )

    return "\n".join(parts) if parts else "—"

def _score_to_color(score: float) -> str:
    """신뢰도 점수 → 화살표 색상 (연속 그라데이션)"""
    if   score >= 0.85: return "#1a9e50"   # 짙은 초록  (Grade A)
    elif score >= 0.65: return "#27ae60"   # 초록       (Grade B 강)
    elif score >= 0.50: return "#f39c12"   # 주황       (Grade B 약)
    elif score >= 0.35: return "#3498db"   # 파랑       (Grade C)
    elif score >= 0.20: return "#e67e22"   # 진한 주황  (D 경계)
    else:               return "#e74c3c"   # 빨강       (Grade D)


def _score_to_thick(score: float) -> int:
    """점수 → 화살표 두께 (px)"""
    return max(2, min(8, int(score * 9)))


def _build_pathway_section(
    pathways: List[List[str]],
    verdicts_map: Dict,
    pathway_scores: List[Dict],
) -> str:
    """
    Per-pathway colored-arrow visualization.
    Shows ALL pathways in their original input order (pathway #1, #2, ...).
    """
    # pathway_id → score record map
    ps_map = {p["pathway_id"]: p for p in pathway_scores}

    def _arrow_svg(score: float, grade: str, gene_a: str, gene_b: str,
                   matched: str, theo: str) -> str:
        color = _score_to_color(score)
        thick = _score_to_thick(score)
        w = 72
        tip = (f"{gene_a}→{gene_b} | Grade {grade} | score={score:.2f}"
               + (f" | matched:{matched}" if matched else "")
               + (f" | theoretical:{theo[:30]}" if theo else ""))
        tip_esc = tip.replace('"', '&quot;').replace("'", "&#39;")
        al = 12; x2 = w - al; mid_y = 14; aw = 10
        return (
            f'<div class="pedge" title="{tip_esc}">'
            f'<svg width="{w}" height="28" viewBox="0 0 {w} 28">'
            f'<line x1="0" y1="{mid_y}" x2="{x2}" y2="{mid_y}" '
            f'stroke="{color}" stroke-width="{thick}" stroke-linecap="round"/>'
            f'<polygon points="{x2},{mid_y-aw//2} {w},{mid_y} {x2},{mid_y+aw//2}" fill="{color}"/>'
            f'</svg>'
            f'<div class="edge-badge" style="background:{color}">{grade} {score:.2f}</div>'
            f'</div>'
        )

    def _grade_bar(ps: Dict) -> str:
        n = ps.get("n_edges", 0)
        if n == 0:
            return ""
        segs = []
        for g, color in [("grade_A","#27ae60"),("grade_B","#f39c12"),
                          ("grade_C","#3498db"),("grade_D","#e74c3c")]:
            cnt = ps.get(g, 0)
            if cnt > 0:
                segs.append(f'<div title="Grade {g[-1]}: {cnt}" '
                            f'style="background:{color};width:{cnt/n*100:.0f}%;height:100%"></div>')
        return (f'<div style="display:flex;height:6px;border-radius:3px;overflow:hidden;'
                f'width:120px;background:#ecf0f1;margin-top:3px">' + "".join(segs) + "</div>")

    # ── Build cards in ORIGINAL pathway order ──────────────────────────
    rows = []
    for pid, path in enumerate(pathways, 1):
        ps       = ps_map.get(pid, {})
        pw_score = ps.get("pathway_score", 0.0)
        pw_color = _score_to_color(pw_score)

        flow_parts = []
        for i, node in enumerate(path):
            flow_parts.append(f'<div class="pnode">{node}</div>')
            if i < len(path) - 1:
                nxt = path[i + 1]
                v   = verdicts_map.get((node, nxt))
                if v:
                    sc    = v.confidence_score
                    grade = v.grade[0]
                    mat   = ", ".join(m["site"] for m in v.matched_sites[:3])
                    theo  = ", ".join(s for s in v.theoretical_sites[:3]
                                      if not s.startswith(("UniProt:","iPTM:")))
                else:
                    sc, grade, mat, theo = 0.0, "?", "", ""
                flow_parts.append(_arrow_svg(sc, grade, node, nxt, mat, theo))

        n_edges = ps.get("n_edges", len(path) - 1)
        rows.append(
            f'<div class="pw-card" data-score="{pw_score:.3f}">'
            f'<div class="pw-header">'
            f'<span class="pw-badge" style="background:{pw_color}">Pathway #{pid}</span>'
            f'<span class="pw-score" style="color:{pw_color}">Score {pw_score:.3f}</span>'
            f'<span class="pw-meta">{len(path)} proteins · {n_edges} edges</span>'
            f'{_grade_bar(ps)}'
            f'</div>'
            f'<div class="pw-flow">' + "".join(flow_parts) + f'</div>'
            f'</div>'
        )

    n_total = len(rows)
    return f"""
<h2 style="margin-top:36px;border-bottom:2px solid #bdc3c7;padding-bottom:6px">
  &#128506; Pathway Confidence Visualization
  <span style="font-size:.7em;color:#7f8c8d;font-weight:normal">
    (Original input order &mdash; all {n_total} pathway(s) shown)
  </span>
</h2>

<div style="background:#fff;padding:12px 16px;border-radius:10px;
            box-shadow:0 2px 8px rgba(0,0,0,.06);margin-bottom:12px;font-size:.88em">
  <b>Arrow color legend:</b>&ensp;
  <span style="display:inline-flex;align-items:center;gap:4px;margin-right:12px">
    <svg width="32" height="12"><line x1="0" y1="6" x2="22" y2="6" stroke="#1a9e50" stroke-width="5"/>
    <polygon points="22,2 32,6 22,10" fill="#1a9e50"/></svg>
    <b style="color:#1a9e50">&ge;0.85</b> Grade A match
  </span>
  <span style="display:inline-flex;align-items:center;gap:4px;margin-right:12px">
    <svg width="32" height="12"><line x1="0" y1="6" x2="22" y2="6" stroke="#f39c12" stroke-width="4"/>
    <polygon points="22,2 32,6 22,10" fill="#f39c12"/></svg>
    <b style="color:#f39c12">0.50&ndash;0.84</b> Grade B
  </span>
  <span style="display:inline-flex;align-items:center;gap:4px;margin-right:12px">
    <svg width="32" height="12"><line x1="0" y1="6" x2="22" y2="6" stroke="#3498db" stroke-width="3"/>
    <polygon points="22,2 32,6 22,10" fill="#3498db"/></svg>
    <b style="color:#3498db">0.35&ndash;0.49</b> Grade C
  </span>
  <span style="display:inline-flex;align-items:center;gap:4px">
    <svg width="32" height="12"><line x1="0" y1="6" x2="22" y2="6" stroke="#e74c3c" stroke-width="2"/>
    <polygon points="22,2 32,6 22,10" fill="#e74c3c"/></svg>
    <b style="color:#e74c3c">&lt;0.35</b> Grade D
  </span>
  &ensp;|&ensp; Hover over an arrow to see edge details.
  &ensp;Arrow thickness is proportional to confidence score.
</div>

<div style="margin-bottom:12px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
  <label style="font-size:.88em">
    Min pathway score filter:
    <input type="range" id="scoreSlider" min="0" max="100" value="0"
           style="vertical-align:middle;width:140px"
           oninput="filterPathways(this.value/100)">
    <span id="sliderVal" style="color:#2980b9;font-weight:bold;
          min-width:34px;display:inline-block">0.00</span>
  </label>
  <button onclick="filterPathways(0)" style="padding:3px 10px;border-radius:4px;
    border:1px solid #bdc3c7;background:#f8f9fa;cursor:pointer;font-size:.85em">Show All</button>
  <span id="pwCount" style="font-size:.85em;color:#7f8c8d">{n_total} pathway(s) shown</span>
</div>

<div id="pathwayContainer">
{"".join(rows)}
</div>

<script>
function filterPathways(minScore) {{
  document.getElementById('sliderVal').textContent = parseFloat(minScore).toFixed(2);
  var cards = document.querySelectorAll('#pathwayContainer .pw-card');
  var visible = 0;
  cards.forEach(function(c) {{
    var s = parseFloat(c.dataset.score);
    if (s >= minScore) {{ c.style.display=''; visible++; }}
    else {{ c.style.display='none'; }}
  }});
  document.getElementById('pwCount').textContent = visible + ' pathway(s) shown';
}}
</script>
"""


def save_html(verdicts: List[EdgeVerdict], output_dir: Path, input_filename: str,
              pathways: Optional[List[List[str]]] = None,
              verdicts_map: Optional[Dict] = None,
              pathway_scores: Optional[List[Dict]] = None) -> str:
    now   = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    count = {g: sum(1 for v in verdicts if v.grade[0] == g) for g in "ABCD"}
    total = len(verdicts)

    def _src_str(sites):
        raw = set()
        for o in sites:
            for s in o.get("sources", []):
                if s: raw.add(s)
        return ", ".join(sorted(raw))[:60] or "—"

    # ── PPI edge detail table rows ──────────────────────────────────────
    rows_html = ""
    for v in sorted(verdicts, key=lambda x: -x.confidence_score):
        g     = v.grade[0]
        color = GRADE_COLORS.get(g, "#95a5a6")
        icon  = GRADE_ICONS.get(g, "")
        obs_a = ", ".join(v.observed_sites_a[:5]) or "—"
        obs_b = ", ".join(v.observed_sites_b[:5]) or "—"

        # Tier 1: B substrate match — primary evidence (determines Grade A)
        matched_b = ", ".join(
            f"{m['site']} (obs:{m['observed']})" for m in v.matched_sites
        ) or "—"

        # Tier 2: A kinase activation match — supplementary
        matched_a = ", ".join(
            f"<span title='{m['source']}'>{m['site']}</span>"
            for m in v.matched_sites_a
        ) or "—"
        matched_a_color = "#c0392b" if v.matched_sites_a else "#bbb"

        theo_db  = [s for s in v.theoretical_sites if not s.startswith(("UniProt:","iPTM:"))]
        theo_up  = [s.replace("UniProt:","") for s in v.theoretical_sites if s.startswith("UniProt:")]
        theo_ipt = [s.replace("iPTM:","") for s in v.theoretical_sites if s.startswith("iPTM:")]
        theo_parts = []
        if theo_db:  theo_parts.append(f'<span style="color:#8e44ad">{", ".join(theo_db[:4])}</span>')
        if theo_up:  theo_parts.append(f'<span style="color:#2980b9" title="UniProt validated">UP:{", ".join(theo_up[:3])}</span>')
        if theo_ipt: theo_parts.append(f'<span style="color:#27ae60" title="iPTMnet">iPTM:{", ".join(theo_ipt[:2])}</span>')
        theo_html = " | ".join(theo_parts) or "—"

        sig_res = ", ".join(s.get("residue","") for s in v.signor_data if s.get("residue")) or "—"
        sig_eff = ", ".join(set(s.get("effect","") for s in v.signor_data if s.get("effect"))) or "—"

        enzsub_sites  = [o for o in v.omnipath_sites if o.get("endpoint") == "enzsub"]
        enzsub_w_res  = [o for o in enzsub_sites if o.get("residue_offset")]
        interactions  = [o for o in v.omnipath_sites if o.get("endpoint") == "interactions"]

        if enzsub_w_res:
            res_labels = ", ".join(
                f"{o.get('residue_type','')}{o.get('residue_offset','')}"
                for o in enzsub_w_res[:4]
            )
            op_enzsub_str = f"{len(enzsub_w_res)} residue(s): {res_labels} ({_src_str(enzsub_w_res)})"
        elif enzsub_sites:
            op_enzsub_str = f"{len(enzsub_sites)} record(s), no residue ({_src_str(enzsub_sites)})"
        else:
            op_enzsub_str = "—"

        if interactions:
            stim = sum(1 for o in interactions if str(o.get("consensus_stimulation","")) in ("1","True","true"))
            op_inter_str = f"{len(interactions)} ({_src_str(interactions)})" + (f" stim:{stim}" if stim else "")
        else:
            op_inter_str = "—"

        # ── v1.5: PSP 기능 주석 HTML 생성 ──────────────────────────────────
        psp_html = _build_psp_function_html(v)

        rows_html += f"""
        <tr>
          <td><b>{v.gene_a}</b></td><td><b>{v.gene_b}</b></td>
          <td style="background:{color};color:white;text-align:center;border-radius:4px;white-space:nowrap">{icon} {g}</td>
          <td style="text-align:center">
            <div style="background:#ecf0f1;border-radius:8px;height:14px;overflow:hidden;margin-bottom:2px">
              <div style="background:{color};width:{v.confidence_score*100:.0f}%;height:100%"></div>
            </div>
            {v.confidence_score:.2f}
          </td>
          <td style="text-align:center">{v.pathway_appearances}</td>
          <td style="color:#2980b9;font-size:.85em">{obs_a}</td>
          <td style="color:#2980b9;font-size:.85em">{obs_b}</td>
          <td style="color:#27ae60;font-weight:bold;font-size:.85em">{matched_b}</td>
          <td style="color:{matched_a_color};font-size:.82em" title="Tier2: kinase A activation sites (UniProt-A)">{matched_a}</td>
          <td style="font-size:.82em">{psp_html}</td>
          <td style="font-size:.82em">{theo_html}</td>
          <td style="font-size:.82em">{sig_res}</td>
          <td style="color:#e67e22;font-size:.82em">{sig_eff}</td>
          <td style="font-size:.78em;color:#7f8c8d">{op_enzsub_str}</td>
          <td style="font-size:.78em;color:#999">{op_inter_str}</td>
          <td style="font-size:.78em;color:#555;max-width:240px">{v.reason[:240]}</td>
        </tr>"""

    # ── Pathway visualization section ──────────────────────────────────
    pathway_viz_html = ""
    if pathways and verdicts_map is not None and pathway_scores is not None:
        pathway_viz_html = _build_pathway_section(pathways, verdicts_map, pathway_scores)

    html = f"""<!DOCTYPE html><html lang="en">
<head><meta charset="UTF-8">
<title>PPI Phosphosite Validation &mdash; {input_filename}</title>
<style>
  body{{font-family:'Segoe UI',sans-serif;background:#f5f6fa;margin:0;padding:20px}}
  h1,h2{{color:#2c3e50}}
  h1{{border-bottom:3px solid #3498db;padding-bottom:8px}}
  .grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:18px 0}}
  .card{{background:white;border-radius:10px;padding:16px;
         box-shadow:0 2px 8px rgba(0,0,0,.08);text-align:center}}
  .card h3{{margin:0 0 6px;font-size:1.9em}}
  .card p{{margin:0;color:#7f8c8d;font-size:.88em}}
  .ga{{border-top:4px solid #27ae60}}.gb{{border-top:4px solid #f39c12}}
  .gc{{border-top:4px solid #3498db}}.gd{{border-top:4px solid #e74c3c}}
  table{{border-collapse:collapse;width:100%;background:white;border-radius:10px;
         overflow:hidden;box-shadow:0 2px 8px rgba(0,0,0,.08)}}
  th{{background:#2c3e50;color:white;padding:9px 7px;font-size:.82em;text-align:left}}
  td{{padding:7px;border-bottom:1px solid #ecf0f1;vertical-align:top}}
  tr:hover{{background:#f8f9fa}}
  .legend{{background:white;padding:14px;border-radius:10px;margin:14px 0;
            box-shadow:0 2px 8px rgba(0,0,0,.08);font-size:.9em;line-height:1.7}}
  .legend2{{background:#eaf4fb;padding:10px 14px;border-radius:8px;margin:8px 0;font-size:.85em}}
  /* Pathway cards */
  .pw-card{{background:white;border-radius:10px;margin-bottom:10px;
            box-shadow:0 2px 6px rgba(0,0,0,.07);overflow:hidden;transition:box-shadow .15s}}
  .pw-card:hover{{box-shadow:0 4px 16px rgba(0,0,0,.13)}}
  .pw-header{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;
              padding:8px 14px;background:#f8f9fa;border-bottom:1px solid #ecf0f1;font-size:.88em}}
  .pw-badge{{color:white;border-radius:12px;padding:2px 10px;font-weight:bold;
             font-size:.9em;white-space:nowrap}}
  .pw-score{{font-weight:bold;font-size:1em}}
  .pw-meta{{color:#95a5a6;font-size:.85em}}
  .pw-flow{{display:flex;align-items:center;flex-wrap:wrap;
            padding:12px 16px;gap:0;row-gap:10px;overflow-x:auto}}
  .pnode{{background:#ecf0f1;border:1.5px solid #bdc3c7;border-radius:6px;
          padding:5px 10px;font-weight:bold;font-size:.88em;white-space:nowrap;
          color:#2c3e50;box-shadow:0 1px 3px rgba(0,0,0,.08);flex-shrink:0}}
  .pedge{{display:flex;flex-direction:column;align-items:center;
          cursor:pointer;flex-shrink:0}}
  .pedge:hover .edge-badge{{opacity:1;transform:translateY(-2px)}}
  .edge-badge{{color:white;border-radius:8px;padding:1px 5px;font-size:.72em;
               font-weight:bold;white-space:nowrap;opacity:.88;
               transition:opacity .12s,transform .12s;line-height:1.4}}
</style>
</head>
<body>
<h1>&#128300; PPI Phosphosite Validation Report</h1>
<p>Generated: {now} &nbsp;|&nbsp; Input: {input_filename} &nbsp;|&nbsp; Edges analyzed: {total}</p>

<div class="legend">
  <b>Grade criteria:</b><br>
  <span style="color:#27ae60">&#10003; <b>A</b>: Observed phosphosite position matches DB validated residue &mdash;
    position-number comparison (S/T/Y prefix optional in source data)</span><br>
  <span style="color:#f39c12">&#9888; <b>B</b>: DB records phosphorylation evidence (residue-level or directional interaction) &mdash;
    observed site absent or positionally mismatched</span><br>
  <span style="color:#3498db">&#9679; <b>C</b>: Non-phospho PPI or UniProt phosphosites only (kinase–substrate link unconfirmed)</span><br>
  <span style="color:#e74c3c">&#10007; <b>D</b>: No phosphorylation evidence in any DB &mdash; potential false positive</span>
</div>
<div class="legend2">
  &#128204; <b>Theoretical site color key:</b>
  <span style="color:#8e44ad">&#9632; Purple: OmniPath enzsub / SIGNOR residue-level evidence (direct phosphorylation)</span> &nbsp;
  <span style="color:#2980b9">&#9632; Blue UP: UniProt validated phosphosite (kinase–substrate link unconfirmed)</span> &nbsp;
  <span style="color:#27ae60">&#9632; Green iPTM: iPTMnet PTM-dependent PPI</span>
  &nbsp;&nbsp;&#9677; <b>PSP effect icons:</b>
  🟢 activation &nbsp; 🔴 inhibition &nbsp; 🔵 binding &nbsp; 🟣 localization &nbsp; 🟠 degradation &nbsp; ⚪ unknown
  &nbsp;&nbsp; ⚕ disease-associated &nbsp; ◈ domain annotation
</div>

<div class="grid">
  <div class="card ga"><h3 style="color:#27ae60">{count['A']}</h3>
    <p>Grade A<br>Phosphosite matched</p></div>
  <div class="card gb"><h3 style="color:#f39c12">{count['B']}</h3>
    <p>Grade B<br>Phospho evidence</p></div>
  <div class="card gc"><h3 style="color:#3498db">{count['C']}</h3>
    <p>Grade C<br>Non-phospho / UniProt only</p></div>
  <div class="card gd"><h3 style="color:#e74c3c">{count['D']}</h3>
    <p>Grade D<br>No evidence</p></div>
</div>

<h2 style="margin-top:36px;border-bottom:2px solid #bdc3c7;padding-bottom:6px">
  &#128203; PPI Edge Detail Table
</h2>
<div style="overflow-x:auto">
<table>
  <thead><tr>
    <th>Protein A</th><th>Protein B</th><th>Grade</th><th>Confidence</th>
    <th>Pathway<br>Count</th>
    <th>Observed sites (A)</th><th>Observed sites (B)</th>
    <th title="Tier1: B observed sites matched to DB substrate residues — determines Grade A">DB-matched sites (B) &#10003;</th>
    <th title="Tier2: A observed sites matched to UniProt-A validated PTM residues — kinase activation evidence" style="color:#c0392b">Kinase A<br>activation sites</th>
    <th title="v1.5 NEW: PhosphoSitePlus site-specific biological function annotation" style="background:#e8f5e9;color:#1a6b3a">&#128200; PSP Site Function<br><span style="font-size:.78em;font-weight:normal">(effect · process · refs)</span></th>
    <th>Theoretical sites<br>(purple=direct / blue=UniProt)</th>
    <th>SIGNOR residue</th><th>SIGNOR effect</th>
    <th>OmniPath<br>enzsub (residue)</th>
    <th>OmniPath<br>interactions</th>
    <th>Verdict rationale</th>
  </tr></thead>
  <tbody>{rows_html}</tbody>
</table></div>

{pathway_viz_html}

<p style="color:#95a5a6;font-size:.82em;margin-top:18px">
  Data sources: OmniPath enzsub &middot; OmniPath interactions &middot;
  SIGNOR 4.0 &middot; iPTMnet &middot; UniProt REST API &middot;
  <b>PhosphoSitePlus / PhosphoELM (v1.5)</b><br>
  PPI Phosphosite Validator v1.5
</p>
</body></html>"""

    path = output_dir / "phosphosite_validation_report.html"
    path.write_text(html, encoding="utf-8")
    return str(path)


def save_tsv(verdicts: List[EdgeVerdict], output_dir: Path) -> str:
    path = output_dir / "phosphosite_validation_results.tsv"
    rows = [v.to_dict() for v in verdicts]
    if _HAS_PANDAS:
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False, encoding="utf-8-sig")
    else:
        if rows:
            header = list(rows[0].keys())
            with open(path, "w", encoding="utf-8-sig") as f:
                f.write("\t".join(header) + "\n")
                for row in rows:
                    f.write("\t".join(str(row.get(k, "")) for k in header) + "\n")
    return str(path)


def save_txt(verdicts: List[EdgeVerdict], output_dir: Path) -> str:
    total = len(verdicts)
    count = {g: sum(1 for v in verdicts if v.grade[0] == g) for g in "ABCD"}
    lines = [
        "=" * 70,
        "  PPI Phosphosite Validation — Text Summary",
        f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 70, "",
        f"  Total analyzed PPI edges : {total}",
        f"  Grade A (exact match)   : {count['A']} ({count['A']/max(total,1)*100:.1f}%)",
        f"  Grade B (theoretical)   : {count['B']} ({count['B']/max(total,1)*100:.1f}%)",
        f"  Grade C (non-phospho)   : {count['C']} ({count['C']/max(total,1)*100:.1f}%)",
        f"  Grade D (no evidence)   : {count['D']} ({count['D']/max(total,1)*100:.1f}%)", "",
    ]
    for grade, label in [("A", "✅ Grade A — observed ↔ DB match"), ("B", "⚠️  Grade B — theoretical site"), ("D", "❌ Grade D — false-positive risk")]:
        lines += ["─" * 70, f"  {label}", "─" * 70]
        for v in verdicts:
            if v.grade[0] == grade:
                extra = ""
                if grade == "A":
                    extra = "Matched sites: " + ", ".join(m["site"] for m in v.matched_sites)
                elif grade == "B":
                    extra = "Theoretical sites: " + ", ".join(v.theoretical_sites[:4])
                else:
                    extra = v.reason[:80]
                lines.append(f"  {v.gene_a:12s} → {v.gene_b:12s}  [score={v.confidence_score:.2f}] {extra}")
        lines.append("")

    path = output_dir / "phosphosite_validation_summary.txt"
    path.write_text("\n".join(lines), encoding="utf-8")
    return str(path)


def save_pathway_tsv(pathway_scores: List[Dict], output_dir: Path) -> str:
    path = output_dir / "pathway_confidence_scores.tsv"
    if pathway_scores:
        if _HAS_PANDAS:
            pd.DataFrame(pathway_scores).to_csv(path, sep="\t", index=False, encoding="utf-8-sig")
        else:
            header = list(pathway_scores[0].keys())
            with open(path, "w", encoding="utf-8-sig") as f:
                f.write("\t".join(header) + "\n")
                for row in pathway_scores:
                    f.write("\t".join(str(row.get(k, "")) for k in header) + "\n")
    return str(path)


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 9: Tkinter UI (v4 스타일 호환)
# ══════════════════════════════════════════════════════════════════════════════

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PPI Phosphosite Validator v1.5")
        self.geometry("1020x820")
        self.minsize(800, 620)

        # ── StringVar / IntVar ──────────────────────────────────────────────
        self.input_path  = tk.StringVar()
        self.output_dir  = tk.StringVar(value=os.getcwd())
        self.status      = tk.StringVar(value="Idle")

        # 실행 옵션
        self.opt_max_edges   = tk.IntVar(value=0)           # 0 = 전체
        self.opt_api_delay   = tk.DoubleVar(value=0.3)
        self.opt_skip_up     = tk.BooleanVar(value=False)   # UniProt PTM
        self.opt_skip_ipt    = tk.BooleanVar(value=False)   # iPTMnet
        self.opt_use_omni    = tk.BooleanVar(value=True)    # OmniPath
        self.opt_use_signor  = tk.BooleanVar(value=True)    # SIGNOR
        self.opt_use_psp     = tk.BooleanVar(value=True)    # v1.5: PhosphoSitePlus 기능 주석

        self._stop_requested = False
        self._build_ui()

    # ── UI 구성 ─────────────────────────────────────────────────────────────
    def _build_ui(self):
        frm = ttk.Frame(self, padding=12)
        frm.pack(fill="both", expand=True)

        row = 0

        # ── 입력 파일 ────────────────────────────────────────────────────────
        ttk.Label(frm, text="Input pathway file (.txt):").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.input_path, width=80).grid(
            row=row, column=1, sticky="we", padx=6)
        ttk.Button(frm, text="Browse", command=self.browse_input).grid(
            row=row, column=2, sticky="e")
        row += 1

        # ── 출력 폴더 ────────────────────────────────────────────────────────
        ttk.Label(frm, text="Output folder:").grid(row=row, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(frm, textvariable=self.output_dir, width=80).grid(
            row=row, column=1, sticky="we", padx=6, pady=(6, 0))
        ttk.Button(frm, text="Select", command=self.browse_output).grid(
            row=row, column=2, sticky="e", pady=(6, 0))
        row += 1

        # ── 구분선 ───────────────────────────────────────────────────────────
        ttk.Separator(frm).grid(row=row, column=0, columnspan=3, sticky="we", pady=8)
        row += 1

        # ── 실행 옵션 패널 ──────────────────────────────────────────────────
        opt_lf = ttk.LabelFrame(frm, text="Run options", padding=8)
        opt_lf.grid(row=row, column=0, columnspan=3, sticky="we", pady=(0, 4))

        # 행 1: 최대 엣지 수 / API 딜레이
        opt_r1 = ttk.Frame(opt_lf)
        opt_r1.pack(fill="x", pady=(0, 4))

        ttk.Label(opt_r1, text="Max edges:").pack(side="left")
        ttk.Spinbox(opt_r1, textvariable=self.opt_max_edges,
                    from_=0, to=9999, width=7).pack(side="left", padx=4)
        ttk.Label(opt_r1, text="(0 = all)", foreground="#888").pack(side="left", padx=(0, 16))

        ttk.Label(opt_r1, text="API delay (sec):").pack(side="left")
        ttk.Spinbox(opt_r1, textvariable=self.opt_api_delay,
                    from_=0.0, to=5.0, increment=0.1, width=6,
                    format="%.1f").pack(side="left", padx=4)
        ttk.Label(opt_r1, text="(avoid overload, recommended 0.3)", foreground="#888").pack(side="left")

        # 행 2: DB 사용 체크박스
        opt_r2 = ttk.Frame(opt_lf)
        opt_r2.pack(fill="x", pady=(0, 2))

        ttk.Label(opt_r2, text="Databases:").pack(side="left")
        ttk.Checkbutton(opt_r2, text="OmniPath", variable=self.opt_use_omni).pack(side="left", padx=8)
        ttk.Checkbutton(opt_r2, text="SIGNOR 4.0", variable=self.opt_use_signor).pack(side="left", padx=4)
        ttk.Checkbutton(opt_r2, text="iPTMnet",
                        variable=self.opt_skip_ipt,
                        onvalue=False, offvalue=True).pack(side="left", padx=4)
        ttk.Checkbutton(opt_r2, text="UniProt PTM",
                        variable=self.opt_skip_up,
                        onvalue=False, offvalue=True).pack(side="left", padx=4)

        # 행 2b: v1.5 PSP 기능 주석 옵션
        opt_r2b = ttk.Frame(opt_lf)
        opt_r2b.pack(fill="x", pady=(0, 4))
        ttk.Checkbutton(opt_r2b, text="📊 PhosphoSitePlus functional annotation (adds biological role and dynamics to Grade A/B matched sites)",
                        variable=self.opt_use_psp).pack(side="left", padx=24)

        # 행 3: 빠른 프리셋 버튼
        opt_r3 = ttk.Frame(opt_lf)
        opt_r3.pack(fill="x")

        ttk.Label(opt_r3, text="Presets:").pack(side="left")
        ttk.Button(opt_r3, text="🚀 Full (slow)",
                   command=lambda: self._apply_preset("full")).pack(side="left", padx=4)
        ttk.Button(opt_r3, text="⚡ Fast (OmniPath + SIGNOR only)",
                   command=lambda: self._apply_preset("fast")).pack(side="left", padx=4)
        ttk.Button(opt_r3, text="🔬 Test (50 edges)",
                   command=lambda: self._apply_preset("test")).pack(side="left", padx=4)
        ttk.Button(opt_r3, text="📵 Offline (cache only)",
                   command=lambda: self._apply_preset("offline")).pack(side="left", padx=4)

        row += 1

        # ── 구분선 ───────────────────────────────────────────────────────────
        ttk.Separator(frm).grid(row=row, column=0, columnspan=3, sticky="we", pady=6)
        row += 1

        # ── 실행/중지 버튼 + 상태 ────────────────────────────────────────────
        btns = ttk.Frame(frm)
        btns.grid(row=row, column=0, columnspan=3, sticky="we", pady=4)

        ttk.Button(btns, text="▶  Run", command=self.run_validation).pack(side="left")
        ttk.Button(btns, text="■  Stop (soft)", command=self.request_stop).pack(
            side="left", padx=8)
        ttk.Label(btns, text="Status:").pack(side="left", padx=(20, 4))
        ttk.Label(btns, textvariable=self.status, foreground="#2980b9",
                  font=("", 9, "bold")).pack(side="left")
        row += 1

        # ── 구분선 ───────────────────────────────────────────────────────────
        ttk.Separator(frm).grid(row=row, column=0, columnspan=3, sticky="we", pady=6)
        row += 1

        # ── 로그 ─────────────────────────────────────────────────────────────
        ttk.Label(frm, text="Log").grid(row=row, column=0, sticky="w")
        row += 1

        self.logbox = tk.Text(frm, height=24, wrap="word", font=("Consolas", 9))
        sb = ttk.Scrollbar(frm, command=self.logbox.yview)
        self.logbox.configure(yscrollcommand=sb.set)
        self.logbox.grid(row=row, column=0, columnspan=2, sticky="nsew")
        sb.grid(row=row, column=2, sticky="ns")

        frm.grid_columnconfigure(1, weight=1)
        frm.grid_rowconfigure(row, weight=1)

    # ── 프리셋 적용 ─────────────────────────────────────────────────────────
    def _apply_preset(self, preset: str):
        if preset == "full":
            self.opt_max_edges.set(0)
            self.opt_api_delay.set(0.3)
            self.opt_use_omni.set(True)
            self.opt_use_signor.set(True)
            self.opt_skip_ipt.set(False)
            self.opt_skip_up.set(False)
            self.opt_use_psp.set(True)
            self.log("[Preset] Full: all databases + PSP functional annotation enabled, validating all edges")
        elif preset == "fast":
            self.opt_max_edges.set(0)
            self.opt_api_delay.set(0.2)
            self.opt_use_omni.set(True)
            self.opt_use_signor.set(True)
            self.opt_skip_ipt.set(True)
            self.opt_skip_up.set(True)
            self.opt_use_psp.set(False)
            self.log("[Preset] Fast: OmniPath + SIGNOR only, skipping PSP/iPTMnet/UniProt")
        elif preset == "test":
            self.opt_max_edges.set(50)
            self.opt_api_delay.set(0.3)
            self.opt_use_omni.set(True)
            self.opt_use_signor.set(True)
            self.opt_skip_ipt.set(False)
            self.opt_skip_up.set(True)
            self.opt_use_psp.set(True)
            self.log("[Preset] Test: first 50 edges only, with PSP functional annotation")
        elif preset == "offline":
            self.opt_max_edges.set(0)
            self.opt_api_delay.set(0.0)
            self.opt_use_omni.set(False)
            self.opt_use_signor.set(False)
            self.opt_skip_ipt.set(True)
            self.opt_skip_up.set(True)
            self.opt_use_psp.set(False)
            self.log("[Preset] Offline: no API calls — using cached data only")

    # ── 파일/폴더 선택 ──────────────────────────────────────────────────────
    def browse_input(self):
        p = filedialog.askopenfilename(
            title="Select input pathway file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if p:
            self.input_path.set(p)

    def browse_output(self):
        d = filedialog.askdirectory(title="Select output folder")
        if d:
            self.output_dir.set(d)

    # ── 로그 ────────────────────────────────────────────────────────────────
    def log(self, msg: str):
        self.logbox.insert("end", f"[{time.strftime('%H:%M:%S')}] {msg}\n")
        self.logbox.see("end")
        self.update_idletasks()

    # ── 중지 요청 ────────────────────────────────────────────────────────────
    def request_stop(self):
        self._stop_requested = True
        self.log("Stop requested — the program will exit after the current edge is processed.")

    # ── 실행 ────────────────────────────────────────────────────────────────
    def run_validation(self):
        inp = self.input_path.get().strip()
        outdir = self.output_dir.get().strip()

        if not inp or not os.path.exists(inp):
            messagebox.showerror("Error", "Please select an input file.")
            return
        if not outdir:
            messagebox.showerror("Error", "Please select an output folder.")
            return

        if not _HAS_REQUESTS:
            messagebox.showwarning("Warning",
                "The requests package is not installed. API calls will be disabled.\n"
                "Only cached data will be used.\n\n"
                "Install with: pip install requests pandas tqdm")

        os.makedirs(outdir, exist_ok=True)
        self._stop_requested = False
        self.status.set("Running...")
        threading.Thread(target=self._worker, args=(inp, outdir), daemon=True).start()

    # ── 워커 스레드 ─────────────────────────────────────────────────────────
    def _worker(self, inp: str, outdir: str):
        try:
            self.log("=" * 55)
            self.log("  PPI Phosphosite Validator v1.5 started")
            self.log("=" * 55)

            # ─ 옵션 읽기 ────────────────────────────────────────────────────
            max_edges    = self.opt_max_edges.get()
            api_delay    = self.opt_api_delay.get()
            skip_uniprot = self.opt_skip_up.get()
            skip_iptmnet = self.opt_skip_ipt.get()
            use_omni     = self.opt_use_omni.get()
            use_signor   = self.opt_use_signor.get()
            use_psp      = self.opt_use_psp.get()   # v1.5

            self.log(f"[Options] max_edges={max_edges or 'all'} | api_delay={api_delay}s")
            self.log(f"[Options] OmniPath={use_omni} | SIGNOR={use_signor} | "
                     f"iPTMnet={'OFF' if skip_iptmnet else 'ON'} | "
                     f"UniProt={'OFF' if skip_uniprot else 'ON'} | "
                     f"PSP_functional_annotation={'ON' if use_psp else 'OFF'}")

            # ─ 파일 파싱 ───────────────────────────────────────────────────
            self.log(f"[Parsing] {os.path.basename(inp)}")
            fp = PathwayFileParser(inp).parse()
            observed = fp.observed_phosphosites
            pathways = fp.pathways
            edges    = list(fp.unique_edges)

            self.log(f"[Parsing] genes: {len(observed)} | phosphosites: {sum(len(v) for v in observed.values())}")
            self.log(f"[Parsing] pathways: {len(pathways)} | unique PPI edges: {len(edges)}")

            if not edges:
                messagebox.showwarning("Warning", "No edges were found. Please check the input file format.")
                self.status.set("Idle")
                return

            if max_edges and max_edges < len(edges):
                edges = edges[:max_edges]
                self.log(f"[Options] validating only the first {max_edges} edges.")

            # ─ 경로별 엣지 등장 횟수 ────────────────────────────────────────
            edge_count: Dict[Tuple, int] = defaultdict(int)
            for path in pathways:
                for i in range(len(path) - 1):
                    a, b = path[i], path[i + 1]
                    if a != b:
                        edge_count[(a, b)] += 1

            # ─ API 클라이언트 초기화 ─────────────────────────────────────────
            out = Path(outdir)
            uniprot_cache = UniProtAccessionCache(str(out / "uniprot_accession_cache.json"))
            omni_client   = OmniPathClient(str(out / "omnipath_cache.json"))
            sig_client    = SIGNORClient(str(out / "signor_cache.json"), uniprot_cache=uniprot_cache)
            ipt_client    = iPTMnetClient(uniprot_cache, str(out / "iptmnet_cache.json"))
            # v1.5: PhosphoSitePlus 기능 주석 클라이언트
            psp_client    = (PhosphoSitePlusClient(uniprot_cache, str(out / "psp_function_cache.json"))
                             if use_psp else None)

            # ─ API 연결 진단 테스트 ──────────────────────────────────────
            if _HAS_REQUESTS:
                self.log("[Diagnostics] testing API connectivity...")
                # OmniPath enzsub
                try:
                    _r = requests.get(OmniPathClient.ENZSUB_URL,
                                      params={"enzymes":"EGFR","substrates":"STAT3",
                                              "modification":"phosphorylation",
                                              "genesymbols":"yes","format":"json",
                                              "organisms":9606,
                                              "fields":"sources,references"},  # v16 FIX
                                      timeout=12)
                    _r.raise_for_status()
                    _d = _r.json()
                    _n_src = sum(1 for it in _d if isinstance(it, dict) and it.get("sources")) if isinstance(_d, list) else 0
                    self.log(f"  [OK] OmniPath enzsub connection successful — test returned "
                             f"{len(_d) if isinstance(_d,list) else '?'} records ({_n_src} with sources populated)")
                except Exception as _te:
                    _status = getattr(getattr(_te, "response", None), "status_code", "?")
                    self.log(f"  [FAIL] OmniPath enzsub connection failed: {type(_te).__name__}: {_te} (status={_status})")
                    self.log(f"         URL: {OmniPathClient.ENZSUB_URL}")
                # OmniPath interactions
                try:
                    _r2 = requests.get(OmniPathClient.INTERACTIONS_URL,
                                       params={"sources":"EGFR","targets":"STAT3",
                                               "genesymbols":"1","organisms":"9606",
                                               "fields":"sources,references",
                                               "format":"json"},
                                       timeout=12)
                    _r2.raise_for_status()
                    _d2 = _r2.json()
                    self.log(f"  [OK] OmniPath interactions connection successful — {len(_d2) if isinstance(_d2,list) else '?'} records")
                except Exception as _te2:
                    _status2 = getattr(getattr(_te2, "response", None), "status_code", "?")
                    self.log(f"  [FAIL] OmniPath interactions connection failed: {type(_te2).__name__}: {_te2} (status={_status2})")
                # SIGNOR
                try:
                    _r3 = requests.get(SIGNORClient.BASE_URL,
                                       params={"id": "P00533", "organism": "9606"},
                                       timeout=12)
                    _r3.raise_for_status()
                    _lines3 = [l for l in _r3.text.strip().split("\n") if l.strip()]
                    self.log(f"  [OK] SIGNOR getData.php connection successful — EGFR interactions: {max(0,len(_lines3)-1)} records")
                except Exception as _te3:
                    self.log(f"  [FAIL] SIGNOR getData.php connection failed: {type(_te3).__name__}: {_te3}")
                    self.log(f"         URL: {SIGNORClient.BASE_URL}?id=P00533&organism=9606")
                # v16 ADD: iPTMnet had NO connectivity test at all before this fix,
                # so a dead/renamed endpoint was indistinguishable from "no data".
                try:
                    _r6 = requests.get(f"{iPTMnetClient.BASE_URL}/P00533/ptmppi/", timeout=12)
                    _r6.raise_for_status()
                    _d6 = _r6.json()
                    self.log(f"  [OK] iPTMnet ptmppi connection successful — EGFR: "
                             f"{len(_d6) if isinstance(_d6, list) else '?'} records")
                except Exception as _te6:
                    _status6 = getattr(getattr(_te6, "response", None), "status_code", "?")
                    self.log(f"  [FAIL] iPTMnet ptmppi connection failed: {type(_te6).__name__}: {_te6} (status={_status6})")
                    self.log(f"         URL: {iPTMnetClient.BASE_URL}/P00533/ptmppi/")
                # v1.5: PSP 기능 주석 진단 (UniProt 주 소스 + PSP/PhosphoELM 보조)
                if use_psp:
                    self.log(f"  [INFO] Functional annotation strategy: UniProt (primary) + PSP/PhosphoELM (supporting)")
                    # PSP AJAX 보조 접근 가능성 확인
                    try:
                        _r4 = requests.get(
                            PhosphoSitePlusClient.PSP_RESIDUE_URL,
                            params={"uniprotAccs": "P00533", "types": "Phosphorylation"},
                            headers={"X-Requested-With": "XMLHttpRequest",
                                     "Referer": "https://www.phosphosite.org/"},
                            timeout=8)
                        _r4.raise_for_status()
                        _d4 = _r4.json()
                        _n4 = len(_d4) if isinstance(_d4, list) else len(_d4.get("sites", _d4.get("data", [])))
                        self.log(f"  [OK] PSP AJAX support endpoint reachable — EGFR sites: {_n4} records (n_refs/disease enrichment enabled)")
                    except Exception as _te4:
                        # PSP 실패는 정상 — UniProt이 주 소스
                        _phelm_ok = False
                        try:
                            _r5 = requests.get(
                                f"{PhosphoSitePlusClient.PHELM_BASE_URL}/P00533/",
                                headers={"Accept": "application/json"}, timeout=8)
                            _r5.raise_for_status()
                            _phelm_ok = True
                            self.log(f"  [OK] PhosphoELM support endpoint reachable — n_refs enrichment enabled")
                        except Exception:
                            pass
                        if not _phelm_ok:
                            self.log(f"  [INFO] PSP/PhosphoELM support endpoints unavailable — using UniProt-only functional annotation")
                            self.log(f"         (activation/inhibition effects and kinase information are still provided from UniProt)")
                self.log("[Diagnostics] complete")

            # 비활성 DB는 빈 결과 반환으로 패치
            if not use_omni:
                omni_client.get_enzyme_substrate = lambda a, b, d=0: []
            if not use_signor:
                sig_client.get_interaction = lambda a, b, d=0: []

            validator = PhosphositeValidator(
                omni_client, sig_client, ipt_client, uniprot_cache, observed,
                api_delay=api_delay,
                skip_uniprot=skip_uniprot,
                skip_iptmnet=skip_iptmnet,
                psp_client=psp_client,   # v1.5
            )

            # ─ 검증 루프 ────────────────────────────────────────────────────
            verdicts: List[EdgeVerdict] = []
            verdicts_map: Dict[Tuple, EdgeVerdict] = {}
            total = len(edges)

            self.log(f"[Validation] starting processing of {total} edges...")
            for idx, (gene_a, gene_b) in enumerate(edges, 1):
                if self._stop_requested:
                    self.log("Stopping as requested.")
                    break
                if idx % 10 == 0 or idx == 1:
                    self.log(f"  [{idx:4d}/{total}] {gene_a} → {gene_b}")
                else:
                    self.logbox.insert("end", f"  [{idx:4d}/{total}] {gene_a} → {gene_b}\n")
                    self.logbox.see("end")
                    self.update_idletasks()

                v = validator.validate_edge(gene_a, gene_b)
                v.pathway_appearances = edge_count.get((gene_a, gene_b), 0)
                verdicts.append(v)
                verdicts_map[(gene_a, gene_b)] = v

            # ─ API 오류 요약 출력 ───────────────────────────────────────────
            # v16 FIX: iPTMnet errors are now tracked too (previously this client
            # had no _api_errors list at all, so its failures were invisible).
            core_errs = (OmniPathClient._api_errors + SIGNORClient._api_errors
                         + iPTMnetClient._api_errors)
            # PSP AJAX 404는 예상된 동작 (UniProt으로 대체됨) — 별도 분류
            psp_errs_all  = PhosphoSitePlusClient._api_errors if use_psp else []
            psp_expected  = [e for e in psp_errs_all if "PSP AJAX" in e and "404" in e]
            psp_real_errs = [e for e in psp_errs_all if e not in psp_expected]
            all_api_errs  = core_errs + psp_real_errs
            if all_api_errs:
                self.log(f"[WARNING] API errors occurred: {len(all_api_errs)}")
                for _em in all_api_errs[:8]:
                    self.log(f"  ▶ {_em}")
                if len(all_api_errs) > 8:
                    self.log(f"  ... and {len(all_api_errs)-8} more errors")
                self.log("  → Check network connectivity, API endpoints, and firewall settings")
            else:
                self.log("[OK] No major API errors (OmniPath + SIGNOR + iPTMnet)")
            if psp_expected:
                self.log(f"[INFO] PSP AJAX unavailable for {len(psp_expected)} cases → automatically replaced by UniProt-based functional annotation")

            # v16 ADD: persist the FULL error list to a file, not just the first 8
            # lines in the on-screen log. This is what makes a persistently-empty
            # DB (e.g. OmniPath interactions or SIGNOR returning 0 for everything)
            # diagnosable after the run ends, instead of only visible while the
            # GUI window is still open.
            try:
                err_log_path = out / "api_errors.log"
                with open(err_log_path, "w", encoding="utf-8") as _ef:
                    _ef.write(f"Generated: {datetime.now().isoformat()}\n")
                    _ef.write(f"Total errors: {len(all_api_errs)}\n\n")
                    _ef.write(f"-- OmniPath ({len(OmniPathClient._api_errors)}) --\n")
                    _ef.write("\n".join(OmniPathClient._api_errors) + "\n\n")
                    _ef.write(f"-- SIGNOR ({len(SIGNORClient._api_errors)}) --\n")
                    _ef.write("\n".join(SIGNORClient._api_errors) + "\n\n")
                    _ef.write(f"-- iPTMnet ({len(iPTMnetClient._api_errors)}) --\n")
                    _ef.write("\n".join(iPTMnetClient._api_errors) + "\n\n")
                    if use_psp:
                        _ef.write(f"-- PhosphoSitePlus/PhosphoELM ({len(psp_real_errs)} real, "
                                  f"{len(psp_expected)} expected-404) --\n")
                        _ef.write("\n".join(psp_real_errs) + "\n")
                self.log(f"[Save] full API error log written to {err_log_path}")
            except Exception as _elog_e:
                self.log(f"[WARNING] could not write api_errors.log: {_elog_e}")

            # ─ 경로 수준 점수 계산 ──────────────────────────────────────────
            self.log("[Pathway score] calculating...")
            pathway_scores = compute_pathway_scores(pathways, verdicts_map)

            # ─ 보고서 저장 ──────────────────────────────────────────────────
            self.log("[Save] generating reports...")
            fname = os.path.basename(inp)
            html_path = save_html(verdicts, out, fname,
                                  pathways=pathways,
                                  verdicts_map=verdicts_map,
                                  pathway_scores=pathway_scores)
            tsv_path  = save_tsv(verdicts, out)
            txt_path  = save_txt(verdicts, out)
            pw_path   = save_pathway_tsv(pathway_scores, out)

            # ─ 최종 요약 ───────────────────────────────────────────────────
            count = {g: sum(1 for v in verdicts if v.grade[0] == g) for g in "ABCD"}
            n = len(verdicts)
            self.log("=" * 55)
            self.log("  Validation complete!")
            self.log(f"  Grade A (exact match)   : {count['A']:4d} ({count['A']/max(n,1)*100:.1f}%)")
            self.log(f"  Grade B (theoretical)   : {count['B']:4d} ({count['B']/max(n,1)*100:.1f}%)")
            self.log(f"  Grade C (non-phospho)   : {count['C']:4d} ({count['C']/max(n,1)*100:.1f}%)")
            self.log(f"  Grade D (no evidence)   : {count['D']:4d} ({count['D']/max(n,1)*100:.1f}%)")
            self.log("─" * 55)
            self.log(f"  HTML  : {html_path}")
            self.log(f"  TSV results: {tsv_path}")
            self.log(f"  TXT   : {txt_path}")
            self.log(f"  Pathway scores: {pw_path}")
            self.log("=" * 55)
            self.status.set("Completed")

        except Exception as e:
            self.status.set("Error")
            self.log(f"[Error] {e}")
            import traceback
            self.log(traceback.format_exc())


# ══════════════════════════════════════════════════════════════════════════════
# 섹션 10: 진입점
# ══════════════════════════════════════════════════════════════════════════════

def main():
    app = App()
    app.mainloop()

if __name__ == "__main__":
    main()