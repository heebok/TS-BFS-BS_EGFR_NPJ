# -*- coding: utf-8 -*-
"""
Fast_PPI_STRING - FINAL Target-Preserving Temporal BFS-BS + Progress UI (STRING API v12.5, Weighted + Regulatory Network)
- STRING v12.5 official pinned API (version-12-5.string-db.org)
- Edge weight = -ln(combined_score/1000)
- Regulatory network mode: supports directionality (A→B) and signal type (activation/inhibition)
- New genome analysis: supports direct input of protein amino-acid sequences
- Path ranking: based on cumulative cost
- Temporal v2: interval-by-interval sequential propagation; baseline is availability only
- Bridge Depth D counts intermediate bridge nodes (max edges per interval = D+1)
- Primary measured-anchor and strict signed-anchor modes
- Legacy v1 time_series backward compatibility retained
- FINAL target-preserving adaptive temporal beam: every reached target keeps at least one route; redundant histories are beam-pruned
- Temporal phosphosite output annotates measured endpoint sites without Cartesian expansion
- Cell-context bridge filter: bulk RNA/proteomics, bridge-only, experimental anchors protected, UNKNOWN retained by default\n- Cell-model identity cross-search: exact name/SIDM/Cellosaurus/DepMap aliases; no similar-name substitution\n- HEK293T aliases include HEK-293T, 293T, 293-T, 293 T, 293tsA1609neo and CVCL_0063; no HEK293 substitution
- First interval: baseline-available proteins may serve as one bridge after two direct STRING-threshold tests; baseline is not activation
- STRING partner retrieval explicitly requests a large threshold-qualified neighborhood; cache is score-aware
- Speed optimization: deterministic neighbor preprocessing, baseline direct-edge lookup, and exact segment memoization are cached
"""

import os, re, json, time, threading, traceback, math, itertools, csv, gzip, zipfile, io, sys, subprocess, html, hashlib
from datetime import datetime, timedelta
from collections import deque, Counter, OrderedDict
from typing import Dict, List, Optional, Set, Tuple
from threading import Lock
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from html.parser import HTMLParser
from urllib.parse import urljoin, urlparse

import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk


# ========== Utility Functions ==========

def canonicalize_isoform(acc: Optional[str]) -> Optional[str]:
    if not acc: return acc
    return re.sub(r"-\d+$", "", acc)


def looks_like_accession(s: Optional[str]) -> bool:
    if not s: return False
    return bool(
        re.fullmatch(r"[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-\d+)?", s) or
        re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9](-\d+)?", s)
    )


def _as_name(x) -> str:
    if isinstance(x, str):
        return x
    if isinstance(x, dict):
        for k in ("name", "gene", "symbol", "accession", "id"):
            v = x.get(k)
            if isinstance(v, str) and v.strip():
                return v.strip()
        return json.dumps(x, ensure_ascii=False)
    return str(x)


def parse_protein_with_site(protein_str: str) -> Tuple[str, Set[str]]:
    """Split site information from a protein name"""
    clean_str = protein_str.lstrip('-')
    if '@' in clean_str:
        parts = clean_str.split('@')
        name = parts[0]
        site_str = parts[1] if len(parts) > 1 else ''
        sites = set(site_str.split('_')) if site_str else set()
        return name, sites
    else:
        return clean_str, set()


def is_fasta_sequence(text: str) -> bool:
    """Check whether the input is in FASTA format"""
    text = text.strip()
    if text.startswith('>'):
        return True
    # Check whether it is a pure amino-acid sequence
    clean = re.sub(r'\s+', '', text)
    return bool(re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY*]+', clean, re.IGNORECASE))


def parse_fasta_sequences(text: str) -> Dict[str, str]:
    """Parse FASTA format -> {name: sequence}"""
    sequences = {}
    current_name = None
    current_seq = []

    for line in text.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_name and current_seq:
                sequences[current_name] = ''.join(current_seq)
            current_name = line[1:].split()[0]  # first word only
            current_seq = []
        elif line:
            current_seq.append(re.sub(r'\s+', '', line))

    if current_name and current_seq:
        sequences[current_name] = ''.join(current_seq)

    # single sequence without FASTA header
    if not sequences:
        clean = re.sub(r'\s+', '', text.strip())
        if re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY*]+', clean, re.IGNORECASE):
            sequences["UserSequence"] = clean

    return sequences


# ========== Persistent Cache ==========

class PersistentCache:
    def __init__(self, cache_file: str = "string_cache_v12_5_fullpartners_cell.json", ttl_days: int = 7):
        self.cache_file = cache_file
        self.ttl = timedelta(days=ttl_days)
        self.lock = Lock()

        self.name_to_gene: Dict[str, str] = {}
        self.gene_to_neighbors: Dict[str, List[Dict]] = {}
        self.timestamps: Dict[str, str] = {}

        self.load()

    def load(self):
        if not os.path.exists(self.cache_file):
            return
        try:
            with open(self.cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.name_to_gene = data.get("name_to_gene", {})
                self.gene_to_neighbors = data.get("gene_to_neighbors", {})
                self.timestamps = data.get("timestamps", {})
            print(f"[Cache] Loaded STRING cache from {self.cache_file}")
        except Exception as e:
            print(f"[Cache] Failed to load: {e}")

    def save(self):
        try:
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump({
                    "name_to_gene": self.name_to_gene,
                    "gene_to_neighbors": self.gene_to_neighbors,
                    "timestamps": self.timestamps
                }, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"[Cache] Failed to save: {e}")

    def is_expired(self, key: str) -> bool:
        if key not in self.timestamps:
            return True
        try:
            cached_time = datetime.fromisoformat(self.timestamps[key])
            return datetime.now() - cached_time > self.ttl
        except:
            return True

    def _update_timestamp(self, key: str):
        self.timestamps[key] = datetime.now().isoformat()

    def get_gene(self, name: str) -> Optional[str]:
        with self.lock:
            if name in self.name_to_gene and not self.is_expired(f"gene_{name}"):
                return self.name_to_gene[name]
            return None

    def put_gene(self, name: str, gene: str):
        if not (isinstance(name, str) and isinstance(gene, str) and name and gene):
            return
        with self.lock:
            self.name_to_gene[name] = gene
            self._update_timestamp(f"gene_{name}")

    def get_neighbors(self, gene: str, network_type: str = "functional",
                      min_score: int = 700) -> Optional[List[Dict]]:
        # STRING threshold is part of the cache identity so a cache produced for
        # another cutoff cannot silently alter the search space.
        key = f"v12.5_official_{gene}_{network_type}_score{int(min_score)}"
        with self.lock:
            if key in self.gene_to_neighbors and not self.is_expired(f"neigh_{key}"):
                return self.gene_to_neighbors[key]
            return None

    def put_neighbors(self, gene: str, neighbors: List[Dict], network_type: str = "functional",
                      min_score: int = 700):
        if not isinstance(gene, str) or not isinstance(neighbors, list):
            return
        key = f"v12.5_official_{gene}_{network_type}_score{int(min_score)}"
        with self.lock:
            self.gene_to_neighbors[key] = neighbors
            self._update_timestamp(f"neigh_{key}")


# ========== STRING API v12.5 ==========

# v12.0 stable URL
STRING_API_BASE_URL = "https://version-12-5.string-db.org/api"
STRING_API_VERSIONED_URL = STRING_API_BASE_URL
# v12.5 PREVIEW URL — as of 2026-09, this is the ONLY endpoint that actually
# returns directional regulatory data (source/target/sign fields). The
# default/pinned v12.0 domains above accept network_type=regulatory without
# error but silently return functional-network data with no directionality
# (verified empirically). This is a pre-release preview; its address/schema
# may change once STRING officially ships v12.5.
STRING_API_PREVIEW_URL = STRING_API_BASE_URL


class StringAPI:
    """STRING v12.5 API with edge weighting + regulatory network support"""

    def __init__(self, logger_print=None, cache_file: str = "string_cache_v12_5_fullpartners_cell.json",
                 network_type: str = "functional", use_versioned_url: bool = False):
        self.sess = self._create_session()
        self.cache = PersistentCache(cache_file=cache_file)
        self._log = logger_print
        if network_type not in ("functional", "physical", "regulatory"):
            raise ValueError("Supported STRING 12.5 networks: functional, physical, regulatory")
        self.network_type = network_type
        self._api_network_type = network_type
        self.base_url = STRING_API_BASE_URL  # Pinned: never silently use latest/preview/12.0
        self.snapshot_neighbors = {}  # all queried neighbors, not just final selected edges
        self.snapshot_lock = Lock()
        if logger_print:
            logger_print(f"[STRING] Official v12.5 pinned: {self.base_url} | {network_type}")

        # user-submitted genome information
        self.custom_species_id = None  # STRG0AXXXXX format
        self.sequence_mode = False

        self.stats = {
            "api_calls": 0,
            "cache_hits": 0,
            "cache_misses": 0,
            "total_time": 0.0,
            "avg_score": 0.0,
            "edges_retrieved": 0,
            "regulatory_edges": 0,
            "activation_edges": 0,
            "inhibition_edges": 0,
            "directional_edges": 0
        }

    def _create_session(self) -> requests.Session:
        session = requests.Session()
        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET", "POST"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        session.mount("http://", adapter)
        return session

    def log(self, msg: str):
        if self._log:
            self._log(msg)

    def _get_json(self, url: str, params: Optional[dict] = None) -> Optional[List]:
        self.stats["api_calls"] += 1
        start_time = time.time()

        try:
            r = self.sess.get(url, params=params, timeout=60)

            if r.status_code == 429:
                retry_after = int(r.headers.get('Retry-After', 60))
                retry_after = max(1, min(retry_after, 120))
                depth = int(getattr(self, "_manual_get_retry_depth", 0))
                if depth >= 4:
                    self.log("[STRING API] Rate limit retry limit reached (5 total attempts).")
                    return None
                self.log(
                    f"[STRING API] Rate limited. "
                    f"Retry {depth + 2}/5 after {retry_after}s..."
                )
                self._manual_get_retry_depth = depth + 1
                try:
                    time.sleep(retry_after)
                    return self._get_json(url, params)
                finally:
                    self._manual_get_retry_depth = depth

            if r.status_code == 503:
                depth = int(getattr(self, "_manual_get_retry_depth", 0))
                if depth >= 4:
                    self.log("[STRING API] 503 retry limit reached (5 total attempts).")
                    return None
                self.log(
                    f"[STRING API] Service temporarily unavailable. "
                    f"Retry {depth + 2}/5 after 30s..."
                )
                self._manual_get_retry_depth = depth + 1
                try:
                    time.sleep(30)
                    return self._get_json(url, params)
                finally:
                    self._manual_get_retry_depth = depth

            if r.status_code != 200:
                self.log(f"[STRING API] Error {r.status_code}: {r.text[:200]}")
                return None

            return r.json()

        except requests.Timeout:
            self.log(f"[STRING API] Timeout after 60s")
            return None
        except Exception as e:
            self.log(f"[STRING API] Error: {e}")
            return None
        finally:
            elapsed = time.time() - start_time
            self.stats["total_time"] += elapsed
            if elapsed > 10:
                self.log(f"[STRING API] Slow request: {elapsed:.1f}s")

    def _post_json(self, url: str, data: Optional[dict] = None) -> Optional[List]:
        """POST request (used for sequence submission, etc.)"""
        self.stats["api_calls"] += 1
        start_time = time.time()

        try:
            r = self.sess.post(url, data=data, timeout=120)

            if r.status_code == 429:
                retry_after = int(r.headers.get('Retry-After', 60))
                retry_after = max(1, min(retry_after, 120))
                depth = int(getattr(self, "_manual_post_retry_depth", 0))
                if depth >= 4:
                    self.log("[STRING API POST] Rate limit retry limit reached (5 total attempts).")
                    return None
                self.log(
                    f"[STRING API POST] Rate limited. "
                    f"Retry {depth + 2}/5 after {retry_after}s..."
                )
                self._manual_post_retry_depth = depth + 1
                try:
                    time.sleep(retry_after)
                    return self._post_json(url, data)
                finally:
                    self._manual_post_retry_depth = depth

            if r.status_code != 200:
                self.log(f"[STRING API POST] Error {r.status_code}: {r.text[:200]}")
                return None

            content_type = r.headers.get('Content-Type', '')
            if 'json' in content_type:
                return r.json()
            else:
                # handle non-JSON formats such as TSV as well
                return [{"raw_response": r.text}]

        except Exception as e:
            self.log(f"[STRING API POST] Error: {e}")
            return None
        finally:
            elapsed = time.time() - start_time
            self.stats["total_time"] += elapsed

    # ─── New Genome Sequence Submission (v12.0+) ────────────────────────────

    def submit_sequences(self, sequences: Dict[str, str],
                         blast_identity_threshold: float = 30.0,
                         blast_evalue_threshold: float = 1e-10) -> Optional[str]:
        """
        Map protein sequences to STRING for PPI analysis (3-step strategy)

        Strategy:
          Step 1: Name-based mapping (STRING get_string_ids)
          Step 2: If that fails, search for the most similar protein with NCBI BLAST
          Step 3: Remap the BLAST hit back to STRING

        Args:
            sequences: {protein_name: amino_acid_sequence}

        Returns:
            Summary string of mapping results or None
        """
        if not sequences:
            self.log("[Sequence] No sequences to submit.")
            return None

        self.log(f"[Sequence] {len(sequences)} protein sequences: processing started...")
        self.log(f"[Sequence] Strategy: 1) name mapping -> 2) BLAST homology search -> 3) STRING remapping")

        mapped_by_name = 0
        mapped_by_blast = 0
        failed = []
        mapping_details = {}  # {input_name: {method, gene, blast_hit, identity, evalue}}

        # -- Step 1: Name-based mapping --
        self.log(f"\n[Step 1] Trying name-based STRING mapping...")
        unmapped_seqs = {}

        for name, seq in sequences.items():
            url = f"{self.base_url}/json/get_string_ids"
            params = {
                "identifiers": name,
                "species": self._get_species(),
                "limit": 1
            }
            data = self._get_json(url, params)
            if data and len(data) > 0:
                gene = data[0].get("preferredName", "")
                if gene:
                    self.cache.put_gene(name, gene)
                    mapped_by_name += 1
                    mapping_details[name] = {
                        "method": "name",
                        "gene": gene,
                        "input_name": name
                    }
                    self.log(f"  ✓ {name} → {gene} (name mapping)")
                else:
                    unmapped_seqs[name] = seq
            else:
                unmapped_seqs[name] = seq
            time.sleep(1)

        self.log(f"[Step 1 Result] Name mapping success: {mapped_by_name}/{len(sequences)}")

        # -- Step 2: Run NCBI BLAST for unmapped sequences --
        if unmapped_seqs:
            self.log(f"\n[Step 2] Unmapped {len(unmapped_seqs)} sequences: running NCBI BLAST search...")
            self.log(f"  (Using the NCBI BLAST REST API; may take 1-3 minutes per sequence)")

            for name, seq in unmapped_seqs.items():
                if self.sequence_mode and hasattr(self, '_stop_flag') and self._stop_flag:
                    break

                blast_result = self._blast_sequence(
                    name, seq,
                    evalue_threshold=blast_evalue_threshold,
                    identity_threshold=blast_identity_threshold
                )

                if blast_result:
                    hit_name = blast_result["hit_name"]
                    identity = blast_result["identity"]
                    evalue = blast_result["evalue"]
                    accession = blast_result.get("accession", "")

                    self.log(f"  [BLAST] {name} -> Best hit: {hit_name} "
                             f"(identity={identity:.1f}%, E={evalue:.1e})")

                    # Step 3: Map the BLAST hit back to STRING
                    # First try accession, and if that fails, try hit_name
                    gene = None
                    for query_id in [accession, hit_name]:
                        if not query_id:
                            continue
                        url = f"{self.base_url}/json/get_string_ids"
                        params = {
                            "identifiers": query_id,
                            "species": self._get_species(),
                            "limit": 1
                        }
                        data = self._get_json(url, params)
                        if data and len(data) > 0:
                            gene = data[0].get("preferredName", "")
                            if gene:
                                break
                        time.sleep(1)

                    if gene:
                        self.cache.put_gene(name, gene)
                        mapped_by_blast += 1
                        mapping_details[name] = {
                            "method": "blast",
                            "gene": gene,
                            "input_name": name,
                            "blast_hit": hit_name,
                            "accession": accession,
                            "identity": identity,
                            "evalue": evalue
                        }
                        self.log(f"  ✓ {name} → {hit_name} → {gene} "
                                 f"(BLAST {identity:.1f}% identity)")
                    else:
                        failed.append(name)
                        self.log(f"  ✗ {name}: BLAST hit ({hit_name}) could not be mapped to STRING")
                else:
                    failed.append(name)
                    self.log(f"  ✗ {name}: BLAST search failed or no significant hit was found")

        # ── Result Summary ──
        total_mapped = mapped_by_name + mapped_by_blast
        self.log(f"\n{'=' * 50}")
        self.log(f"[Sequence Mapping Summary]")
        self.log(f"  Name-based mapping: {mapped_by_name}")
        self.log(f"  BLAST-based mapping: {mapped_by_blast}")
        self.log(f"  Mapping failed: {len(failed)}")
        self.log(f"  Total mapped successfully: {total_mapped}/{len(sequences)}")

        if mapping_details:
            self.log(f"\n[Mapping Details]")
            for name, info in mapping_details.items():
                if info["method"] == "name":
                    self.log(f"  {name} → {info['gene']} (name)")
                else:
                    self.log(f"  {name} → {info['blast_hit']} → {info['gene']} "
                             f"(BLAST {info['identity']:.1f}%, E={info['evalue']:.1e})")

        if failed:
            self.log(f"\n[Proteins that failed to map]")
            for name in failed:
                self.log(f"  ✗ {name}")
            self.log(f"\n  ⚠️  Failed proteins will be excluded from the analysis.")
            self.log(f"  Alternative: upload the full proteome to STRING Web (https://string-db.org) and then")
            self.log(f"       enter the assigned Species ID (STRG0AXXXXX) into the program.")

        self.log(f"{'=' * 50}")

        # store mapping details in the instance (for writing to result files)
        self.sequence_mapping_details = mapping_details

        self.sequence_mode = True
        return f"mapped_{total_mapped}_of_{len(sequences)}"

    def _blast_sequence(self, name: str, sequence: str,
                        evalue_threshold: float = 1e-10,
                        identity_threshold: float = 30.0) -> Optional[Dict]:
        """
        Search protein sequences using the NCBI BLAST REST API

        Args:
            name: protein name (for logging)
            sequence: amino-acid sequence
            evalue_threshold: E-value threshold (default 1e-10)
            identity_threshold: minimum identity % (default 30%)

        Returns:
            {"hit_name", "accession", "identity", "evalue", "score"} or None
        """
        self.log(f"    [BLAST] {name} ({len(sequence)}aa)  search in progress...")

        try:
            # ── Submit BLAST job ──
            put_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            put_params = {
                "CMD": "Put",
                "PROGRAM": "blastp",
                "DATABASE": "swissprot",  # curated DB → fast and highly reliable
                "QUERY": sequence,
                "FORMAT_TYPE": "JSON2",
                "HITLIST_SIZE": 5,
                "EXPECT": str(evalue_threshold),
                "FILTER": "L",  # low complexity filter
            }

            put_resp = self.sess.post(put_url, data=put_params, timeout=30)
            if put_resp.status_code != 200:
                self.log(f"    [BLAST] Submission failed: HTTP {put_resp.status_code}")
                return None

            # Extract RID
            rid_match = re.search(r"RID = (\S+)", put_resp.text)
            if not rid_match:
                self.log(f"    [BLAST] RID was not received")
                return None
            rid = rid_match.group(1)
            self.log(f"    [BLAST] Job submitted (RID: {rid})")

            # ── Waiting for results (up to 3 minutes) ──
            get_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            max_wait = 180  #  sec
            poll_interval = 10
            waited = 0

            while waited < max_wait:
                time.sleep(poll_interval)
                waited += poll_interval

                status_params = {
                    "CMD": "Get",
                    "RID": rid,
                    "FORMAT_OBJECT": "SearchInfo"
                }
                status_resp = self.sess.get(get_url, params=status_params, timeout=30)

                if "Status=WAITING" in status_resp.text:
                    if waited % 30 == 0:
                        self.log(f"    [BLAST] Waiting... ({waited}s elapsed)")
                    continue
                elif "Status=FAILED" in status_resp.text:
                    self.log(f"    [BLAST] Search failed")
                    return None
                elif "Status=READY" in status_resp.text:
                    break
                else:
                    continue

            if waited >= max_wait:
                self.log(f"    [BLAST] Timeout ({max_wait} sec)")
                return None

            # ── Retrieve results ──
            result_params = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
                "HITLIST_SIZE": 5
            }
            result_resp = self.sess.get(get_url, params=result_params, timeout=60)

            if result_resp.status_code != 200:
                self.log(f"    [BLAST] Failed to retrieve results")
                return None

            # Simple XML parsing (extract best hit)
            text = result_resp.text

            # Hit extract information
            hit_def_match = re.search(r"<Hit_def>(.+?)</Hit_def>", text)
            hit_acc_match = re.search(r"<Hit_accession>(.+?)</Hit_accession>", text)
            identity_match = re.search(r"<Hsp_identity>(\d+)</Hsp_identity>", text)
            align_len_match = re.search(r"<Hsp_align-len>(\d+)</Hsp_align-len>", text)
            evalue_match = re.search(r"<Hsp_evalue>(.+?)</Hsp_evalue>", text)
            score_match = re.search(r"<Hsp_bit-score>(.+?)</Hsp_bit-score>", text)

            if not (hit_def_match and identity_match and align_len_match and evalue_match):
                # when no hit is present
                if "<Iteration_message>No hits found</Iteration_message>" in text:
                    self.log(f"    [BLAST] No significant hit found")
                else:
                    self.log(f"    [BLAST] Failed to parse results")
                return None

            hit_def = hit_def_match.group(1)
            accession = hit_acc_match.group(1) if hit_acc_match else ""
            identity_count = int(identity_match.group(1))
            align_len = int(align_len_match.group(1))
            evalue = float(evalue_match.group(1))
            bit_score = float(score_match.group(1)) if score_match else 0

            identity_pct = (identity_count / align_len * 100) if align_len > 0 else 0

            # Try to extract the gene name from Hit_def
            # SwissProt format: "RecName: Full=xxx; AltName: ..." or "gene_SPECIES"
            hit_name = hit_def.split(";")[0]
            # From the GN= tag, try to extract the gene symbol
            gn_match = re.search(r"GN=(\S+)", text)
            if gn_match:
                hit_name = gn_match.group(1)

            if identity_pct < identity_threshold:
                self.log(f"    [BLAST] Best-hit identity({identity_pct:.1f}%) is below "
                         f"threshold ({identity_threshold}%)")
                return None

            return {
                "hit_name": hit_name,
                "accession": accession,
                "identity": identity_pct,
                "evalue": evalue,
                "score": bit_score,
                "hit_def": hit_def
            }

        except Exception as e:
            self.log(f"    [BLAST] Error: {e}")
            return None

    def set_custom_species(self, species_id: str):
        """Set the species ID for the user-submitted genome (STRG0AXXXXX format)"""
        if species_id and species_id.startswith("STRG"):
            self.custom_species_id = species_id
            self.log(f"[Species] Custom species ID set: {species_id}")
        else:
            self.custom_species_id = None

    def _get_species(self) -> int | str:
        """Return the current species"""
        if self.custom_species_id:
            return self.custom_species_id
        return 9606  # default: Human

    # ─── Existing API Methods (v12.5 upgraded) ────────────────────

    def resolve_gene_symbol(self, protein_name: str) -> Optional[str]:
        """Convert a protein name to a gene symbol"""
        if not protein_name:
            return None

        name = protein_name.strip().lstrip('-')
        if '@' in name:
            name = name.split('@')[0]

        if not name:
            return None

        # check cache
        gene_cached = self.cache.get_gene(name)
        if gene_cached:
            self.stats["cache_hits"] += 1
            return gene_cached

        self.stats["cache_misses"] += 1

        # STRING API: map identifiers
        url = f"{self.base_url}/json/get_string_ids"
        params = {
            "identifiers": name,
            "species": self._get_species(),
            "limit": 1
        }

        data = self._get_json(url, params)

        if not data or len(data) == 0:
            return None

        gene = data[0].get("preferredName", "")
        if gene:
            self.cache.put_gene(name, gene)
            return gene

        return None

    def get_interacting_proteins(
            self,
            gene_symbol: str,
            min_score: int = 700,
            return_quality_info: bool = False
    ) -> List:
        """
        Query interacting proteins with STRING API v12.5
        - uses functional/physical/regulatory networks according to network_type
        - regulatory mode: includes directionality and regulatory type (activation/inhibition) information

        Returns:
            return_quality_info=False: list of gene symbols
            return_quality_info=True: list of dictionaries
              (gene, score, weight, confidence + regulatory fields)
        """
        if not gene_symbol:
            return []

        # check cache
        cached = self.cache.get_neighbors(gene_symbol, f"{self.network_type}_{self._get_species()}", min_score)
        if cached is not None:
            self.stats["cache_hits"] += 1
            with self.snapshot_lock:
                self.snapshot_neighbors[gene_symbol] = {"origin": "v12.5_cache", "neighbors": [dict(x) for x in cached]}
            filtered = [inter for inter in cached if inter["score"] >= min_score]
            if return_quality_info:
                return filtered
            else:
                return [inter["gene"] for inter in filtered]

        self.stats["cache_misses"] += 1

        # Call STRING API - add network_type parameter (v12.5)
        url = f"{self.base_url}/json/interaction_partners"
        params = {
            "identifiers": gene_symbol,
            "species": self._get_species(),
            # Ask STRING for the actual analysis threshold and an explicitly large
            # partner set.  Omitting limit previously yielded only a small default
            # neighborhood, making H=50/300 and baseline screening ineffective.
            "required_score": int(min_score),
            "limit": 10000,
            "network_type": self._api_network_type  # functional, physical, regulatory
        }

        data = self._get_json(url, params)
        if data is None:
            raise RuntimeError(f"STRING v12.5 API failed for {gene_symbol}; cannot create a trustworthy snapshot")
        if not isinstance(data, list):
            raise ValueError(f"Unexpected STRING response for {gene_symbol}: {type(data).__name__}")

        neighbors = []

        if data:
            for inter in data:
                if self.network_type == "regulatory":
                    # ─── STRING v12.5 official directed regulatory schema ───
                    # Empirically verified (2026-09) against
                    # version-12-5.string-db.org/api/json/interaction_partners:
                    # fields are source_string_id, target_string_id,
                    # source_preferred_name, target_preferred_name,
                    # sign ("pos"/"neg"), combined_score, experimental_score,
                    # database_score, textmining_score. There is NO "score",
                    # "preferredName_A/B", "action", "directionality",
                    # "is_directed", "score_A_to_B" or "score_B_to_A" field —
                    # the previous implementation assumed field names that do
                    # not exist in the actual API response, which silently
                    # dropped every regulatory edge.
                    src_name = inter.get("source_preferred_name", "") or inter.get("source_string_id", "")
                    tgt_name = inter.get("target_preferred_name", "") or inter.get("target_string_id", "")

                    if not src_name or not tgt_name or "combined_score" not in inter or "sign" not in inter:
                        raise ValueError("Official STRING 12.5 regulatory response missing required directed fields; refusing functional fallback")

                    # Forward BFS requires OUTGOING regulatory neighbors only.
                    # Incoming source->query rows must not consume the H cap or be
                    # mistaken for query->source traversal candidates.
                    if src_name.casefold() != gene_symbol.casefold():
                        continue
                    partner = tgt_name

                    score = inter.get("combined_score", 0)
                    score_1000 = int(score * 1000)
                    weight = -math.log(score) if score > 0 else float('inf')

                    entry = {
                        "gene": partner,
                        "score": score_1000,
                        "weight": weight,
                        "confidence": "high" if score_1000 >= 700 else "medium" if score_1000 >= 400 else "low"
                    }

                    entry["source"] = src_name
                    entry["target"] = tgt_name
                    entry["is_directed"] = True  # official regulatory edges are always directed

                    sign = inter.get("sign", "")
                    if sign == "pos":
                        action_type = "activation"
                    elif sign == "neg":
                        action_type = "inhibition"
                    else:
                        action_type = sign or "unknown"
                    entry["action"] = action_type
                    entry["directionality"] = f"{src_name}->{tgt_name}"

                    # official regulatory has no separate A->B / B->A scores;
                    # combined_score already reflects this directed edge.
                    entry["score_A_to_B"] = score
                    entry["score_B_to_A"] = 0

                    # update statistics
                    self.stats["regulatory_edges"] += 1
                    if action_type == "activation":
                        self.stats["activation_edges"] += 1
                    elif action_type == "inhibition":
                        self.stats["inhibition_edges"] += 1
                    self.stats["directional_edges"] += 1

                else:
                    # ─── functional / physical (unchanged) ───
                    partner = inter.get("preferredName_B", "") or inter.get("stringId_B", "")
                    score = inter.get("score", 0)

                    if not partner:
                        continue

                    score_1000 = int(score * 1000)
                    weight = -math.log(score) if score > 0 else float('inf')

                    entry = {
                        "gene": partner,
                        "score": score_1000,
                        "weight": weight,
                        "confidence": "high" if score_1000 >= 700 else "medium" if score_1000 >= 400 else "low"
                    }

                self.stats["edges_retrieved"] += 1
                self.stats["avg_score"] += score_1000

                neighbors.append(entry)

        neighbors.sort(key=lambda x: x["score"], reverse=True)
        self.cache.put_neighbors(gene_symbol, neighbors, f"{self.network_type}_{self._get_species()}", min_score)
        with self.snapshot_lock:
            self.snapshot_neighbors[gene_symbol] = {"origin": "v12.5_API", "neighbors": [dict(x) for x in neighbors]}

        filtered = [inter for inter in neighbors if inter["score"] >= min_score]
        self.log(
            f"[STRING Neighborhood] {gene_symbol}: {len(filtered)} partners at score >= {int(min_score)} "
            f"(requested limit=10000, network={self.network_type})"
        )

        if return_quality_info:
            return filtered
        else:
            return [inter["gene"] for inter in filtered]

    def get_edge_weight(self, gene_a: str, gene_b: str) -> float:
        """Return the edge weight between two genes"""
        neighbors = self.get_interacting_proteins(gene_a, min_score=0, return_quality_info=True)

        for inter in neighbors:
            if inter["gene"] == gene_b:
                return inter["weight"]

        return float('inf')

    def get_regulatory_info(self, gene_a: str, gene_b: str) -> Optional[Dict]:
        """
        Return regulatory relationship information between two genes (regulatory mode only)
        Returns: {directionality, action, is_directed, score_A_to_B, score_B_to_A}
        """
        if self.network_type not in ("regulatory",):
            return None

        neighbors = self.get_interacting_proteins(gene_a, min_score=0, return_quality_info=True)
        for inter in neighbors:
            if (inter["gene"] == gene_b and inter.get("is_directed")
                    and inter.get("source", "").upper() == gene_a.upper()
                    and inter.get("target", "").upper() == gene_b.upper()):
                return {
                    "source": inter.get("source", gene_a),
                    "target": inter.get("target", gene_b),
                    "directionality": inter.get("directionality", "unknown"),
                    "action": inter.get("action", ""),
                    "is_directed": inter.get("is_directed", False),
                    "score_A_to_B": inter.get("score_A_to_B", 0),
                    "score_B_to_A": inter.get("score_B_to_A", 0),
                    "score": inter.get("score", 0)
                }
        return None

    def print_stats_brief(self):
        total_requests = self.stats["cache_hits"] + self.stats["cache_misses"]
        hit_rate = 0 if total_requests == 0 else self.stats["cache_hits"] / total_requests * 100

        msg = (
            f"[API Status] Calls: {self.stats['api_calls']}, "
            f"Cache hit rate: {hit_rate:.1f}%, "
            f"Edges: {self.stats['edges_retrieved']}"
        )
        if self.network_type in ("regulatory",):
            msg += (f", Regulatory: {self.stats['regulatory_edges']}"
                    f" (activation: {self.stats['activation_edges']}, "
                    f"inhibition: {self.stats['inhibition_edges']}, "
                    f"Directional: {self.stats['directional_edges']})")
        self.log(msg)

    def print_stats(self):
        total_requests = self.stats["cache_hits"] + self.stats["cache_misses"]
        hit_rate = 0 if total_requests == 0 else self.stats["cache_hits"] / total_requests * 100

        avg_score = 0
        if self.stats["edges_retrieved"] > 0:
            avg_score = self.stats["avg_score"] / self.stats["edges_retrieved"]

        msg = (
            f"\n{'=' * 60}\n"
            f"STRING API v12.5 Statistics (Network Type: {self.network_type}):\n"
            f"  API Base URL: {self.base_url}\n"
            f"  Total API calls: {self.stats['api_calls']}\n"
            f"  Cache hits: {self.stats['cache_hits']}\n"
            f"  Cache misses: {self.stats['cache_misses']}\n"
            f"  Cache hit rate: {hit_rate:.1f}%\n"
            f"  Total API time: {self.stats['total_time']:.2f} sec\n"
            f"  Retrieved edges: {self.stats['edges_retrieved']}\n"
            f"  Average Combined Score: {avg_score:.1f}/1000\n"
        )

        if self.network_type in ("regulatory",):
            msg += (
                f"\n  === Regulatory Network Statistics ===\n"
                f"  Regulatory edges: {self.stats['regulatory_edges']}\n"
                f"  Activation edges: {self.stats['activation_edges']}\n"
                f"  Inhibition edges: {self.stats['inhibition_edges']}\n"
                f"  Directional edges: {self.stats['directional_edges']}\n"
            )

        if self.sequence_mode:
            msg += f"\n  Sequence input mode: enabled\n"
        if self.custom_species_id:
            msg += f"  Custom Species ID: {self.custom_species_id}\n"

        msg += f"{'=' * 60}"
        self.log(msg)




# ========== Cell-Context Bridge Filter ==========


def _norm_cell_name(x: str) -> str:
    """Normalize cell-line names only for exact alias matching (no fuzzy guessing)."""
    return re.sub(r"[^A-Z0-9]+", "", str(x or "").upper())


def _norm_gene_symbol(x: str) -> str:
    s = str(x or "").strip()
    if not s:
        return ""
    if '@' in s:
        s = s.split('@', 1)[0]
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s).strip()  # e.g. TP53 (7157)
    return s.upper()


def _open_delimited_text(path: str):
    """Open plain/gz/zip delimited files in text mode without adding dependencies."""
    low = path.lower()
    if low.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8-sig', errors='replace', newline='')
    if low.endswith('.zip'):
        zf = zipfile.ZipFile(path)
        members = [n for n in zf.namelist() if not n.endswith('/')]
        if not members:
            zf.close()
            raise ValueError(f"ZIP contains no files: {path}")
        preferred = sorted(members, key=lambda n: (0 if n.lower().endswith(('.tsv','.csv','.txt')) else 1, len(n)))[0]
        raw = zf.open(preferred, 'r')
        wrapper = io.TextIOWrapper(raw, encoding='utf-8-sig', errors='replace', newline='')
        wrapper._cellctx_zipfile = zf
        return wrapper
    return open(path, 'r', encoding='utf-8-sig', errors='replace', newline='')


def _detect_delimiter(path: str, sample: str = "") -> str:
    if path.lower().endswith(('.tsv', '.txt', '.tsv.gz', '.tsv.zip')):
        return '\t'
    if sample:
        try:
            return csv.Sniffer().sniff(sample, delimiters=',\t;').delimiter
        except Exception:
            pass
    return ','


class CellContextFilter:
    """
    Cell-line-aware bridge feasibility filter.

    Core safety rule:
      * Experimental anchors are NEVER removed by cell-context bulk data.
      * Only STRING-only intermediate bridge candidates are filtered.
      * Conservative Union (default) keeps UNKNOWN genes.
      * A bridge is blocked in Conservative Union only if RNA is explicitly
        below threshold AND proteomics is explicitly negative.

    Proteomics numeric matrices are interpreted conservatively:
      any finite numeric value = detected; missing/NA = unknown.
      Zero is treated as negative only for clearly binary 0/1 matrices.
    """

    MODES = ('off', 'conservative_union', 'rna_only', 'proteome_only', 'strict_intersection')

    def __init__(self, config=None, logger_print=None):
        self.config = dict(config or {})
        self._log = logger_print
        self.enabled = bool(self.config.get('enabled', False))
        self.mode = str(self.config.get('mode', 'conservative_union') or 'conservative_union')
        if self.mode not in self.MODES:
            self.mode = 'conservative_union'
        self.cell_line = str(self.config.get('cell_line', '') or '').strip()
        self.rna_threshold = float(self.config.get('rna_threshold', 1.0) or 1.0)
        self.rna_file = str(self.config.get('rna_file', '') or '').strip()
        self.proteome_file = str(self.config.get('proteome_file', '') or '').strip()
        self.metadata_file = str(self.config.get('metadata_file', '') or '').strip()

        self.aliases = set()
        self.identity = {'requested_cell_line': self.cell_line, 'verified': False, 'aliases': []}
        self.experimental_genes = set()

        self.rna_present = set()
        self.rna_negative = set()
        self.rna_values = {}
        self.rna_source_info = {}

        self.protein_present = set()
        self.protein_negative = set()
        self.protein_values = {}
        self.protein_source_info = {}

        self.decisions = {}
        self.stats = Counter()

    def log(self, msg):
        if self._log:
            self._log(msg)

    @staticmethod
    def _header_map(header):
        return {re.sub(r'[^a-z0-9]+', '', str(h).lower()): i for i, h in enumerate(header)}

    @staticmethod
    def _find_col(header, candidates):
        norm = [re.sub(r'[^a-z0-9]+', '', str(h).lower()) for h in header]
        cand = {re.sub(r'[^a-z0-9]+', '', c.lower()) for c in candidates}
        for i, n in enumerate(norm):
            if n in cand:
                return i
        return None

    @staticmethod
    def _find_col_preferred(header, candidates):
        """Return the first available column by candidate priority order."""
        norm_to_index = {
            re.sub(r'[^a-z0-9]+', '', str(h).lower()): i
            for i, h in enumerate(header)
        }
        for c in candidates:
            key = re.sub(r'[^a-z0-9]+', '', str(c).lower())
            if key in norm_to_index:
                return norm_to_index[key]
        return None

    def _resolve_aliases(self):
        if not self.cell_line:
            raise ValueError('Cell-context filtering is enabled, but no cell line was specified.')
        self.aliases = {_norm_cell_name(self.cell_line)}
        cmp_id = str(self.config.get('cmp_model_id') or '').strip()
        if cmp_id:
            self.aliases.add(_norm_cell_name(cmp_id))
        matched_rows = []
        if self.metadata_file:
            if not os.path.isfile(self.metadata_file):
                raise ValueError(f'Cell-line metadata file not found: {self.metadata_file}')
            fh = _open_delimited_text(self.metadata_file)
            try:
                sample = fh.read(8192); fh.seek(0)
                delim = _detect_delimiter(self.metadata_file, sample)
                reader = csv.reader(fh, delimiter=delim)
                header = next(reader, [])
                if not header:
                    raise ValueError('Cell-line metadata file has no header.')
                name_cols=[]
                for i,h in enumerate(header):
                    nh=re.sub(r'[^a-z0-9]+','',h.lower())
                    if nh in {
                        'cellline','celllinename','strippedcelllinename','modelname','modelid',
                        'depmapid','depmapmodelid','ccle_name','cclename','sangerid','modelid'
                    }:
                        name_cols.append(i)
                if not name_cols:
                    raise ValueError('Could not identify cell-line/model columns in metadata file.')
                q=_norm_cell_name(self.cell_line)
                for row in reader:
                    vals=[row[i].strip() for i in name_cols if i < len(row) and row[i].strip()]
                    if any(_norm_cell_name(v)==q for v in vals):
                        matched_rows.append(vals)
                if len(matched_rows) > 1:
                    # Distinct duplicate matches are unsafe; require user correction.
                    uniq={tuple(v) for v in matched_rows}
                    if len(uniq) > 1:
                        raise ValueError(f"Cell-line identity is ambiguous in metadata for '{self.cell_line}' ({len(uniq)} matches).")
                if matched_rows:
                    for v in matched_rows[0]:
                        self.aliases.add(_norm_cell_name(v))
            finally:
                try: fh.close()
                except Exception: pass
        self.aliases.discard('')
        self.identity = {
            'requested_cell_line': self.cell_line,
            'verified': bool(self.aliases),
            'aliases': sorted(self.aliases),
            'metadata_file': self.metadata_file or None,
        }

    def _classify_rna_value(self, raw, value_col_name=''):
        s=str(raw or '').strip()
        if not s or s.upper() in {'NA','NAN','NULL','NONE','-'}:
            return None, None
        try:
            val=float(s)
        except Exception:
            return None, None
        # Auto-handle common log2(TPM+1) matrices by column/file naming.
        name=(value_col_name + ' ' + os.path.basename(self.rna_file)).lower()
        threshold=self.rna_threshold
        if 'logp1' in name or 'log2' in name:
            threshold=math.log2(self.rna_threshold + 1.0)
        return val, (val > threshold)

    @staticmethod
    def _scan_table_header(fh, delimiter, required_groups, max_lines=40):
        """Find a real header row even when downloaded tables start with metadata lines."""
        reader = csv.reader(fh, delimiter=delimiter)
        groups = [
            {re.sub(r'[^a-z0-9]+', '', str(x).lower()) for x in grp}
            for grp in required_groups
        ]
        rows = []
        for idx, row in enumerate(reader):
            rows.append(row)
            norm = {re.sub(r'[^a-z0-9]+', '', str(x).lower()) for x in row}
            if all(norm & grp for grp in groups):
                return row, idx
            if idx + 1 >= max_lines:
                break
        return (rows[0] if rows else []), 0

    def _load_rna(self):
        if not self.rna_file:
            return
        if not os.path.isfile(self.rna_file):
            raise ValueError(f'RNA bulk file not found: {self.rna_file}')
        fh=_open_delimited_text(self.rna_file)
        try:
            sample=fh.read(16384); fh.seek(0)
            delim=_detect_delimiter(self.rna_file, sample)

            header, skipped = self._scan_table_header(
                fh, delim,
                [
                    ['Gene','Gene name','gene_name','gene_symbol','Symbol'],
                    ['rsem_tpm','RSEM TPM','nTPM','TPM','TPM value','TPM_value','pTPM',
                     'expression','rna_expression','log2(TPM+1)','TPM_log2','rsem_fpkm','htseq_fpkm']
                ],
                max_lines=40
            )
            if not header:
                raise ValueError('RNA file has no readable header.')

            try: fh.close()
            except Exception: pass
            fh=_open_delimited_text(self.rna_file)
            reader=csv.reader(fh, delimiter=delim)
            for _ in range(skipped + 1):
                next(reader, None)

            cell_i=self._find_col_preferred(header, ['Cell line','Cell_line','cell_line','cellLine',
                                                    'model_name','modelName','CellLineName','ModelID','model_id','Sanger ID'])
            gene_i=self._find_col_preferred(header, ['gene_symbol','Symbol','Gene','Gene name','gene_name'])
            val_i=self._find_col_preferred(header, ['rsem_tpm','RSEM TPM','nTPM','TPM','TPM value','TPM_value',
                                                   'pTPM','expression','rna_expression','log2(TPM+1)','TPM_log2',
                                                   'rsem_fpkm','htseq_fpkm'])
            self.log(
                f"[Cell Context] RNA parser columns: "
                f"gene={header[gene_i] if gene_i is not None else None}, "
                f"value={header[val_i] if val_i is not None else None}, "
                f"cell/model={header[cell_i] if cell_i is not None else None}"
            )
            matched=0
            cmp_model_specific = bool(
                self.config.get('cmp_model_id')
                or str(self.config.get('auto_source') or '').lower().startswith('cell model passports')
                or re.search(r'SIDM\d+.*rnaseq', os.path.basename(self.rna_file), re.I)
            )
            if gene_i is not None and val_i is not None:
                # Long table. Official Cell Model Passports downloads are already
                # scoped to one SIDM. In that case do NOT re-filter rows by model_name;
                # the selected SIDM itself is the identity guarantee.
                for row in reader:
                    if max(gene_i,val_i) >= len(row): continue
                    if (not cmp_model_specific) and cell_i is not None:
                        if cell_i >= len(row) or _norm_cell_name(row[cell_i]) not in self.aliases:
                            continue
                    g=_norm_gene_symbol(row[gene_i])
                    if not g: continue
                    val,pos=self._classify_rna_value(row[val_i], header[val_i])
                    if val is None: continue
                    matched += 1
                    self.rna_values[g]=val
                    (self.rna_present if pos else self.rna_negative).add(g)
                self.rna_source_info={
                    'layout':'cmp_model_specific_long' if cmp_model_specific else ('single_model_long' if cell_i is None else 'long'),
                    'cell_column':header[cell_i] if cell_i is not None else None,
                    'gene_column':header[gene_i],
                    'value_column':header[val_i],
                    'matched_rows':matched
                }
            else:
                # Matrix layout. First try cell line as a column (genes in rows).
                header_alias_idx=None
                for i,h in enumerate(header[1:],1):
                    if _norm_cell_name(h) in self.aliases:
                        header_alias_idx=i; break
                if header_alias_idx is not None:
                    for row in reader:
                        if len(row) <= header_alias_idx: continue
                        g=_norm_gene_symbol(row[0])
                        if not g: continue
                        val,pos=self._classify_rna_value(row[header_alias_idx], header[header_alias_idx])
                        if val is None: continue
                        matched += 1; self.rna_values[g]=val
                        (self.rna_present if pos else self.rna_negative).add(g)
                    self.rna_source_info={'layout':'matrix_gene_rows','selected_column':header[header_alias_idx],'matched_genes':matched}
                else:
                    # Cell line as a row (DepMap-style; genes in columns).
                    found_row=None
                    for row in reader:
                        if row and _norm_cell_name(row[0]) in self.aliases:
                            found_row=row; break
                    if found_row is None:
                        is_cmp_specific = bool(
                            self.config.get('cmp_model_id')
                            or str(self.config.get('auto_source') or '').lower().startswith('cell model passports')
                            or re.search(r'SIDM\d+.*rnaseq', os.path.basename(self.rna_file), re.I)
                        )
                        if is_cmp_specific:
                            raise ValueError(
                                "Official Cell Model Passports RNA file is model-specific, "
                                "but the gene/expression columns could not be recognized. "
                                "Expected CMP fields include gene_symbol + rsem_tpm. "
                                f"Detected header: {header[:12]}"
                            )
                        raise ValueError(f"Exact cell-line/model match not found in RNA file for '{self.cell_line}'. Provide the correct metadata file or cell-line ID.")
                    for i in range(1, min(len(header),len(found_row))):
                        g=_norm_gene_symbol(header[i])
                        if not g: continue
                        val,pos=self._classify_rna_value(found_row[i], header[i])
                        if val is None: continue
                        matched += 1; self.rna_values[g]=val
                        (self.rna_present if pos else self.rna_negative).add(g)
                    self.rna_source_info={'layout':'matrix_cell_rows','selected_row':found_row[0],'matched_genes':matched}
            if matched == 0:
                raise ValueError(
                    f"RNA file produced zero usable gene measurements for '{self.cell_line}'. "
                    f"Detected columns: gene={header[gene_i] if gene_i is not None else None}, "
                    f"value={header[val_i] if val_i is not None else None}, "
                    f"cell/model={header[cell_i] if cell_i is not None else None}, "
                    f"CMP_model_specific={cmp_model_specific}."
                )
        finally:
            try: fh.close()
            except Exception: pass

    @staticmethod
    def _prot_token(raw):
        s=str(raw or '').strip()
        if not s or s.upper() in {'NA','NAN','NULL','NONE','-'}:
            return None
        sl=s.lower()
        if sl in {'true','yes','present','detected','expressed','positive'}: return True
        if sl in {'false','no','absent','not detected','not_detected','undetected','negative'}: return False
        try:
            return float(s)
        except Exception:
            return None

    def _load_proteome(self):
        if not self.proteome_file:
            return
        if not os.path.isfile(self.proteome_file):
            raise ValueError(f'Proteomics bulk file not found: {self.proteome_file}')
        fh=_open_delimited_text(self.proteome_file)
        try:
            sample=fh.read(16384); fh.seek(0)
            delim=_detect_delimiter(self.proteome_file, sample)

            header, skipped = self._scan_table_header(
                fh, delim,
                [
                    ['Gene','Gene name','gene_name','gene_symbol','Symbol'],
                    ['detected','present','protein_detected','protein_expression',
                     'Protein Intensity','protein_intensity','intensity','abundance','protein']
                ],
                max_lines=40
            )
            if not header:
                raise ValueError('Proteomics file has no readable header.')

            try: fh.close()
            except Exception: pass
            fh=_open_delimited_text(self.proteome_file)
            reader=csv.reader(fh, delimiter=delim)
            for _ in range(skipped + 1):
                next(reader, None)

            cell_i=self._find_col(header, ['Cell line','Cell_line','cell_line','cellLine',
                                          'model_name','modelName','CellLineName','ModelID','model_id','Sanger ID'])
            gene_i=self._find_col(header, ['Gene','Gene name','gene_name','gene_symbol','Symbol'])
            val_i=self._find_col(header, ['detected','present','protein_detected','protein_expression',
                                         'Protein Intensity','protein_intensity','intensity','abundance','protein'])
            matched=0
            cmp_model_specific = bool(
                self.config.get('cmp_model_id')
                or str(self.config.get('auto_source') or '').lower().startswith('cell model passports')
                or re.search(r'SIDM\d+.*proteomics', os.path.basename(self.proteome_file), re.I)
            )

            def apply(g, token, binary_hint=False):
                nonlocal matched
                if token is None: return
                matched += 1
                if isinstance(token, bool):
                    if token: self.protein_present.add(g)
                    else: self.protein_negative.add(g)
                    self.protein_values[g]=token
                elif isinstance(token, (int,float)) and math.isfinite(token):
                    self.protein_values[g]=float(token)
                    if binary_hint and token in (0,0.0): self.protein_negative.add(g)
                    else: self.protein_present.add(g)

            if gene_i is not None and val_i is not None:
                # Long table; model-specific Cell Model Passports proteomics files
                # may have no cell-line column because the whole file is one SIDM.
                for row in reader:
                    if max(gene_i,val_i) >= len(row): continue
                    if (not cmp_model_specific) and cell_i is not None:
                        if cell_i >= len(row) or _norm_cell_name(row[cell_i]) not in self.aliases:
                            continue
                    g=_norm_gene_symbol(row[gene_i]); tok=self._prot_token(row[val_i])
                    if g:
                        apply(g,tok, binary_hint=('detect' in header[val_i].lower() or 'present' in header[val_i].lower()))
                self.protein_source_info={
                    'layout':'cmp_model_specific_long' if cmp_model_specific else ('single_model_long' if cell_i is None else 'long'),
                    'cell_column':header[cell_i] if cell_i is not None else None,
                    'gene_column':header[gene_i],
                    'value_column':header[val_i],
                    'matched_rows':matched
                }
            else:
                header_alias_idx=None
                for i,h in enumerate(header[1:],1):
                    if _norm_cell_name(h) in self.aliases:
                        header_alias_idx=i; break
                if header_alias_idx is not None:
                    vals=[]
                    rows=[]
                    for row in reader:
                        if len(row) <= header_alias_idx: continue
                        g=_norm_gene_symbol(row[0]); tok=self._prot_token(row[header_alias_idx])
                        if g and tok is not None:
                            vals.append(tok); rows.append((g,tok))
                    nums=[v for v in vals if isinstance(v,(int,float)) and not isinstance(v,bool)]
                    binary_hint=bool(nums) and set(nums).issubset({0.0,1.0})
                    for g,tok in rows: apply(g,tok,binary_hint)
                    self.protein_source_info={'layout':'matrix_gene_rows','selected_column':header[header_alias_idx],'matched_genes':matched,'binary_matrix':binary_hint}
                else:
                    found_row=None
                    for row in reader:
                        if row and _norm_cell_name(row[0]) in self.aliases:
                            found_row=row; break
                    if found_row is None:
                        is_cmp_specific = bool(
                            self.config.get('cmp_model_id')
                            or str(self.config.get('auto_source') or '').lower().startswith('cell model passports')
                            or re.search(r'SIDM\d+.*proteomics', os.path.basename(self.proteome_file), re.I)
                        )
                        if is_cmp_specific:
                            raise ValueError(
                                "Official Cell Model Passports proteomics file is model-specific, "
                                "but the Symbol/Protein Intensity columns could not be recognized. "
                                f"Detected header: {header[:12]}"
                            )
                        raise ValueError(f"Exact cell-line/model match not found in proteomics file for '{self.cell_line}'. Provide the correct metadata file or cell-line ID.")
                    raw=[]
                    for i in range(1,min(len(header),len(found_row))):
                        g=_norm_gene_symbol(header[i]); tok=self._prot_token(found_row[i])
                        if g and tok is not None: raw.append((g,tok))
                    nums=[v for _,v in raw if isinstance(v,(int,float)) and not isinstance(v,bool)]
                    binary_hint=bool(nums) and set(nums).issubset({0.0,1.0})
                    for g,tok in raw: apply(g,tok,binary_hint)
                    self.protein_source_info={'layout':'matrix_cell_rows','selected_row':found_row[0],'matched_genes':matched,'binary_matrix':binary_hint}
            if matched == 0:
                raise ValueError(f"Proteomics file matched the cell line but produced zero usable protein measurements for '{self.cell_line}'.")
        finally:
            try: fh.close()
            except Exception: pass

    def build(self):
        if not self.enabled or self.mode == 'off':
            self.enabled=False
            self.log('[Cell Context] Disabled.')
            return self
        if not self.rna_file and not self.proteome_file:
            raise ValueError('Cell-context filtering is enabled, but neither RNA nor proteomics bulk file was selected.')
        self._resolve_aliases()
        self.log(f"[Cell Context] Cell line requested: {self.cell_line}")
        if self.metadata_file:
            self.log(f"[Cell Context] Identity aliases resolved from metadata: {', '.join(sorted(self.aliases))}")
        if self.config.get('cmp_model_id'):
            self.log(
                f"[Cell Context] Official model-specific data: {self.config.get('cmp_model_id')}. "
                "The RNA/proteomics file is already scoped to one SIDM; a cell-line-name column "
                "inside the table is therefore not required."
            )
        else:
            self.log('[Cell Context] No metadata file supplied; exact cell-line/model name match is required in each bulk file.')
        if self.rna_file:
            self._load_rna()
            self.log(f"[Cell Context] RNA loaded: expressed={len(self.rna_present)}, below-threshold={len(self.rna_negative)}; source={self.rna_source_info}")
        if self.proteome_file:
            self._load_proteome()
            self.log(f"[Cell Context] Proteomics loaded: detected={len(self.protein_present)}, explicit-negative={len(self.protein_negative)}; source={self.protein_source_info}")
        self.identity['verified']=True
        self.log(f"[Cell Context] Mode={self.mode}; RNA threshold={self.rna_threshold:g}; UNKNOWN bridge candidates are retained in conservative mode.")
        return self

    def set_experimental_genes(self, genes):
        self.experimental_genes={_norm_gene_symbol(g) for g in (genes or []) if _norm_gene_symbol(g)}
        if self.enabled:
            self.log(f"[Cell Context] Experimental-anchor protection list: {len(self.experimental_genes)} genes.")

    def bridge_allowed(self, gene):
        g=_norm_gene_symbol(gene)
        if not self.enabled or self.mode == 'off' or not g:
            return True
        if g in self.experimental_genes:
            reason='experimental_gene_protected'; allow=True
        else:
            rp=g in self.rna_present; rn=g in self.rna_negative
            pp=g in self.protein_present; pn=g in self.protein_negative
            if self.mode == 'conservative_union':
                if rp or pp:
                    allow=True; reason='positive_RNA_or_protein'
                elif rn and pn:
                    allow=False; reason='RNA_below_threshold_AND_explicit_protein_negative'
                else:
                    allow=True; reason='unknown_retained'
            elif self.mode == 'rna_only':
                if rp: allow=True; reason='RNA_positive'
                elif rn: allow=False; reason='RNA_below_threshold'
                else: allow=True; reason='RNA_unknown_retained'
            elif self.mode == 'proteome_only':
                if pp: allow=True; reason='protein_detected'
                elif pn: allow=False; reason='protein_explicit_negative'
                else: allow=True; reason='protein_unknown_retained'
            else:  # strict_intersection, sensitivity analysis only
                allow=bool(rp and pp)
                reason='RNA_AND_protein_positive' if allow else 'strict_requires_RNA_AND_protein_positive'
        self.stats['bridge_candidates_checked'] += 1
        self.stats['allowed' if allow else 'blocked'] += 1
        self.stats[reason] += 1
        if g not in self.decisions:
            self.decisions[g]={
                'gene':g,
                'allowed':allow,
                'reason':reason,
                'rna_status':'present' if g in self.rna_present else ('below_threshold' if g in self.rna_negative else 'unknown'),
                'rna_value':self.rna_values.get(g),
                'protein_status':'detected' if g in self.protein_present else ('explicit_negative' if g in self.protein_negative else 'unknown'),
                'protein_value':self.protein_values.get(g),
                'experimental_protected':g in self.experimental_genes,
            }
        return allow

    def summary(self):
        return {
            'enabled': self.enabled,
            'mode': self.mode,
            'cell_line': self.cell_line,
            'cmp_model_id': self.config.get('cmp_model_id') or None,
            'auto_source': self.config.get('auto_source') or None,
            'identity': self.identity,
            'rna_threshold': self.rna_threshold,
            'rna_file': self.rna_file or None,
            'proteome_file': self.proteome_file or None,
            'metadata_file': self.metadata_file or None,
            'rna_source_info': self.rna_source_info,
            'protein_source_info': self.protein_source_info,
            'rna_present_count': len(self.rna_present),
            'rna_below_threshold_count': len(self.rna_negative),
            'protein_detected_count': len(self.protein_present),
            'protein_explicit_negative_count': len(self.protein_negative),
            'experimental_protected_count': len(self.experimental_genes),
            'decision_stats': dict(self.stats),
            'policy': (
                'bridge-only; experimental anchors always protected; conservative_union keeps UNKNOWN; '
                'numeric proteomics values are positive evidence, missing values are UNKNOWN unless the source is explicitly binary/categorical'
            ),
        }

    def write_audit_csv(self, path):
        fields=['gene','allowed','reason','rna_status','rna_value','protein_status','protein_value','experimental_protected']
        with open(path,'w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
            for g in sorted(self.decisions): w.writerow(self.decisions[g])


# ========== Temporal Input Adapter (Stage 3) ==========

class TemporalInputAdapter:
    """Normalize legacy v1 and Temporal v2 JSON inputs without altering source data."""

    def __init__(self, data: dict):
        self.data = data or {}
        self.is_temporal_v2 = bool(self.data.get("intervals")) and str(self.data.get("schema_version", "")).startswith("2")
        self.schema_label = "Temporal v2 interval model" if self.is_temporal_v2 else "Legacy v1 time-point model"

    @staticmethod
    def tokens_to_gene_map(tokens, gene_resolver):
        out = {}
        for tok in tokens or []:
            raw = _as_name(tok).strip()
            if not raw:
                continue
            gene = gene_resolver(raw)
            if gene:
                out.setdefault(gene, []).append(raw)
        return out

    def all_site_tokens(self):
        if not self.is_temporal_v2:
            for arr in (self.data.get("time_series") or {}).values():
                for tok in arr or []:
                    yield tok
            return
        baseline = self.data.get("baseline") or {}
        for tok in baseline.get("available", []) or []:
            yield tok
        bc = self.data.get("baseline_contrast") or {}
        for key in ("increased", "decreased", "response_anchors", "response_supporting", "response_opposing", "response_unresolved"):
            for tok in bc.get(key, []) or []:
                yield tok
        for interval in self.data.get("intervals", []) or []:
            for key in ("increased", "decreased", "response_anchors", "response_supporting", "response_opposing", "response_unresolved"):
                for tok in interval.get(key, []) or []:
                    yield tok

    def condition_functional_pool(self, gene_resolver):
        """Build a gene-level condition-specific pool from pretreatment-vs-reference Q states.

        Only genes with opposing evidence and NO supporting evidence are excluded.
        Mixed supporting/opposing genes remain eligible.  No drug-name inference is used.
        """
        bc = self.data.get("baseline_contrast") or {}
        supporting_map = self.tokens_to_gene_map(bc.get("response_supporting", []) or [], gene_resolver)
        opposing_map = self.tokens_to_gene_map(bc.get("response_opposing", []) or [], gene_resolver)
        unresolved_map = self.tokens_to_gene_map(bc.get("response_unresolved", []) or [], gene_resolver)
        s = set(supporting_map)
        o = set(opposing_map)
        u = set(unresolved_map)
        return {
            "excluded_low_activity": o - s,
            "supporting_only": s - o,
            "mixed": s & o,
            "unresolved_only": u - s - o,
            "supporting_map": supporting_map,
            "opposing_map": opposing_map,
            "unresolved_map": unresolved_map,
            "has_baseline_contrast": bool(bc),
            "baseline_contrast_type": bc.get("type", "") if isinstance(bc, dict) else ""
        }

    def intervals(self):
        return self.data.get("intervals", []) or [] if self.is_temporal_v2 else []


# ========== PhosphoSite Manager ==========

class PhosphoSiteManager:
    """Manage phosphorylation site information"""

    def __init__(self):
        self.protein_sites: Dict[str, Set[str]] = {}
        self.protein_site_combos: Dict[str, List[str]] = {}
        self.stats = {"total_sites": 0, "proteins_with_sites": 0}

    def add_site(self, protein: str, site: str):
        if protein not in self.protein_sites:
            self.protein_sites[protein] = set()
        self.protein_sites[protein].add(site)

    def add_site_combo(self, protein: str, combo_str: str):
        if not combo_str:
            return
        if protein not in self.protein_site_combos:
            self.protein_site_combos[protein] = []
        if combo_str not in self.protein_site_combos[protein]:
            self.protein_site_combos[protein].append(combo_str)

    def get_sites(self, protein: str) -> Set[str]:
        return self.protein_sites.get(protein, set())

    def get_site_combos(self, protein: str) -> List[str]:
        return self.protein_site_combos.get(protein, [])

    def has_sites(self, protein: str) -> bool:
        return protein in self.protein_sites

    def load_from_json(self, json_data: dict, logger_print=None):
        if logger_print:
            logger_print("[PhosphoSite] Loading site information from JSON...")

        adapter = TemporalInputAdapter(json_data)
        for protein_str in adapter.all_site_tokens():
            protein_name, sites = parse_protein_with_site(_as_name(protein_str))
            if sites:
                for site in sites:
                    self.add_site(protein_name, site)
                clean_str = _as_name(protein_str).lstrip('-')
                if '@' in clean_str:
                    combo_str = clean_str.split('@', 1)[1]
                    self.add_site_combo(protein_name, combo_str)

        self.stats["proteins_with_sites"] = len(self.protein_sites)
        self.stats["total_sites"] = sum(len(s) for s in self.protein_sites.values())

        if logger_print:
            logger_print(f"[PhosphoSite] Input schema: {adapter.schema_label}")
            logger_print(f"[PhosphoSite] {self.stats['proteins_with_sites']} proteins, "
                         f"{self.stats['total_sites']} sites loaded")
            total_combos = sum(len(v) for v in self.protein_site_combos.values())
            logger_print(f"[PhosphoSite] {total_combos} site combinations tracked")


# ========== PhosphoPath Generator ==========

def generate_phospho_paths(
        base_paths: List[Tuple],
        site_manager,
        max_combos_per_path: int = 500
) -> List[Tuple[str, ...]]:
    if not site_manager:
        return []

    result_set: Set[Tuple[str, ...]] = set()

    for path in base_paths:
        if not path:
            continue

        gene_options: List[List[str]] = []
        for gene in path:
            combos = site_manager.get_site_combos(gene)
            if combos:
                gene_options.append([f"{gene}@{c}" for c in combos])
            else:
                gene_options.append([gene])

        total_combos = 1
        for opts in gene_options:
            total_combos *= len(opts)

        if total_combos > max_combos_per_path:
            gene_options = [opts[:1] if len(opts) > 3 else opts for opts in gene_options]

        for combo_tuple in itertools.islice(itertools.product(*gene_options), max_combos_per_path):
            result_set.add(combo_tuple)

    return list(result_set)


def generate_temporal_endpoint_phospho_paths(temporal_records: List[Dict]) -> List[Tuple[str, ...]]:
    """
    FINAL phosphosite output for Temporal v2.

    Do NOT form a Cartesian product of every phosphosite on every gene.  Only
    measured temporal target anchors receive a representative @site annotation;
    prior-network bridge nodes remain gene-level.  All target_anchor_tokens are
    still preserved losslessly in the structured temporal JSON.
    """
    out: Set[Tuple[str, ...]] = set()
    for rec in temporal_records or []:
        genes = list(rec.get("genes", []) or [])
        if not genes:
            continue
        labels = list(genes)
        # Map each segment's target occurrence in the ordered gene chain.
        search_from = 0
        for seg in rec.get("segments", []) or []:
            target = str(seg.get("target_anchor", ""))
            toks = sorted(set(str(x) for x in (seg.get("target_anchor_tokens", []) or []) if x))
            if not target or not toks:
                continue
            try:
                idx = genes.index(target, search_from)
            except ValueError:
                try:
                    idx = genes.index(target)
                except ValueError:
                    continue
            # Prefer a token whose gene prefix matches the target.  One representative
            # token prevents combinatorial expansion; the full token list remains JSON.
            rep = next((t for t in toks if t.split('@',1)[0].lstrip('-') == target), toks[0])
            labels[idx] = rep
            search_from = idx
        out.add(tuple(labels))
    return list(out)


def write_phospho_results(
        out_path: str,
        complete_phospho: List[Tuple[str, ...]],
        intermediate_phospho: List[Tuple[str, ...]],
        elapsed: float
):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("*** Path results including phosphorylation sites (all_Phospho_results) ***\n")
        f.write("*** (Base paths are gene-symbol based; even for the same path, different @XXX combinations are written on separate lines) ***\n")
        f.write(f"(Total elapsed time: {elapsed:.1f} sec)\n\n")

        f.write(f"=== Complete paths with phosphorylation sites ({len(complete_phospho)} variants) ===\n\n")
        if complete_phospho:
            for i, path_tuple in enumerate(complete_phospho, 1):
                f.write(f"#{i}: {' -> '.join(path_tuple)}\n")
        else:
            f.write("No phosphorylation-site information or no complete paths are available.\n")

        f.write(f"\n\n=== Intermediate paths with phosphorylation sites ({len(intermediate_phospho)} variants) ===\n\n")
        if intermediate_phospho:
            by_length: Dict[int, List] = {}
            for pt in intermediate_phospho:
                by_length.setdefault(len(pt), []).append(pt)
            for length in sorted(by_length.keys(), reverse=True):
                grp = by_length[length]
                f.write(f"--- Length {length} ({len(grp)} variants) ---\n")
                for i, path_tuple in enumerate(grp[:200], 1):
                    f.write(f"  {i}. {' -> '.join(path_tuple)}\n")
                if len(grp) > 200:
                    f.write(f"  ... additional {len(grp) - 200} more\n")
                f.write("\n")
        else:
            f.write("No phosphorylation-site information or no intermediate paths are available.\n")

        f.write("\n*** Phosphorylation Path Statistics ***\n")
        f.write(f"Complete path variants: {len(complete_phospho)}\n")
        f.write(f"Intermediate path variants: {len(intermediate_phospho)}\n")


# ========== PathwaySimulator (STRING v12.5 - Regulatory Network Support) ==========

class PathwaySimulator:
    def __init__(self, api_client, max_depth=2, advanced=None,
                 logger_print=None, site_manager=None,
                 progress_callback=None, cell_context=None):
        self.api = api_client
        self.max_depth = max_depth
        self._log = logger_print
        self._progress_callback = progress_callback
        self.interrupt = False
        self.name_cache: Dict[str, str] = {}
        self.site_manager = site_manager
        self.temporal_path_records: List[Dict] = []
        self.temporal_interval_summary: List[Dict] = []
        self.temporal_pruning_summary: List[Dict] = []
        self.input_schema_label = "Legacy v1 time-point model"
        self.cell_context = cell_context

        # Collect regulatory information (regulatory mode)
        self.regulatory_info: Dict[Tuple[str, str], Dict] = {}

        # Performance caches. These only reuse deterministic computations; they do
        # not change STRING threshold, H, D, Q/cell-context filtering, direction,
        # beam width, route ranking, or target-preservation semantics.
        self._ranked_neighbor_cache: Dict[Tuple, Tuple[Dict, ...]] = {}
        self._neighbor_lookup_cache: Dict[Tuple, Dict[str, Dict]] = {}
        self._baseline_first_hop_cache: Dict[Tuple, Tuple[Tuple[str, Dict, float], ...]] = {}
        self._segment_memo = OrderedDict()
        self._segment_memo_max = 50000
        self._perf_stats = {
            "ranked_neighbor_hits": 0, "ranked_neighbor_misses": 0,
            "neighbor_lookup_hits": 0, "neighbor_lookup_misses": 0,
            "segment_memo_hits": 0, "segment_memo_misses": 0,
        }

        self.adv = {
            "degree_cap": 50,
            "beam_width": 150,
            "max_paths_per_target": 300,
            "hub_blacklist": [],
            "collapse_gene_level": True,
            "no_immediate_backtrack": True,
            "dedupe_gene_in_path": True,
            "verbose": False,
            "min_string_score": 700,
            "regulatory_weight_bonus": 0.1, # weight bonus for direction-confirmed regulatory edges
            "temporal_anchor_mode": "measured",  # measured | strict_signed
            "exclude_opposing_as_bridge": True,

            # Stage 3.3 temporal safety controls
            "temporal_segment_route_cap": 5,
            "temporal_beam_width": 0,       # 0 = reuse ordinary Beam Width K
            "temporal_routes_per_target": 3,
            "temporal_full_path_dedupe": True
        }
        if isinstance(advanced, dict):
            self.adv.update({k: v for k, v in advanced.items() if k in self.adv})

    def log(self, msg: str):
        if self._log:
            self._log(msg)

    def set_stop_flag(self, flag: bool):
        self.interrupt = flag

    def report_progress(self, fraction: float, phase: str = "",
                        completed: int = 0, total: int = 0):
        # Progress reporting must never alter/interrupt the scientific run.
        if not self._progress_callback:
            return
        try:
            fraction = max(0.0, min(1.0, float(fraction)))
            self._progress_callback(
                fraction,
                str(phase or ""),
                int(completed or 0),
                int(total or 0)
            )
        except Exception:
            pass

    def _format_pathway(self, pathway: List[Tuple[str, float]], show_sites: bool = True) -> str:
        if not pathway or len(pathway) < 2:
            return ""

        result = []
        total_cost = 0.0

        for i in range(len(pathway)):
            gene, cost = pathway[i]

            if show_sites and self.site_manager:
                sites = self.site_manager.get_sites(gene)
                if sites:
                    site_str = '_'.join(sorted(list(sites))[:3])
                    result.append(f"{gene}@{site_str}")
                else:
                    result.append(gene)
            else:
                result.append(gene)

            if i < len(pathway) - 1:
                # In regulatory mode, show the regulatory type in the arrow
                if self.api.network_type in ("regulatory",):
                    next_gene = pathway[i + 1][0]
                    reg = self.regulatory_info.get((gene, next_gene))
                    if reg:
                        action = reg.get("action", "")
                        if action in ("activation", "catalysis", "expression"):
                            result.append(" ──▶ ")  # activation
                        elif action in ("inhibition",):
                            result.append(" ──⊣ ")  # inhibition
                        elif reg.get("is_directed"):
                            result.append(" ──→ ")  # direction confirmed only
                        else:
                            result.append(" ── ")   # direction uncertain
                    else:
                        result.append(" --> ")
                else:
                    result.append(" --> ")
                total_cost += pathway[i + 1][1]

        path_str = ''.join(result)
        return f"{path_str} [cost: {total_cost:.3f}]"

    def _format_pathway_arrow(self, genes: List[str]) -> str:
        """Format a gene list using arrows with regulatory information"""
        if len(genes) < 2:
            return ' -> '.join(genes)

        parts = [genes[0]]
        for i in range(1, len(genes)):
            prev, curr = genes[i-1], genes[i]
            reg = self.regulatory_info.get((prev, curr))
            if reg and self.api.network_type in ("regulatory",):
                action = reg.get("action", "")
                if action in ("activation", "catalysis", "expression"):
                    parts.append(f" ──▶ {curr}")
                elif action in ("inhibition",):
                    parts.append(f" ──⊣ {curr}")
                elif reg.get("is_directed"):
                    parts.append(f" ──→ {curr}")
                else:
                    parts.append(f" -> {curr}")
            else:
                parts.append(f" -> {curr}")
        return ''.join(parts)

    def _gene(self, name: str, include_site: bool = False) -> str:
        if not name:
            return ""

        name_clean, sites = parse_protein_with_site(name)

        if name_clean in self.name_cache:
            gene = self.name_cache[name_clean]
        else:
            g = self.api.resolve_gene_symbol(name_clean) or name_clean
            self.name_cache[name_clean] = g
            gene = g

        if include_site and sites:
            site_str = '_'.join(sorted(list(sites))[:3])
            return f"{gene}@{site_str}"

        return gene

    def _eq_relaxed(self, a: str, b: str) -> bool:
        if not a or not b:
            return False
        if a == b:
            return True
        return self._gene(a) == self._gene(b)

    def _neighbors(self, gene: str) -> List[Dict]:
        return self.api.get_interacting_proteins(
            gene,
            min_score=self.adv.get("min_string_score", 700),
            return_quality_info=True
        ) or []

    def _ranked_neighbors(self, gene: str) -> Tuple[Dict, ...]:
        """Return the exact H-truncated neighbor list used by ordinary BFS."""
        H = int(self.adv.get("degree_cap", 50)) or 50
        hubs = tuple(sorted(set(self.adv.get("hub_blacklist", []))))
        min_score = int(self.adv.get("min_string_score", 700)) or 700
        key = (gene, self.api.network_type, H, min_score, hubs)
        cached = self._ranked_neighbor_cache.get(key)
        if cached is not None:
            self._perf_stats["ranked_neighbor_hits"] += 1
            return cached
        self._perf_stats["ranked_neighbor_misses"] += 1
        hubset=set(hubs)
        neigh=[x for x in list(self._neighbors(gene)) if x.get("gene") and x.get("gene") not in hubset]
        if self.api.network_type == "regulatory":
            directed=sorted([x for x in neigh if x.get("is_directed")], key=lambda x:x["weight"])
            undirected=sorted([x for x in neigh if not x.get("is_directed")], key=lambda x:x["weight"])
            ranked=tuple((directed+undirected)[:H])
        else:
            neigh.sort(key=lambda x:x["weight"])
            ranked=tuple(neigh[:H])
        self._ranked_neighbor_cache[key]=ranked
        return ranked

    def _neighbor_lookup(self, gene: str) -> Dict[str, Dict]:
        min_score=int(self.adv.get("min_string_score",700)) or 700
        key=(gene,self.api.network_type,min_score)
        cached=self._neighbor_lookup_cache.get(key)
        if cached is not None:
            self._perf_stats["neighbor_lookup_hits"] += 1
            return cached
        self._perf_stats["neighbor_lookup_misses"] += 1
        lookup={}
        for info in self._neighbors(gene):
            g=info.get("gene")
            if g and g not in lookup: lookup[g]=info
        self._neighbor_lookup_cache[key]=lookup
        return lookup

    def _segment_memo_get(self,key):
        val=self._segment_memo.get(key)
        if val is None:
            self._perf_stats["segment_memo_misses"] += 1
            return None
        self._perf_stats["segment_memo_hits"] += 1
        self._segment_memo.move_to_end(key)
        return val

    def _segment_memo_put(self,key,value):
        self._segment_memo[key]=value
        self._segment_memo.move_to_end(key)
        while len(self._segment_memo)>self._segment_memo_max:
            self._segment_memo.popitem(last=False)

    def find_first_interval_baseline_segments(
            self, start_gene: str, end_gene: str, baseline_genes: Set[str],
            blocked_genes: Set[str], bridge_blocked_genes: Optional[Set[str]] = None,
            max_paths_override: Optional[int] = None) -> List[List[Tuple[str, float]]]:
        """First-interval direct/direct baseline bridge screen with no H pre-truncation.

        Baseline detection is availability evidence only, never activation. A baseline
        protein b can serve as one bridge only if start->b and b->target both satisfy
        the configured STRING threshold, Q restrictions, AND the cell-context bridge
        filter. The bridge counts toward D, so D>=1 is required.
        """
        if self.interrupt or self.max_depth < 1 or not start_gene or not end_gene:
            return []
        if self._eq_relaxed(start_gene,end_gene): return []
        hard_blocked=set(blocked_genes or set())
        bridge_blocked=set(bridge_blocked_genes or set())
        baseline=set(baseline_genes or set())
        if not baseline or end_gene in hard_blocked: return []
        is_regulatory=(self.api.network_type == "regulatory")
        reg_bonus=float(self.adv.get("regulatory_weight_bonus",0.1))
        hubs=set(self.adv.get("hub_blacklist",[]))

        # Cache first-hop candidates after all bridge-eligibility gates, including
        # cell context. Baseline MS detection itself is valid cell-presence evidence;
        # experimentally protected genes therefore remain eligible by policy.
        cell_key = None
        if self.cell_context and self.cell_context.enabled:
            cell_key=(self.cell_context.mode,self.cell_context.cell_line,
                      len(self.cell_context.rna_present),len(self.cell_context.protein_present))
        first_key=(start_gene,self.api.network_type,int(self.adv.get("min_string_score",700)) or 700,
                   tuple(sorted(baseline)),tuple(sorted(hard_blocked)),tuple(sorted(bridge_blocked)),
                   tuple(sorted(hubs)),float(reg_bonus),cell_key)
        eligible=self._baseline_first_hop_cache.get(first_key)
        if eligible is None:
            tmp=[]
            for info1 in self._neighbors(start_gene):
                b=info1.get("gene")
                if not b or b not in baseline or b==start_gene: continue
                if b in hubs or b in hard_blocked or b in bridge_blocked: continue
                if self.cell_context and self.cell_context.enabled and not self.cell_context.bridge_allowed(b):
                    continue
                if is_regulatory and (not info1.get("is_directed") or
                    info1.get("source","").casefold()!=start_gene.casefold() or
                    info1.get("target","").casefold()!=b.casefold()):
                    continue
                w1=info1.get("weight",float('inf'))
                if is_regulatory: w1=max(0.01,w1-reg_bonus)
                tmp.append((b,info1,w1))
            eligible=tuple(tmp)
            self._baseline_first_hop_cache[first_key]=eligible

        candidates=[]
        for b,info1,w1 in eligible:
            if self._eq_relaxed(b,end_gene): continue
            lookup=self._neighbor_lookup(b)
            info2=lookup.get(end_gene)
            if info2 is None:
                for cand in lookup.values():
                    g=cand.get("gene")
                    if g and self._eq_relaxed(g,end_gene): info2=cand; break
            if info2 is None: continue
            t=info2.get("gene")
            if not t or t in hard_blocked: continue
            if is_regulatory and (not info2.get("is_directed") or
                info2.get("source","").casefold()!=b.casefold() or
                info2.get("target","").casefold()!=t.casefold()):
                continue
            w2=info2.get("weight",float('inf'))
            if is_regulatory: w2=max(0.01,w2-reg_bonus)
            if is_regulatory:
                for a,z,inf in ((start_gene,b,info1),(b,t,info2)):
                    self.regulatory_info[(a,z)]={"source":a,"target":z,"action":inf.get("action",""),
                        "directionality":inf.get("directionality","unknown"),"is_directed":inf.get("is_directed",False),
                        "score_A_to_B":inf.get("score_A_to_B",0),"score_B_to_A":inf.get("score_B_to_A",0)}
            candidates.append([(start_gene,0.0),(b,w1),(t,w2)])
        unique={}
        for path in candidates:
            k=tuple(g for g,_ in path); c=sum(x for _,x in path)
            if k not in unique or c < sum(x for _,x in unique[k]): unique[k]=path
        out=list(unique.values())
        out.sort(key=lambda p:(sum(c for _,c in p),'->'.join(g for g,_ in p)))
        return out

    def find_path_segment(self, start_gene: str, end_gene: str, blocked_genes: Set[str],
                          bridge_blocked_genes: Optional[Set[str]] = None,
                          max_paths_override: Optional[int] = None) -> List[List[Tuple[str, float]]] :
        """
        Find candidate prior-network routes between two measured anchors.

        Stage 3 bridge-depth semantics:
            D = maximum number of INTERMEDIATE bridge nodes per interval.
            D=0: A -> B (1 edge)
            D=1: A -> B or A -> K1 -> B (<=2 edges)
            D=2: A -> ... -> K1 -> K2 -> B (<=3 edges)

        The target is tested before a neighbor is counted as a bridge.  D resets for
        every temporal interval.  A search-global visited set is intentionally not
        used, so alternative convergent paths are not erased; path-local cycle
        prevention and a per-layer beam cap control combinatorial growth.
        """
        if self.interrupt or not start_gene or not end_gene:
            return []
        if self._eq_relaxed(start_gene, end_gene):
            return [[(start_gene, 0.0), (end_gene, 0.0)]]

        hard_blocked = set(blocked_genes or set())
        bridge_blocked = set(bridge_blocked_genes or set())

        found: List[List[Tuple[str, float]]] = []
        H = int(self.adv.get("degree_cap", 50)) or 50
        per_cap = int(max_paths_override if max_paths_override is not None else self.adv.get("max_paths_per_target", 300)) or 1
        segment_beam = int(self.adv.get("beam_width", 150)) or 150
        verbose = bool(self.adv.get("verbose", False))

        # Include cell-context identity in the memo key. bridge_allowed() is deterministic
        # within one run, so identical source/target/blocker searches can be safely reused.
        cell_key=None
        if self.cell_context and self.cell_context.enabled:
            cell_key=(self.cell_context.mode,self.cell_context.cell_line,
                      len(self.cell_context.rna_present),len(self.cell_context.rna_negative),
                      len(self.cell_context.protein_present),len(self.cell_context.protein_negative))
        memo_key=(start_gene,end_gene,tuple(sorted(hard_blocked)),tuple(sorted(bridge_blocked)),
                  self.max_depth,H,per_cap,segment_beam,self.api.network_type,
                  int(self.adv.get("min_string_score",700)) or 700,
                  tuple(sorted(set(self.adv.get("hub_blacklist",[])))),
                  bool(self.adv.get("no_immediate_backtrack",True)),bool(self.adv.get("dedupe_gene_in_path",True)),
                  float(self.adv.get("regulatory_weight_bonus",0.1)),cell_key)
        memoized=self._segment_memo_get(memo_key)
        if memoized is not None: return memoized
        if verbose:
            self.log(f"    [BFS Start] {start_gene} → {end_gene}  (Bridge Depth D={self.max_depth}, max edges={self.max_depth + 1})")
        bfs_start_time = time.time()
        is_regulatory = (self.api.network_type in ("regulatory",))
        reg_bonus = float(self.adv.get("regulatory_weight_bonus", 0.1))
        hubs = set(self.adv.get("hub_blacklist", []))

        # frontier item = (node, path, cumulative_cost), where bridge_count is the layer index.
        frontier = [(start_gene, [(start_gene, 0.0)], 0.0)]
        expanded_total = 0

        for bridge_count in range(self.max_depth + 1):
            if not frontier or self.interrupt:
                break
            if verbose:
                self.log(f"    [BFS Layer] bridge_count={bridge_count}, frontier={len(frontier)}, found={len(found)}")

            next_frontier = []
            for node, path, cum_cost in frontier:
                if self.interrupt:
                    break
                neigh = self._ranked_neighbors(node)
                expanded_total += 1

                genes_in_path = {g for g, _ in path}
                prev_gene = path[-2][0] if len(path) >= 2 else None

                for nxt_info in neigh:
                    nxt = nxt_info["gene"]
                    edge_weight = nxt_info["weight"]
                    if is_regulatory:
                        # Directed regulatory prior: only regulator -> target may be traversed.
                        if not nxt_info.get("is_directed") or nxt_info.get("source", "").upper() != node.upper() or nxt_info.get("target", "").upper() != nxt.upper():
                            continue
                        edge_weight = max(0.01, edge_weight - reg_bonus)

                    # Explicit inhibitors are hard blockers, including as endpoints.
                    if nxt in hard_blocked:
                        continue
                    if self.adv.get("no_immediate_backtrack", True) and prev_gene and nxt == prev_gene:
                        continue
                    if self.adv.get("dedupe_gene_in_path", True) and nxt in genes_in_path:
                        continue

                    if is_regulatory:
                        self.regulatory_info[(node, nxt)] = {
                            "source": node, "target": nxt,
                            "action": nxt_info.get("action", ""),
                            "directionality": nxt_info.get("directionality", "unknown"),
                            "is_directed": nxt_info.get("is_directed", False),
                            "score_A_to_B": nxt_info.get("score_A_to_B", 0),
                            "score_B_to_A": nxt_info.get("score_B_to_A", 0)
                        }

                    newp = path + [(nxt, edge_weight)]
                    new_cost = cum_cost + edge_weight

                    # IMPORTANT: target is checked before bridge counting/blocking.
                    if self._eq_relaxed(nxt, end_gene):
                        found.append(newp)
                        if verbose:
                            self.log(f"    [FOUND D={bridge_count}] {' -> '.join(g for g, _ in newp)} (cost: {new_cost:.3f})")
                        if len(found) >= per_cap:
                            found.sort(key=lambda p: (sum(c for _, c in p), '->'.join(g for g, _ in p)))
                            result=found[:per_cap]
                            if not self.interrupt: self._segment_memo_put(memo_key,result)
                            return result
                        continue

                    # Non-target neighbor becomes one additional bridge node.
                    if nxt in bridge_blocked:
                        continue
                    # Cell-context filter is deliberately BRIDGE-ONLY.
                    # Experimental target anchors were already accepted above and are never
                    # removed because a bulk RNA/proteome file lacks them.
                    if self.cell_context and self.cell_context.enabled:
                        if not self.cell_context.bridge_allowed(nxt):
                            continue
                    if bridge_count >= self.max_depth:
                        continue
                    next_frontier.append((nxt, newp, new_cost))

            # Keep multiple convergent alternatives, but cap the state frontier by cost.
            if next_frontier:
                next_frontier.sort(key=lambda x: (x[2], '->'.join(g for g, _ in x[1])))
                frontier = next_frontier[:segment_beam]
            else:
                frontier = []

        found.sort(key=lambda p: (sum(c for _, c in p), '->'.join(g for g, _ in p)))
        bfs_elapsed = time.time() - bfs_start_time
        result=found[:per_cap]
        if not self.interrupt: self._segment_memo_put(memo_key,result)
        if verbose or bfs_elapsed >= 1.0:
            self.log(f"    [BFS Done] {start_gene} → {end_gene}: {len(found)} paths, expanded={expanded_total} ({bfs_elapsed:.2f} sec)")
        return result

    def trace_full_pathway(self, start_proteins, input_data, inhibitor_drugs=None):
        adapter = TemporalInputAdapter(input_data if isinstance(input_data, dict) else {})
        self.input_schema_label = adapter.schema_label
        self.log(f"[Input Schema] {adapter.schema_label}")
        if adapter.is_temporal_v2:
            return self._trace_temporal(start_proteins, input_data, inhibitor_drugs or [])
        ts = (input_data or {}).get("time_series", {}) if isinstance(input_data, dict) else (input_data or {})
        return self._trace_legacy(start_proteins, ts, inhibitor_drugs or [])

    def _trace_temporal(self, start_proteins, temporal_data, inhibitor_drugs):
        """FINAL temporal tracing with a Q-derived condition-specific functional pool.

        Key principles:
          * No gene is excluded because of a drug/condition name.
          * Pretreatment-vs-reference baseline_contrast is interpreted through Stage 2.1 Q states.
          * A gene is excluded from the starting functional pool only when it is opposing-only
            (>=1 response_opposing site and no response_supporting site).
          * Mixed supporting/opposing genes remain eligible.
          * Within each interval, opposing-only genes are not active targets/bridges/sources.
          * Alternative source-target routes are preserved up to the original BFS limits
            (D, STRING threshold, degree cap H, beam width K, max paths/target N).
          * A target-preserving adaptive temporal beam guarantees >=1 route per reached target.
        """
        is_regulatory = (self.api.network_type in ("regulatory",))
        mode_str = f" ({self.api.network_type} mode)" if is_regulatory else ""
        anchor_mode = str(self.adv.get("temporal_anchor_mode", "measured"))
        exclude_opp_bridge = bool(self.adv.get("exclude_opposing_as_bridge", True))
        intervals = temporal_data.get("intervals", []) or []
        self.temporal_path_records = []
        self.temporal_interval_summary = []
        self.temporal_pruning_summary = []

        self.log(f"=== Start TEMPORAL interval pathway tracing{mode_str} ===")
        self.log(f"  [Temporal] Target mode: {anchor_mode}")
        self.log(f"  [Temporal] Bridge Depth D={self.max_depth} means intermediate bridge nodes; max edges/interval={self.max_depth + 1}")
        self.log("  [Temporal] Baseline availability is NOT treated as activation.")
        self.log("  [Temporal] D resets independently for every interval.")
        self.log("  [FINAL] No drug-name hard blocking; target-preserving adaptive temporal beam enabled.")

        adapter = TemporalInputAdapter(temporal_data)
        all_names: Set[str] = set(_as_name(x).lstrip('-') for x in (start_proteins or []))
        for tok in adapter.all_site_tokens():
            all_names.add(_as_name(tok).lstrip('-'))
        for inh in inhibitor_drugs or []:
            all_names.add(_as_name(inh).lstrip('-'))

        name_to_gene = {}
        names_to_map = sorted(x for x in all_names if x)
        n_names_to_map = max(1, len(names_to_map))
        self.report_progress(
            0.0,
            "Mapping protein/site identifiers",
            0,
            len(names_to_map)
        )
        for map_idx, n in enumerate(names_to_map, 1):
            if self.interrupt:
                break
            g = self._gene(n)
            name_to_gene[n] = g if g else "None"
            if (
                map_idx == 1
                or map_idx == len(names_to_map)
                or map_idx % max(1, n_names_to_map // 100) == 0
            ):
                self.report_progress(
                    0.08 * (map_idx / n_names_to_map),
                    "Mapping protein/site identifiers",
                    map_idx,
                    len(names_to_map)
                )

        if self.cell_context and self.cell_context.enabled:
            experimental_genes = set()
            for sp in start_proteins or []:
                g = self._gene(sp)
                if g: experimental_genes.add(g)
            for tok in adapter.all_site_tokens():
                g = self._gene(_as_name(tok))
                if g: experimental_genes.add(g)
            self.cell_context.set_experimental_genes(experimental_genes)

        # Baseline availability is used only as first-interval bridge evidence.
        # It is never interpreted as activation or as an extra temporal layer.
        baseline_tokens = (temporal_data.get("baseline") or {}).get("available", []) or []
        baseline_map = TemporalInputAdapter.tokens_to_gene_map(baseline_tokens, self._gene)
        baseline_genes = set(baseline_map)
        self.log(
            f"  [Baseline Bridge Pool] {len(baseline_genes)} baseline-available genes loaded; "
            "used only as first-interval bridge candidates (not activation targets)."
        )

        # Explicit/manual blockers remain supported for backward compatibility; the FINAL model
        # never infers a blocker from the condition/drug name.
        hard_blocked: Set[str] = set()
        for x in inhibitor_drugs or []:
            g = self._gene(x)
            if g:
                hard_blocked.add(g)

        pool_info = adapter.condition_functional_pool(self._gene)
        condition_excluded = set(pool_info["excluded_low_activity"])
        condition_mixed = set(pool_info["mixed"])
        condition_supporting = set(pool_info["supporting_only"])
        condition_unresolved = set(pool_info["unresolved_only"])
        self.condition_functional_pool_info = {
            "basis": "baseline_contrast_Q_opposing_only",
            "has_baseline_contrast": bool(pool_info["has_baseline_contrast"]),
            "baseline_contrast_type": pool_info["baseline_contrast_type"],
            "excluded_low_activity_genes": sorted(condition_excluded),
            "supporting_only_genes": sorted(condition_supporting),
            "mixed_genes": sorted(condition_mixed),
            "unresolved_only_genes": sorted(condition_unresolved),
            "drug_name_inference_used": False,
            "rule": "exclude gene only if baseline_contrast has response_opposing site(s) and no response_supporting site"
        }
        if pool_info["has_baseline_contrast"]:
            self.log(f"  [Q-Pool] Pretreatment baseline contrast detected ({pool_info['baseline_contrast_type'] or 'unspecified'}).")
            self.log(f"  [Q-Pool] supporting-only={len(condition_supporting)}, mixed={len(condition_mixed)}, unresolved-only={len(condition_unresolved)}, excluded opposing-only={len(condition_excluded)}")
            if condition_excluded:
                preview = ', '.join(sorted(condition_excluded)[:25])
                more = f" ... (+{len(condition_excluded)-25} more)" if len(condition_excluded) > 25 else ""
                self.log(f"  [Q-Pool] Excluded from condition-specific functional pool: {preview}{more}")
        else:
            self.log("  [Q-Pool] No pretreatment baseline contrast; no condition-specific Q exclusion applied.")

        start_genes = []
        for sp in start_proteins or []:
            g = self._gene(sp)
            if g and g not in start_genes:
                start_genes.append(g)
        if not start_genes:
            self.log("[Error] Failed to map the start proteins.")
            return [], name_to_gene, []
        if not intervals:
            self.log("[Error] Temporal v2 input contains no intervals.")
            return [], name_to_gene, []

        # Explicit experimental source takes precedence over pool exclusion; warn if a source
        # itself is classified opposing-only at pretreatment.
        for sg in start_genes:
            if sg in condition_excluded:
                self.log(f"  [Q-Pool Warning] Explicit start protein {sg} is opposing-only at baseline contrast; retained ONLY to initiate interval 1 because it is an explicit experimental source. It remains excluded as an ordinary bridge/target/later source.")

        all_complete: Set[Tuple[str, ...]] = set()
        all_intermediate: Set[Tuple[str, ...]] = set()
        final_records = []

        total_temporal_units = max(
            1, len(start_genes) * len(intervals)
        )

        for sidx, sg in enumerate(start_genes, 1):
            if self.interrupt:
                break
            self.log(f"--- Temporal Start {sidx}/{len(start_genes)}: {sg} ---")
            states = [{"path": [(sg, 0.0)], "segments": [], "source_anchor": sg}]
            survived_all = True

            for iidx, interval in enumerate(intervals):
                if self.interrupt:
                    survived_all = False
                    break
                iid = str(interval.get("id") or f"{interval.get('from','?')}-{interval.get('to','?')}")
                is_last = (iidx == len(intervals) - 1)

                measured = interval.get("response_anchors")
                if measured is None:
                    measured = list(interval.get("increased", []) or []) + list(interval.get("decreased", []) or [])
                supporting = interval.get("response_supporting", []) or []
                opposing = interval.get("response_opposing", []) or []
                unresolved = interval.get("response_unresolved", []) or []

                measured_map = TemporalInputAdapter.tokens_to_gene_map(measured, self._gene)
                supporting_map = TemporalInputAdapter.tokens_to_gene_map(supporting, self._gene)
                opposing_map = TemporalInputAdapter.tokens_to_gene_map(opposing, self._gene)
                unresolved_map = TemporalInputAdapter.tokens_to_gene_map(unresolved, self._gene)

                measured_genes = set(measured_map)
                supporting_genes = set(supporting_map)
                opposing_genes = set(opposing_map)
                unresolved_genes = set(unresolved_map)
                interval_opposing_only = opposing_genes - supporting_genes
                interval_mixed = opposing_genes & supporting_genes

                # Target construction preserves measured evidence but removes genes whose
                # Q evidence indicates low functional state. Mixed genes remain eligible.
                if anchor_mode == "strict_signed":
                    targets = set(supporting_genes)
                else:
                    targets = set(measured_genes) - interval_opposing_only
                targets -= condition_excluded
                target_map = {g: measured_map.get(g, []) for g in targets}

                self.log(f"--- Interval {iidx+1}/{len(intervals)}: {iid} ---")
                self.log(f"  Measured={len(measured_genes)}, supporting={len(supporting_genes)}, opposing={len(opposing_genes)}, unresolved={len(unresolved_genes)}, mixed={len(interval_mixed)}")
                self.log(f"  Q-low opposing-only this interval={len(interval_opposing_only)}; condition-pool excluded={len(condition_excluded)}")
                self.log(f"  Eligible targets ({anchor_mode})={len(targets)}, incoming source paths={len(states)}")

                if not targets:
                    self.log("  [Temporal Stop] No eligible functional anchors in this interval; temporal chain cannot advance.")
                    survived_all = False
                    break
                if not states:
                    survived_all = False
                    break

                next_states = []
                reached_genes = set()
                endpoint_records = []
                pair_count = len(states) * len(targets)
                processed_pairs = 0
                pair_progress_step = max(1, pair_count // 100)
                temporal_unit_index = (
                    (sidx - 1) * len(intervals) + iidx
                )
                self.log(f"  [Estimate] {len(states)} source paths x {len(targets)} target genes = {pair_count} pairs")

                base_search_fraction = (
                    temporal_unit_index / total_temporal_units
                )
                self.report_progress(
                    0.08 + 0.86 * base_search_fraction,
                    f"Temporal interval {iidx + 1}/{len(intervals)}: {iid}",
                    0,
                    pair_count
                )

                # Q-derived pool restrictions: condition baseline low-state genes are excluded
                # throughout the run; current interval opposing-only genes cannot be bridges.
                bridge_blocked = set(condition_excluded)
                if exclude_opp_bridge:
                    bridge_blocked |= interval_opposing_only
                effective_hard_blocked = set(hard_blocked) | set(condition_excluded)

                for state in states:
                    p = state["path"]
                    if len(p) > 1:
                        all_intermediate.add(tuple(g for g, _ in p))
                    source_gene = p[-1][0]
                    # Explicit experimental start proteins are allowed to initiate the FIRST
                    # interval even when pretreatment Q classifies them as opposing-only.
                    # This is essential for ligand-driven experiments such as EGF -> EGFR:
                    # the experiment explicitly imposes the source, while Q-derived pool
                    # exclusion still prevents that gene from being used later as an
                    # ordinary bridge/target/propagating source.
                    is_explicit_initial_source = (iidx == 0 and any(self._eq_relaxed(source_gene, x) for x in start_genes))
                    if source_gene in condition_excluded and not is_explicit_initial_source:
                        continue
                    prior_full_path_genes = {g for g, _ in p[:-1]}
                    full_path_dedupe = bool(
                        self.adv.get("temporal_full_path_dedupe", True)
                    )
                    temporal_pair_cap = max(
                        1,
                        int(self.adv.get("temporal_segment_route_cap", 5)) or 1
                    )

                    for tgt in targets:
                        processed_pairs += 1

                        if (
                            processed_pairs == 1
                            or processed_pairs == pair_count
                            or processed_pairs % pair_progress_step == 0
                        ):
                            interval_fraction = (
                                processed_pairs / max(1, pair_count)
                            )
                            search_fraction = (
                                temporal_unit_index + interval_fraction
                            ) / total_temporal_units
                            self.report_progress(
                                0.08 + 0.86 * search_fraction,
                                (
                                    f"Temporal interval "
                                    f"{iidx + 1}/{len(intervals)}: {iid}"
                                ),
                                processed_pairs,
                                pair_count
                            )

                        path_specific_blocked = set(effective_hard_blocked)

                        # Prevent re-entry into genes already traversed in earlier
                        # temporal intervals.  Current source is excluded by p[:-1].
                        # Same-anchor persistence remains allowed.
                        if (
                            full_path_dedupe
                            and not self._eq_relaxed(source_gene, tgt)
                        ):
                            path_specific_blocked |= prior_full_path_genes

                        segs = self.find_path_segment(
                            source_gene,
                            tgt,
                            path_specific_blocked,
                            bridge_blocked,
                            max_paths_override=temporal_pair_cap
                        )

                        # FIRST INTERVAL ONLY: separately screen all baseline-available
                        # bridge candidates without H pre-truncation. Both direct STRING
                        # legs, Q eligibility and cell-context eligibility must pass.
                        baseline_segs=[]
                        if iidx == 0 and self.max_depth >= 1:
                            baseline_segs=self.find_first_interval_baseline_segments(
                                source_gene,tgt,baseline_genes,path_specific_blocked,
                                bridge_blocked,max_paths_override=temporal_pair_cap)

                        merged={}
                        for _seg in list(segs)+list(baseline_segs):
                            _k=tuple(g for g,_ in _seg); _c=sum(c for _,c in _seg)
                            if _k not in merged or _c < sum(c for _,c in merged[_k]): merged[_k]=_seg
                        _all=list(merged.values())
                        _all.sort(key=lambda q:(sum(c for _,c in q),'->'.join(g for g,_ in q)))
                        if iidx == 0:
                            _mandatory=[]; _mandatory_keys=set()
                            for _seg in baseline_segs:
                                if len(_seg)==3 and _seg[1][0] in measured_genes:
                                    _k=tuple(g for g,_ in _seg)
                                    if _k not in _mandatory_keys:
                                        _mandatory.append(_seg); _mandatory_keys.add(_k)
                            _mandatory.sort(key=lambda q:(sum(c for _,c in q),'->'.join(g for g,_ in q)))
                            _selected=list(_mandatory); _selected_keys={tuple(g for g,_ in q) for q in _selected}
                            _effective=max(temporal_pair_cap,len(_mandatory))
                            for _seg in _all:
                                _k=tuple(g for g,_ in _seg)
                                if _k in _selected_keys: continue
                                if len(_selected)>=_effective: break
                                _selected.append(_seg); _selected_keys.add(_k)
                            segs=_selected
                        else:
                            segs=_all[:temporal_pair_cap]

                        for seg in segs:
                            if not seg or len(seg) < 2:
                                continue
                            same_anchor = self._eq_relaxed(source_gene, tgt)
                            combined = p if same_anchor else p + seg[1:]
                            bridge_nodes = [] if same_anchor else [g for g, _ in seg[1:-1]]
                            segment_record = {
                                "interval": iid,
                                "from": interval.get("from"),
                                "to": interval.get("to"),
                                "source_anchor": source_gene,
                                "target_anchor": tgt,
                                "target_anchor_tokens": measured_map.get(tgt, []),
                                "target_role": ("mixed" if tgt in interval_mixed else
                                                "supporting" if tgt in supporting_genes else
                                                "unresolved" if tgt in unresolved_genes else "measured"),
                                "bridge_nodes": bridge_nodes,
                                "baseline_supported_bridge": bool(
                                    iidx == 0 and len(bridge_nodes) == 1 and bridge_nodes[0] in baseline_genes
                                ),
                                "baseline_bridge_measured_in_interval": bool(
                                    iidx == 0 and len(bridge_nodes) == 1 and bridge_nodes[0] in baseline_genes
                                    and bridge_nodes[0] in measured_genes
                                ),
                                "bridge_depth": len(bridge_nodes),
                                "edge_count": 0 if same_anchor else len(seg) - 1,
                                "segment_cost": 0.0 if same_anchor else round(sum(c for _, c in seg[1:]), 6)
                            }
                            new_segments = state["segments"] + [segment_record]
                            ns = {"path": combined, "segments": new_segments, "source_anchor": tgt}
                            endpoint_records.append(ns)
                            reached_genes.add(tgt)
                            # Because opposing-only targets were removed before searching, every
                            # reached target is eligible to seed the next interval.
                            next_states.append(ns)

                # FINAL: exact gene-chain dedupe precedes target-preserving adaptive
                # temporal beam selection.
                unique = {}
                for st in next_states:
                    key = tuple(g for g, _ in st["path"])
                    cost = sum(c for _, c in st["path"])
                    if key not in unique or cost < sum(c for _, c in unique[key]["path"]):
                        unique[key] = st
                next_states = list(unique.values())
                next_states.sort(
                    key=lambda st: (
                        sum(c for _, c in st["path"]),
                        '->'.join(g for g, _ in st["path"])
                    )
                )

                exact_unique_count = len(next_states)

                # Base temporal beam.  0 means reuse the segment Beam K
                # only as the BASE temporal width.
                configured_temporal_k = int(
                    self.adv.get("temporal_beam_width", 0) or 0
                )
                if configured_temporal_k <= 0:
                    configured_temporal_k = (
                        int(self.adv.get("beam_width", 150)) or 150
                    )

                per_target_cap = max(
                    1,
                    int(self.adv.get("temporal_routes_per_target", 3)) or 1
                )

                def _state_key(st):
                    return (
                        sum(c for _, c in st["path"]),
                        '->'.join(g for g, _ in st["path"])
                    )

                # Group by the reached target that seeds the next interval.
                # Each non-empty group is one actually reached target gene.
                by_target = {}
                for st in next_states:
                    by_target.setdefault(st["source_anchor"], []).append(st)
                for group in by_target.values():
                    group.sort(key=_state_key)

                post_dedup_reached_target_count = len(by_target)

                # First-interval measured-baseline bridge histories are protected
                # through temporal beam pruning, while baseline remains non-activating.
                mandatory_baseline_states=[]
                if iidx == 0:
                    for st in next_states:
                        _segs=st.get("segments",[]) or []
                        if _segs and _segs[-1].get("baseline_bridge_measured_in_interval",False):
                            mandatory_baseline_states.append(st)
                    mandatory_baseline_states.sort(key=_state_key)

                # Preserve mandatory measured-baseline histories AND guarantee
                # one representative for every reached target.  The effective
                # temporal beam must therefore be large enough for the UNION of
                # those two protected sets, not merely the maximum of their
                # individual sizes.
                representatives=[group[0] for group in by_target.values() if group]
                representatives.sort(key=_state_key)

                _required_ids = {id(st) for st in mandatory_baseline_states}
                _required_ids.update(id(st) for st in representatives)
                protected_union_count = len(_required_ids)

                effective_temporal_k = max(
                    configured_temporal_k,
                    protected_union_count
                )
                temporal_beam_expanded_for_targets = (
                    effective_temporal_k > configured_temporal_k
                )

                selected=list(mandatory_baseline_states)
                selected_ids={id(st) for st in selected}
                for _rep in representatives:
                    if id(_rep) not in selected_ids:
                        selected.append(_rep); selected_ids.add(id(_rep))

                # Preferentially retain additional histories per target
                # according to the target route quota Q.
                alternatives = []
                for group in by_target.values():
                    for st in group[1:per_target_cap]:
                        if id(st) not in selected_ids:
                            alternatives.append(st)
                alternatives.sort(key=_state_key)
                selected.extend(
                    alternatives[
                        :max(0, effective_temporal_k - len(selected))
                    ]
                )

                # Fill any remaining capacity globally by cumulative cost.
                if len(selected) < effective_temporal_k:
                    selected_ids = {id(st) for st in selected}
                    leftovers = [
                        st for st in next_states
                        if id(st) not in selected_ids
                    ]
                    leftovers.sort(key=_state_key)
                    selected.extend(
                        leftovers[
                            :effective_temporal_k - len(selected)
                        ]
                    )

                next_states = selected[:effective_temporal_k]
                next_states.sort(key=_state_key)

                propagated_target_count = len(
                    {st["source_anchor"] for st in next_states}
                )
                all_reached_targets_preserved = (
                    propagated_target_count
                    == post_dedup_reached_target_count
                )

                if not all_reached_targets_preserved:
                    raise RuntimeError(
                        "Target-preserving temporal beam invariant failed: "
                        f"{propagated_target_count}/"
                        f"{post_dedup_reached_target_count} reached targets "
                        "were preserved."
                    )

                if exact_unique_count != len(endpoint_records):
                    self.log(
                        f"  [Exact Dedupe] "
                        f"{len(endpoint_records)} candidate paths -> "
                        f"{exact_unique_count} unique gene chains."
                    )

                if temporal_beam_expanded_for_targets:
                    self.log(
                        f"  [Adaptive Temporal Beam] Base K="
                        f"{configured_temporal_k} expanded to "
                        f"K_eff={effective_temporal_k} because "
                        f"{post_dedup_reached_target_count} target genes "
                        f"were reached and {protected_union_count} protected "
                        f"target/baseline histories were required."
                    )

                self.log(
                    f"  [Target-Preserving Temporal Beam] "
                    f"{exact_unique_count} unique routes -> "
                    f"{len(next_states)} propagated states; "
                    f"targets preserved="
                    f"{propagated_target_count}/"
                    f"{post_dedup_reached_target_count}; "
                    f"K_base={configured_temporal_k}, "
                    f"K_eff={effective_temporal_k}, "
                    f"protected_union={protected_union_count}, "
                    f"routes/target quota={per_target_cap}, "
                    f"segment routes/pair="
                    f"{int(self.adv.get('temporal_segment_route_cap', 5))}."
                )

                self.temporal_pruning_summary.append({
                    "start_protein": sg,
                    "interval": iid,
                    "candidate_endpoint_paths": len(endpoint_records),
                    "after_exact_gene_path_dedup": exact_unique_count,
                    "after_target_preserving_temporal_beam": len(next_states),
                    "best_per_source_target_pruning": False,
                    "best_per_reached_target_pruning": False,
                    "target_stratified_interval_beam": True,
                    "target_preserving_adaptive_temporal_beam": True,
                    "fixed_interval_beam_pruning": False,
                    "configured_temporal_beam_width":
                        configured_temporal_k,
                    "effective_temporal_beam_width":
                        effective_temporal_k,
                    "adaptive_beam_expanded_for_targets":
                        temporal_beam_expanded_for_targets,
                    "reached_target_genes_after_exact_dedup":
                        post_dedup_reached_target_count,
                    "propagated_target_genes":
                        propagated_target_count,
                    "all_reached_targets_preserved":
                        all_reached_targets_preserved,
                    "temporal_routes_per_target": int(
                        self.adv.get("temporal_routes_per_target", 3)
                    ),
                    "temporal_full_path_dedupe": bool(
                        self.adv.get("temporal_full_path_dedupe", True)
                    ),
                    "segment_route_limit":
                        "temporal_segment_route_cap per source-target pair",
                    "temporal_segment_route_cap": int(
                        self.adv.get("temporal_segment_route_cap", 5)
                    ),
                    "legacy_max_paths_per_target_N": int(
                        self.adv.get("max_paths_per_target", 300)
                    )
                })

                self.temporal_interval_summary.append({
                    "start_protein": sg,
                    "interval": iid,
                    "from": interval.get("from"),
                    "to": interval.get("to"),
                    "anchor_mode": anchor_mode,
                    "measured_anchor_genes": len(measured_genes),
                    "eligible_target_genes": len(targets),
                    "supporting_genes": len(supporting_genes),
                    "opposing_genes": len(opposing_genes),
                    "opposing_only_genes": len(interval_opposing_only),
                    "mixed_genes": len(interval_mixed),
                    "unresolved_genes": len(unresolved_genes),
                    "condition_pool_excluded_genes": len(condition_excluded),
                    "reached_target_genes": len(reached_genes),
                    "candidate_endpoint_paths": len(endpoint_records),
                    "after_exact_gene_path_dedup": exact_unique_count,
                    "reached_target_genes_after_exact_dedup":
                        post_dedup_reached_target_count,
                    "propagating_unique_paths": len(next_states),
                    "propagated_target_genes":
                        propagated_target_count,
                    "all_reached_targets_preserved":
                        all_reached_targets_preserved,
                    "configured_temporal_beam_width":
                        configured_temporal_k,
                    "effective_temporal_beam_width":
                        effective_temporal_k,
                    "adaptive_beam_expanded_for_targets":
                        temporal_beam_expanded_for_targets,
                    "bridge_depth_D": self.max_depth
                })

                if not endpoint_records:
                    self.log(f"  [Temporal Stop] No route reached any eligible anchor in {iid}.")
                    survived_all = False
                    break
                if not is_last and not next_states:
                    self.log(f"  [Temporal Stop] Anchors were reached in {iid}, but none can propagate to the next interval.")
                    survived_all = False
                    break

                self.log(f"  [Interval Done] reached genes={len(reached_genes)}, endpoint candidates={len(endpoint_records)}, unique propagating routes={len(next_states)}")
                self.log(
                    "  [Speed Cache] "
                    f"ranked-neighbor hit/miss={self._perf_stats['ranked_neighbor_hits']}/{self._perf_stats['ranked_neighbor_misses']}, "
                    f"direct-lookup hit/miss={self._perf_stats['neighbor_lookup_hits']}/{self._perf_stats['neighbor_lookup_misses']}, "
                    f"segment hit/miss={self._perf_stats['segment_memo_hits']}/{self._perf_stats['segment_memo_misses']}"
                )
                states = next_states

            if survived_all and states:
                for st in states:
                    genes = tuple(g for g, _ in st["path"])
                    if len(genes) > 1:
                        all_complete.add(genes)
                        final_records.append(st)

        # Exact full-path gene-chain deduplication only.
        uniq_final = {}
        for st in final_records:
            key = tuple(g for g, _ in st["path"])
            cost = sum(c for _, c in st["path"])
            if key not in uniq_final or cost < sum(c for _, c in uniq_final[key]["path"]):
                uniq_final[key] = st
        final_records = list(uniq_final.values())
        final_records.sort(key=lambda st: (sum(c for _, c in st["path"]), '->'.join(g for g, _ in st["path"])))

        self.temporal_path_records = []
        for idx, st in enumerate(final_records, 1):
            rec = {
                "path_id": idx,
                "genes": [g for g, _ in st["path"]],
                "total_cost": round(sum(c for _, c in st["path"]), 6),
                "segments": st["segments"]
            }
            self.temporal_path_records.append(rec)

        self.log(f"=== Temporal tracing completed: {len(all_complete)} unique complete paths, {len(all_intermediate)} intermediate paths ===")
        self.report_progress(
            0.94,
            "Temporal search complete; preparing result files",
            len(intervals),
            len(intervals)
        )
        return list(all_complete), name_to_gene, list(all_intermediate)

    def _trace_legacy(self, start_proteins, time_series_data, inhibitor_drugs):
        is_regulatory = (self.api.network_type in ("regulatory",))
        mode_str = f" ({self.api.network_type} mode)" if is_regulatory else ""
        self.log(f"=== Start full pathway tracing (STRING API v12.5 - weighted){mode_str} ===")

        if is_regulatory:
            self.log("  [Regulatory Network] Apply directionality (A->B) and regulatory type (activation/inhibition) information")
            self.log("  [Regulatory Network] Prioritize forward direction from the previous to the next time point")

        # Name -> Gene Mapping
        all_names: Set[str] = set()
        for n in (start_proteins or []):
            all_names.add(_as_name(n).lstrip('-'))

        if isinstance(time_series_data, dict):
            for tp, arr in (time_series_data or {}).items():
                for tok in (arr or []):
                    all_names.add(_as_name(tok).lstrip('-'))

        for inh in (inhibitor_drugs or []):
            all_names.add(_as_name(inh).lstrip('-'))

        name_to_gene = {}
        for n in sorted(all_names):
            g = self._gene(n)
            name_to_gene[n] = g if g else "None"

        self.log("=== Name→Gene Mapping (Summary) ===")
        shown = 0
        for k, v in sorted(name_to_gene.items()):
            self.log(f"  {k} -> {v}")
            shown += 1
            if shown >= 30 and len(name_to_gene) > 30:
                self.log(f"  ... ({len(name_to_gene) - 30} more)")
                break

        if self.cell_context and self.cell_context.enabled:
            experimental_genes = set()
            for sp in start_proteins or []:
                g = self._gene(sp)
                if g: experimental_genes.add(g)
            if isinstance(time_series_data, dict):
                for arr in time_series_data.values():
                    for tok in arr or []:
                        g = self._gene(_as_name(tok).lstrip('-'))
                        if g: experimental_genes.add(g)
            self.cell_context.set_experimental_genes(experimental_genes)

        blocked: Set[str] = set()
        for x in (inhibitor_drugs or []):
            g = self._gene(x)
            if g:
                blocked.add(g)

        start_genes: List[str] = []
        for sp in (start_proteins or []):
            g = self._gene(sp)
            if g:
                start_genes.append(g)
            else:
                self.log(f"[Warning] Failed to map start protein: {sp}")

        if not start_genes:
            self.log("[Error] Failed to map the start proteins.")
            return [], name_to_gene, []

        def tkey(k):
            m = re.search(r'(\d+)', str(k))
            return int(m.group(1)) if m else float('inf')

        tps = sorted((time_series_data or {}).keys(), key=tkey)
        all_paths: Set[Tuple] = set()
        all_intermediate_paths: Set[Tuple] = set()
        K = int(self.adv.get("beam_width", 150)) or 150

        total_legacy_units = max(
            1, len(start_genes) * max(1, len(tps))
        )

        for idx, sg in enumerate(start_genes, 1):
            if self.interrupt:
                break

            self.log(f"--- Start {idx}/{len(start_genes)}: {sg} ---")
            current = [[(sg, 0.0)]]
            finished = []

            self.log(f"  [Init] Start paths: {current}")
            self.log(f"  [Init] Path count: {len(current)}, first path length: {len(current[0])}")

            for step, tp in enumerate(tps, 1):
                if self.interrupt:
                    break

                total_steps = len(tps)
                self.log(f"--- Time point {step}/{total_steps}: {tp} ---")

                # save intermediate paths
                for p in current:
                    if len(p) > 1:
                        path_genes = tuple(g for g, c in p)
                        all_intermediate_paths.add(path_genes)

                arr = time_series_data.get(tp, []) or []
                positives, negatives = [], []

                for tok in arr:
                    s = _as_name(tok).strip()
                    if not s:
                        continue
                    if s.startswith('-'):
                        negatives.append(s[1:])
                    else:
                        positives.append(s)

                targets = set()
                for n in positives:
                    g = self._gene(n)
                    if g:
                        targets.add(g)

                blockers = set(blocked)
                for n in negatives:
                    g = self._gene(n)
                    if g:
                        blockers.add(g)

                self.log(f"  [+] Targets: {len(targets)} / [-] Blocked: {len(blockers)}")

                if not targets:
                    self.log("  No targets -> next time point")
                    continue

                target_list = list(targets)
                target_names = [t for t in target_list[:5]]
                if len(targets) > 5:
                    self.log(f"  Target examples: {', '.join(target_names)}... ({len(targets) - 5} items)")
                else:
                    self.log(f"  Targets: {', '.join(target_names)}")

                total_pairs = len(current) * len(targets)
                legacy_unit_index = (
                    (idx - 1) * max(1, len(tps)) + (step - 1)
                )
                self.log(f"  [Estimate] {len(current)} paths x {len(targets)} targets = {total_pairs} source-target pairs to process")
                self.report_progress(
                    0.94 * (
                        legacy_unit_index / total_legacy_units
                    ),
                    f"Legacy time point {step}/{len(tps)}: {tp}",
                    0,
                    total_pairs
                )
                if total_pairs > 1000:
                    self.log(f"  ⚠️  Large workload expected (approx. {total_pairs // 100} * 10 sec = {total_pairs // 10} sec may be required)")

                next_paths = []
                processed_targets = 0
                total_current_paths = len(current)

                if not current:
                    self.log(f"  [Warning] Current path list is empty -> ending this start point")
                    break

                last_report_time = time.time()

                for p_idx, p in enumerate(current, 1):
                    if not p or not isinstance(p, list):
                        self.log(f"  [Warning] Skipping invalid path format: {p}")
                        continue

                    if not isinstance(p[-1], tuple) or len(p[-1]) != 2:
                        self.log(f"  [Warning] Invalid node format: {p[-1]}, full path: {p}")
                        continue

                    last_gene, _ = p[-1]

                    if total_current_paths > 5 and p_idx % 5 == 0:
                        elapsed = time.time() - last_report_time
                        self.log(f"  [Progress] Path {p_idx}/{total_current_paths} processing... "
                                 f"(elapsed: {elapsed:.1f} sec, paths found: {len(next_paths)} items)")
                        last_report_time = time.time()

                    for tgt_idx, tgt in enumerate(targets, 1):
                        processed_targets += 1

                        legacy_progress_step = max(
                            1, min(50, max(1, total_pairs // 100))
                        )
                        if (
                            processed_targets % legacy_progress_step == 0
                            or processed_targets == total_pairs
                        ):
                            progress_pct = (
                                processed_targets / max(1, total_pairs)
                            ) * 100
                            self.log(f"  [Progress] {processed_targets}/{total_pairs}  pairs processed "
                                     f"({progress_pct:.1f}%), paths found: {len(next_paths)}")

                            within_fraction = (
                                processed_targets / max(1, total_pairs)
                            )
                            overall_fraction = (
                                legacy_unit_index + within_fraction
                            ) / total_legacy_units
                            self.report_progress(
                                0.94 * overall_fraction,
                                f"Legacy time point {step}/{len(tps)}: {tp}",
                                processed_targets,
                                total_pairs
                            )

                        segs = self.find_path_segment(last_gene, tgt, blockers)
                        for seg in segs:
                            if seg and len(seg) > 1:
                                next_paths.append(p + seg[1:])
                            elif seg and len(seg) == 1:
                                self.log(f"    [Warning] Single-node segment: {seg}")

                self.log(f"  [Done] Total {processed_targets} source-target pairs processed")

                self.log(f"  Connected candidate paths: {len(next_paths)}")

                if not next_paths:
                    self.log(f"  Result: {tp} no connecting path -> ending this start point")
                    finished.extend(current)
                    break

                def path_cost(path):
                    return sum(c for g, c in path)

                def key_fn(path):
                    return (path_cost(path), '->'.join(g for g, c in path))

                next_paths.sort(key=key_fn)

                if len(next_paths) > K:
                    # ── Adaptive beam width: automatically enlarge K in proportion to the number of targets ──
                    # Base K uses the user setting and expands when there are many targets
                    # Keep cost-based truncation unchanged (for fairness)
                    unique_targets_in_paths = len(set(p[-1][0] for p in next_paths))
                    adaptive_K = max(K, unique_targets_in_paths)
                    # Upper limit: 5x K or 5000 (memory/speed safety guard)
                    max_K = min(K * 5, 5000)
                    adaptive_K = min(adaptive_K, max_K)

                    if adaptive_K > K:
                        self.log(f"  [Beam] Targets detected: {unique_targets_in_paths} -> "
                                 f"adaptive beam width: {K} → {adaptive_K}")

                    next_paths = next_paths[:adaptive_K]

                    # Report the number of targets actually covered
                    covered = len(set(p[-1][0] for p in next_paths))
                    self.log(f"  [Beam] Keeping top {len(next_paths)} "
                             f"(target coverage: {covered}/{unique_targets_in_paths} items)")
                    self.log(f"  [Cost Range] Min: {path_cost(next_paths[0]):.3f}, "
                             f"Max: {path_cost(next_paths[-1]):.3f}")

                current = next_paths

                self.log(f"  [Time Point Done] {tp}: current paths {len(current)} kept")

                if hasattr(self.api, 'print_stats_brief'):
                    self.api.print_stats_brief()
            else:
                finished.extend(current)

            for p in finished:
                path_genes = tuple(g for g, c in p)
                all_paths.add(path_genes)

        self.log(f"=== Full pathway tracing completed: total {len(all_paths)} complete paths, {len(all_intermediate_paths)} intermediate paths ===")

        if is_regulatory and self.regulatory_info:
            self.log(f"\n=== Regulatory Relationship Summary ===")
            act_count = sum(1 for v in self.regulatory_info.values()
                          if v.get("action") in ("activation", "catalysis", "expression"))
            inh_count = sum(1 for v in self.regulatory_info.values()
                          if v.get("action") in ("inhibition",))
            dir_count = sum(1 for v in self.regulatory_info.values() if v.get("is_directed"))
            self.log(f"  Collected regulatory relationships: {len(self.regulatory_info)}")
            self.log(f"  activation: {act_count}, inhibition: {inh_count}, Directional edges: {dir_count}")

        return list(all_paths), name_to_gene, list(all_intermediate_paths)


# ========== Advanced Dialog ==========

class AdvancedDialog(tk.Toplevel):
    def __init__(self, master, initial):
        super().__init__(master)
        self.title("Deep Search Profile (STRING v12.5)")
        self.resizable(False, False)
        self.result = dict(initial or {})
        self.v_degree = tk.IntVar(value=int(initial.get("degree_cap", 50)))
        self.v_beam = tk.IntVar(value=int(initial.get("beam_width", 150)))
        self.v_maxpt = tk.IntVar(value=int(initial.get("max_paths_per_target", 300)))
        self.v_score = tk.IntVar(value=int(initial.get("min_string_score", 700)))
        self.v_verbose = tk.BooleanVar(value=bool(initial.get("verbose", False)))
        self.v_collapse = tk.BooleanVar(value=bool(initial.get("collapse_gene_level", True)))
        self.v_back = tk.BooleanVar(value=bool(initial.get("no_immediate_backtrack", True)))
        self.v_dedupe = tk.BooleanVar(value=bool(initial.get("dedupe_gene_in_path", True)))
        self.v_hubs = tk.StringVar(value=",".join(initial.get("hub_blacklist", [])))
        self.v_reg_bonus = tk.DoubleVar(value=float(initial.get("regulatory_weight_bonus", 0.1)))
        self.v_temporal_pair_cap = tk.IntVar(
            value=int(initial.get("temporal_segment_route_cap", 5))
        )
        self.v_temporal_beam = tk.IntVar(
            value=int(initial.get("temporal_beam_width", 0))
        )
        self.v_temporal_per_target = tk.IntVar(
            value=int(initial.get("temporal_routes_per_target", 3))
        )
        self.v_temporal_full_dedupe = tk.BooleanVar(
            value=bool(initial.get("temporal_full_path_dedupe", True))
        )

        frm = tk.Frame(self, padx=10, pady=10)
        frm.pack(fill=tk.BOTH, expand=True)

        def row(label, w):
            r = tk.Frame(frm)
            r.pack(fill=tk.X, pady=3)
            tk.Label(r, text=label, width=35, anchor="w").pack(side=tk.LEFT)
            w.pack(side=tk.LEFT, fill=tk.X, expand=True)

        row("Degree cap H (per node):", tk.Entry(frm, textvariable=self.v_degree, width=8))
        row("Segment beam K (BFS frontier):", tk.Entry(frm, textvariable=self.v_beam, width=8))
        row("Max paths per target N:", tk.Entry(frm, textvariable=self.v_maxpt, width=8))
        row("Min STRING score (0-1000):", tk.Entry(frm, textvariable=self.v_score, width=8))
        row("Hub blacklist (genes, comma):", tk.Entry(frm, textvariable=self.v_hubs, width=40))
        row("Regulatory weight bonus:", tk.Entry(frm, textvariable=self.v_reg_bonus, width=8))
        row(
            "Temporal segment routes/pair:",
            tk.Entry(frm, textvariable=self.v_temporal_pair_cap, width=8)
        )
        row(
            "Temporal base K (0 = segment K):",
            tk.Entry(frm, textvariable=self.v_temporal_beam, width=8)
        )
        row(
            "Temporal routes per reached target:",
            tk.Entry(frm, textvariable=self.v_temporal_per_target, width=8)
        )
        tk.Label(
            frm,
            text=(
                "FINAL Temporal v2: target-preserving adaptive beam; "
                "every reached target retains >=1 route."
            ),
            fg="gray",
            anchor="w"
        ).pack(fill=tk.X, pady=(4, 2))

        opt = tk.Frame(frm)
        opt.pack(fill=tk.X, pady=5)
        tk.Checkbutton(opt, text="Verbose logging", variable=self.v_verbose).pack(anchor="w")
        tk.Checkbutton(opt, text="Collapse by gene (dedupe isoforms)", variable=self.v_collapse).pack(anchor="w")
        tk.Checkbutton(opt, text="No immediate backtracking (A→B→A ban)", variable=self.v_back).pack(anchor="w")
        tk.Checkbutton(opt, text="Disallow duplicate genes in a path", variable=self.v_dedupe).pack(anchor="w")
        tk.Checkbutton(
            opt,
            text="Disallow gene revisit across temporal intervals",
            variable=self.v_temporal_full_dedupe
        ).pack(anchor="w")

        btns = tk.Frame(frm)
        btns.pack(fill=tk.X, pady=8)
        tk.Button(btns, text="OK", width=10, command=self.on_ok).pack(side=tk.RIGHT, padx=5)
        tk.Button(btns, text="Cancel", width=10, command=self.destroy).pack(side=tk.RIGHT)

    def on_ok(self):
        try:
            self.result["degree_cap"] = max(1, int(self.v_degree.get()))
            self.result["beam_width"] = max(1, int(self.v_beam.get()))
            self.result["max_paths_per_target"] = max(1, int(self.v_maxpt.get()))
            self.result["min_string_score"] = max(0, min(1000, int(self.v_score.get())))
            hubs = [h.strip() for h in self.v_hubs.get().split(",") if h.strip()]
            self.result["hub_blacklist"] = hubs
            self.result["verbose"] = bool(self.v_verbose.get())
            self.result["collapse_gene_level"] = bool(self.v_collapse.get())
            self.result["no_immediate_backtrack"] = bool(self.v_back.get())
            self.result["dedupe_gene_in_path"] = bool(self.v_dedupe.get())
            self.result["regulatory_weight_bonus"] = float(self.v_reg_bonus.get())
            self.result["temporal_segment_route_cap"] = max(
                1, int(self.v_temporal_pair_cap.get())
            )
            self.result["temporal_beam_width"] = max(
                0, int(self.v_temporal_beam.get())
            )
            self.result["temporal_routes_per_target"] = max(
                1, int(self.v_temporal_per_target.get())
            )
            self.result["temporal_full_path_dedupe"] = bool(
                self.v_temporal_full_dedupe.get()
            )
        except Exception as e:
            messagebox.showerror("Error", f"Invalid advanced parameters: {e}")
            return
        self.destroy()


# ========== Sequence Input Dialog ==========

class SequenceInputDialog(tk.Toplevel):
    """Dialog for direct protein sequence input"""

    def __init__(self, master):
        super().__init__(master)
        self.title("Protein Sequence Input (FASTA Format)")
        self.geometry("600x500")
        self.result = None

        frm = tk.Frame(self, padx=10, pady=10)
        frm.pack(fill=tk.BOTH, expand=True)

        tk.Label(frm, text="Enter protein amino-acid sequences in FASTA format:",
                 font=("Arial", 11, "bold")).pack(anchor="w", pady=5)

        tk.Label(frm, text="Example:\n>EGFR\nMRPSGTAGAALLALLAALCPA...\n>BRAF\nMAAPLSSPLMS...",
                 justify=tk.LEFT, fg="gray", font=("Courier", 9)).pack(anchor="w", pady=3)

        self.text_area = scrolledtext.ScrolledText(frm, height=20, width=70, wrap=tk.WORD,
                                                    font=("Courier", 10))
        self.text_area.pack(fill=tk.BOTH, expand=True, pady=5)

        # Species ID option
        species_frm = tk.Frame(frm)
        species_frm.pack(fill=tk.X, pady=5)
        tk.Label(species_frm, text="Species ID (leave blank for Human/9606):").pack(side=tk.LEFT)
        self.species_var = tk.StringVar(value="")
        tk.Entry(species_frm, textvariable=self.species_var, width=20).pack(side=tk.LEFT, padx=5)
        tk.Label(species_frm, text="(STRG0AXXXXX format allowed)", fg="gray").pack(side=tk.LEFT)

        # BLAST Settings
        blast_frm = tk.LabelFrame(frm, text="BLAST Settings (used when name mapping fails)", padx=5, pady=3)
        blast_frm.pack(fill=tk.X, pady=5)

        blast_row1 = tk.Frame(blast_frm)
        blast_row1.pack(fill=tk.X, pady=2)
        tk.Label(blast_row1, text="Minimum Identity (%):").pack(side=tk.LEFT)
        self.identity_var = tk.DoubleVar(value=30.0)
        tk.Entry(blast_row1, textvariable=self.identity_var, width=8).pack(side=tk.LEFT, padx=5)
        tk.Label(blast_row1, text="E-value threshold:").pack(side=tk.LEFT, padx=(15,0))
        self.evalue_var = tk.StringVar(value="1e-10")
        tk.Entry(blast_row1, textvariable=self.evalue_var, width=10).pack(side=tk.LEFT, padx=5)

        tk.Label(blast_frm,
                 text="Strategy: 1) try STRING name mapping -> 2) if that fails, use NCBI BLAST (SwissProt)\n"
                      "to search for the most similar existing protein -> 3) use that protein's STRING PPI as a proxy",
                 fg="gray", font=("Arial", 9), justify=tk.LEFT).pack(anchor="w", pady=2)

        # Load from file button
        btn_frm = tk.Frame(frm)
        btn_frm.pack(fill=tk.X, pady=5)
        tk.Button(btn_frm, text="Load FASTA File", command=self.load_fasta).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frm, text="Set as Start Proteins", command=self.on_ok,
                  bg="lightgreen").pack(side=tk.RIGHT, padx=5)
        tk.Button(btn_frm, text="Cancel", command=self.on_cancel).pack(side=tk.RIGHT)

    def load_fasta(self):
        path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa"), ("All files", "*.*")]
        )
        if path:
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    self.text_area.delete(1.0, tk.END)
                    self.text_area.insert(tk.END, f.read())
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")

    def on_ok(self):
        text = self.text_area.get(1.0, tk.END).strip()
        if not text:
            messagebox.showwarning("Warning", "Please enter sequence data.")
            return

        sequences = parse_fasta_sequences(text)
        if not sequences:
            messagebox.showwarning("Warning", "No valid FASTA sequences were recognized.")
            return

        self.result = {
            "sequences": sequences,
            "species_id": self.species_var.get().strip(),
            "blast_identity_threshold": self.identity_var.get(),
            "blast_evalue_threshold": float(self.evalue_var.get())
        }
        self.destroy()

    def on_cancel(self):
        self.result = None
        self.destroy()




# ========== Automatic Cell-Context Data Manager ==========

CMP_API_BASE = "https://api.cellmodelpassports.sanger.ac.uk"
CMP_WEB_BASE = "https://cellmodelpassports.sanger.ac.uk"
CMP_MODEL_CACHE_TTL_DAYS = 30

# High-confidence cross-references used only to identify the SAME cell model.
# Never substitute a similarly named but different cell line.
# MDA-MB-468: Cell Model Passports SIDM00628; Cellosaurus CVCL_0419;
# DepMap ACH-000849.
CMP_KNOWN_XREFS = {
    "SIDM00628": {
        "name": "MDA-MB-468",
        "aliases": [
            "MDA-MB-468", "MDA-MB 468", "MDA-MB468", "MDAMB468",
            "MDA-468", "MDA468", "MB468",
            "SIDM00628", "CVCL_0419", "ACH-000849", "908123"
        ],
        "cellosaurus": "CVCL_0419",
        "depmap": "ACH-000849",
        "cosmic": "908123",
    },
    "__HEK293T__": {
        "name": "HEK293T",
        "aliases": [
            "HEK293T", "HEK-293T", "HEK 293T", "293T", "293-T", "293 T",
            "Human Embryonic Kidney 293T", "293tsA1609neo", "CVCL_0063"
        ],
        "cellosaurus": "CVCL_0063",
        "depmap": "",
        "cosmic": "",
        "requires_cmp_exact_match": True,
    },
}

CMP_KNOWN_QUERY_TO_SIDM = {}
for _sidm, _xref in CMP_KNOWN_XREFS.items():
    _tokens = [_xref.get("name", "")] + list(_xref.get("aliases") or [])
    if not str(_sidm).startswith("__"):
        _tokens = [_sidm] + _tokens
    for _token in _tokens:
        _k = _norm_cell_name(_token)
        if _k:
            CMP_KNOWN_QUERY_TO_SIDM[_k] = _sidm


def _cell_context_cache_dir() -> str:
    try:
        base = os.path.dirname(os.path.abspath(__file__))
    except Exception:
        base = os.getcwd()
    candidate = os.path.join(base, "cell_context_data")
    try:
        os.makedirs(candidate, exist_ok=True)
        test = os.path.join(candidate, ".__write_test__")
        with open(test, "w", encoding="utf-8") as f:
            f.write("ok")
        os.remove(test)
        return candidate
    except Exception:
        fallback = os.path.join(os.path.expanduser("~"), ".ts_bfs_bs_cell_context")
        os.makedirs(fallback, exist_ok=True)
        return fallback


def _cmp_model_cache_path():
    return os.path.join(_cell_context_cache_dir(), "cmp_models.json")


def _cmp_model_dir(sidm: str):
    p = os.path.join(_cell_context_cache_dir(), str(sidm).strip())
    os.makedirs(p, exist_ok=True)
    return p


def _is_usable_auto_file(path: str, min_bytes: int = 256) -> bool:
    try:
        return os.path.isfile(path) and os.path.getsize(path) > min_bytes
    except Exception:
        return False


def _cmp_session():
    s = requests.Session()
    retry = Retry(
        total=3, connect=3, read=3,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.headers.update({
        "User-Agent": "TS-BFS-BS-CellContext/2.0 (local academic research client)",
        "Accept": "application/vnd.api+json, application/json;q=0.9, */*;q=0.5",
    })
    return s


def _cmp_get_json(path: str, params=None, timeout=60):
    url = path if str(path).startswith("http") else CMP_API_BASE + str(path)
    s = _cmp_session()
    r = s.get(url, params=params or {}, timeout=timeout)
    r.raise_for_status()
    ctype = str(r.headers.get("content-type") or "").lower()
    if "html" in ctype:
        raise RuntimeError("Cell Model Passports API returned HTML instead of JSON.")
    try:
        doc = r.json()
    except Exception as e:
        raise RuntimeError(f"Cell Model Passports API returned invalid JSON: {e}")
    if isinstance(doc, dict) and doc.get("errors"):
        raise RuntimeError(f"Cell Model Passports API error: {doc['errors']}")
    return doc


def _cmp_name_strings(value):
    out = []
    if value is None:
        return out
    if isinstance(value, str):
        s = value.strip()
        if s:
            out.append(s)
        return out
    if isinstance(value, (list, tuple)):
        for x in value:
            out.extend(_cmp_name_strings(x))
        return out
    if isinstance(value, dict):
        for k in ("name", "model_name", "cell_line_name", "value", "synonym"):
            if k in value:
                out.extend(_cmp_name_strings(value.get(k)))
        if not out:
            for v in value.values():
                if isinstance(v, (str, list, tuple, dict)):
                    out.extend(_cmp_name_strings(v))
        return out
    return out


def _cmp_enrich_model_record(model):
    """Attach trusted aliases/identifiers without merging distinct cell lines."""
    m = dict(model or {})
    sidm = str(m.get("id") or "").strip()
    aliases = [str(x).strip() for x in (m.get("aliases") or []) if str(x).strip()]
    identifiers = [str(x).strip() for x in (m.get("identifiers") or []) if str(x).strip()]

    xref = CMP_KNOWN_XREFS.get(sidm)
    if xref:
        # Prefer the verified canonical name for a known SIDM.
        m["name"] = xref["name"]
        aliases = [xref["name"]] + list(xref.get("aliases") or []) + aliases
        identifiers += [
            sidm,
            xref.get("cellosaurus", ""),
            xref.get("depmap", ""),
            xref.get("cosmic", ""),
        ]

    seen = set()
    clean_aliases = []
    for x in aliases:
        x = str(x or "").strip()
        k = _norm_cell_name(x)
        if x and k and k not in seen:
            seen.add(k)
            clean_aliases.append(x)

    seen_ids = set()
    clean_ids = []
    for x in identifiers:
        x = str(x or "").strip()
        k = _norm_cell_name(x)
        if x and k and k not in seen_ids:
            seen_ids.add(k)
            clean_ids.append(x)

    primary = str(m.get("name") or (clean_aliases[0] if clean_aliases else sidm)).strip()
    m.update({
        "id": sidm,
        "name": primary,
        "aliases": clean_aliases,
        "identifiers": clean_ids,
        "display": f"{primary}  [{sidm}]" if sidm else primary,
    })
    return m


def _cmp_model_record(row):
    sidm = str((row or {}).get("id") or "").strip()
    attrs = (row or {}).get("attributes") or {}
    names = _cmp_name_strings(attrs.get("names"))
    if not names:
        for key in ("name", "model_name", "cell_line_name"):
            names.extend(_cmp_name_strings(attrs.get(key)))

    # Capture identifier-like strings if the API includes them. The current API
    # may return only names in sparse-field mode, so trusted cross-references are
    # also injected by _cmp_enrich_model_record().
    identifiers = []
    for key in ("identifiers", "identifier", "external_ids", "external_identifiers",
                "cellosaurus", "depmap_id", "broad_id", "cosmic_id", "ccle_id"):
        if key in attrs:
            identifiers.extend(_cmp_name_strings(attrs.get(key)))

    cleaned, seen = [], set()
    for n in names:
        n = str(n).strip()
        if not n or n.startswith(("http://", "https://")):
            continue
        key = _norm_cell_name(n)
        if not key or key in seen:
            continue
        seen.add(key)
        cleaned.append(n)

    primary = cleaned[0] if cleaned else sidm
    return _cmp_enrich_model_record({
        "id": sidm,
        "name": primary,
        "aliases": cleaned,
        "identifiers": identifiers,
    })


def _cmp_ensure_known_models(models):
    """Ensure verified SIDMs remain selectable even if the display-name feed omits them."""
    by_id = {str(m.get("id") or "").strip(): _cmp_enrich_model_record(m) for m in (models or [])}
    for sidm, xref in CMP_KNOWN_XREFS.items():
        if str(sidm).startswith("__") or xref.get("requires_cmp_exact_match"):
            continue
        if sidm not in by_id:
            by_id[sidm] = _cmp_enrich_model_record({
                "id": sidm,
                "name": xref["name"],
                "aliases": list(xref.get("aliases") or []),
                "identifiers": [xref.get("cellosaurus", ""), xref.get("depmap", ""), xref.get("cosmic", "")],
            })
    out = list(by_id.values())
    out.sort(key=lambda x: (_norm_cell_name(x.get("name")), x.get("id", "")))
    return out


def _cmp_fetch_models(progress_cb=None):
    models = []
    page = 1
    page_size = 100
    total = None
    while True:
        params = {
            "fields[model]": "names",
            "page[size]": str(page_size),
            "page[number]": str(page),
        }
        doc = _cmp_get_json("/models", params=params, timeout=90)
        rows = doc.get("data") or []
        if total is None:
            try:
                total = int((doc.get("meta") or {}).get("count") or 0)
            except Exception:
                total = 0
        for row in rows:
            rec = _cmp_model_record(row)
            if rec["id"]:
                models.append(rec)
        if progress_cb:
            progress_cb(len(models), total or 0)
        if not rows:
            break
        if total and len(models) >= total:
            break
        if len(rows) < page_size and not total:
            break
        page += 1
        time.sleep(0.55)
    models = _cmp_ensure_known_models(models)
    cache = {
        "source": CMP_API_BASE + "/models",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "models": models,
    }
    Path(_cmp_model_cache_path()).write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")
    return models


def _cmp_load_cached_models(allow_stale=True):
    p = _cmp_model_cache_path()
    if not os.path.isfile(p):
        return []
    try:
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        models = d.get("models") or []
        if not isinstance(models, list):
            return []
        if not allow_stale:
            age = time.time() - os.path.getmtime(p)
            if age > CMP_MODEL_CACHE_TTL_DAYS * 86400:
                return []
        return _cmp_ensure_known_models(models)
    except Exception:
        return []


def _flatten_strings(obj, key_path=""):
    if isinstance(obj, str):
        yield key_path, obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{key_path}.{k}" if key_path else str(k)
            yield from _flatten_strings(v, p)
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            yield from _flatten_strings(v, f"{key_path}[{i}]")


def _best_download_url_from_file_row(row, wanted_token):
    strings = list(_flatten_strings(row))
    token = wanted_token.lower()
    filename = ""
    for _, s in strings:
        ss = str(s).strip()
        if token in ss.lower() and (".csv" in ss.lower() or ".zip" in ss.lower()):
            filename = os.path.basename(urlparse(ss).path) or os.path.basename(ss)
            break
    urls = []
    for key, s in strings:
        ss = html.unescape(str(s).strip())
        if ss.startswith(("http://", "https://", "/")):
            full = urljoin(CMP_API_BASE + "/", ss)
            lo, kp = full.lower(), key.lower()
            score = 0
            if token in lo: score += 100
            if "download" in lo or "download" in kp: score += 40
            if lo.endswith((".zip", ".csv", ".gz", ".tsv")): score += 30
            if "cellmodelpassports.sanger.ac.uk" in lo and "api." not in lo: score += 10
            if "/models/" in lo and "/files" in lo: score -= 20
            urls.append((score, full))
    urls.sort(reverse=True)
    return filename, (urls[0][1] if urls and urls[0][0] > 0 else "")


class _PassportLinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self._href = None
        self._text = []
    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            self._href = dict(attrs).get("href")
            self._text = []
    def handle_data(self, data):
        if self._href is not None:
            self._text.append(data)
    def handle_endtag(self, tag):
        if tag.lower() == "a" and self._href is not None:
            self.links.append((self._href, " ".join(self._text).strip()))
            self._href = None
            self._text = []


def _cmp_scrape_passport_file(sidm: str, wanted_token: str):
    page = f"{CMP_WEB_BASE}/passports/{sidm}"
    s = _cmp_session()
    r = s.get(page, timeout=60, headers={"Accept": "text/html,*/*"})
    r.raise_for_status()
    parser = _PassportLinkParser()
    parser.feed(r.text)
    token = wanted_token.lower()
    candidates = []
    for href, label in parser.links:
        combined = f"{href} {label}".lower()
        if token not in combined:
            continue
        full = urljoin(page, html.unescape(href))
        score = 100 if token in full.lower() else 0
        if full.lower().endswith((".zip", ".csv", ".gz")): score += 50
        if "processed" in combined: score += 20
        candidates.append((score, full))
    if candidates:
        candidates.sort(reverse=True)
        return candidates[0][1]
    pat = re.compile(r'''href=["']([^"']*''' + re.escape(wanted_token) + r'''[^"']*)["']''', re.I)
    m = pat.search(r.text)
    return urljoin(page, html.unescape(m.group(1))) if m else ""


def _cmp_model_files(sidm: str):
    doc = _cmp_get_json(f"/models/{sidm}/files", params={"page[size]": "100", "page[number]": "1"}, timeout=90)
    rows = doc.get("data") or []
    found = {
        "rna": {"filename": f"{sidm}_rnaseq.csv.zip", "url": ""},
        "proteomics": {"filename": f"{sidm}_proteomics.csv.zip", "url": ""},
    }
    for row in rows:
        blob = " ".join(s for _, s in _flatten_strings(row)).lower()
        if "rnaseq" in blob or "rna_seq" in blob or "rna-seq" in blob:
            fn, u = _best_download_url_from_file_row(row, "rnaseq")
            if fn: found["rna"]["filename"] = fn
            if u: found["rna"]["url"] = u
        if "proteomics" in blob:
            fn, u = _best_download_url_from_file_row(row, "proteomics")
            if fn: found["proteomics"]["filename"] = fn
            if u: found["proteomics"]["url"] = u
    if not found["rna"]["url"]:
        found["rna"]["url"] = _cmp_scrape_passport_file(sidm, "rnaseq")
    if not found["proteomics"]["url"]:
        found["proteomics"]["url"] = _cmp_scrape_passport_file(sidm, "proteomics")
    return found


def _download_stream(url, dest: str, progress_cb=None, timeout=120):
    if not url:
        raise RuntimeError("No download URL was found.")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    try:
        if os.path.exists(tmp): os.remove(tmp)
    except Exception:
        pass
    s = _cmp_session()
    with s.get(url, stream=True, timeout=timeout, allow_redirects=True, headers={"Accept": "application/zip,text/csv,application/octet-stream,*/*"}) as r:
        r.raise_for_status()
        ctype = str(r.headers.get("content-type") or "").lower()
        total = int(r.headers.get("content-length") or 0)
        done, first = 0, b""
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if not chunk: continue
                if not first: first = chunk[:512]
                f.write(chunk)
                done += len(chunk)
                if progress_cb: progress_cb(done, total)
    size = os.path.getsize(tmp) if os.path.isfile(tmp) else 0
    low = first.lstrip().lower()
    if "html" in ctype or low.startswith(b"<!doctype html") or low.startswith(b"<html"):
        try: os.remove(tmp)
        except Exception: pass
        raise RuntimeError("The server returned an HTML page instead of the processed data file.")
    if size < 128:
        try: os.remove(tmp)
        except Exception: pass
        raise RuntimeError(f"Downloaded file is unexpectedly small ({size} bytes).")
    os.replace(tmp, dest)
    return dest


def _cmp_write_local_metadata(model):
    sidm = str(model.get("id") or "").strip()
    name = str(model.get("name") or sidm).strip()
    aliases = [str(x).strip() for x in (model.get("aliases") or []) if str(x).strip()]
    folder = _cmp_model_dir(sidm)
    p = os.path.join(folder, f"{sidm}_metadata.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["CellLineName", "ModelID", "Alias1", "Alias2", "Alias3"])
        w.writeheader()
        w.writerow({
            "CellLineName": name,
            "ModelID": sidm,
            "Alias1": aliases[1] if len(aliases) > 1 else "",
            "Alias2": aliases[2] if len(aliases) > 2 else "",
            "Alias3": aliases[3] if len(aliases) > 3 else "",
        })
    return p


def _cmp_cached_model_paths(model):
    sidm = str(model.get("id") or "").strip()
    folder = _cmp_model_dir(sidm)
    rna = os.path.join(folder, f"{sidm}_rnaseq.csv.zip")
    prot = os.path.join(folder, f"{sidm}_proteomics.csv.zip")
    meta = os.path.join(folder, f"{sidm}_metadata.csv")
    return {
        "rna": rna if _is_usable_auto_file(rna) else "",
        "proteomics": prot if _is_usable_auto_file(prot) else "",
        "metadata": meta if _is_usable_auto_file(meta, min_bytes=20) else "",
    }


# ========== Cell Context Dialog ==========

class CellContextDialog(tk.Toplevel):
    def __init__(self, master, initial=None):
        super().__init__(master)
        initial = dict(initial or {})
        self.title("Cell-Context Bridge Filter — Cell Model Passports")
        self.geometry("980x640")
        self.minsize(920, 600)
        self.resizable(True, True)
        self.result = None
        self.v_enabled = tk.BooleanVar(value=bool(initial.get("enabled", False)))
        self.v_cell = tk.StringVar(value=str(initial.get("cell_line", "")))
        self.v_rna = tk.StringVar(value=str(initial.get("rna_file", "")))
        self.v_prot = tk.StringVar(value=str(initial.get("proteome_file", "")))
        self.v_meta = tk.StringVar(value=str(initial.get("metadata_file", "")))
        self.v_mode = tk.StringVar(value=str(initial.get("mode", "conservative_union")))
        self.v_thr = tk.DoubleVar(value=float(initial.get("rna_threshold", 1.0) or 1.0))
        self.v_status = tk.StringVar(value="Checking Cell Model Passports model cache...")
        self.v_model_search = tk.StringVar(value="")
        self._models = []
        self._display_to_model = {}
        self._selected_model = None
        self._busy = False

        frm = tk.Frame(self, padx=12, pady=10)
        frm.pack(fill=tk.BOTH, expand=True)
        tk.Checkbutton(frm, text="Enable cell-context bridge filtering", variable=self.v_enabled, font=("Arial", 10, "bold")).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 8))

        auto = tk.LabelFrame(frm, text="Cell Model Passports automatic setup", padx=8, pady=7)
        auto.grid(row=1, column=0, columnspan=5, sticky="ew", pady=(0, 10))
        tk.Label(auto, text="Uses the official Sanger Cell Model Passports API. Model-specific processed files are cached locally.", fg="#333").grid(row=0, column=0, columnspan=4, sticky="w")
        tk.Button(auto, text="Refresh Model List", width=18, command=lambda: self.refresh_model_list(force=True), bg="lightcyan").grid(row=1, column=0, sticky="w", pady=(7, 2))
        tk.Button(auto, text="Load Selected Model", width=19, command=self.download_selected_model).grid(row=1, column=1, sticky="w", padx=5, pady=(7, 2))
        tk.Button(auto, text="Use Existing Cache", width=18, command=self.use_selected_cache).grid(row=1, column=2, sticky="w", padx=5, pady=(7, 2))
        tk.Button(auto, text="Open Cache Folder", width=18, command=self.open_cache_folder).grid(row=1, column=3, sticky="w", padx=5, pady=(7, 2))

        tk.Label(auto, text="Find model:").grid(row=2, column=0, sticky="w", pady=(6, 2))
        search_entry = tk.Entry(auto, textvariable=self.v_model_search, width=38)
        search_entry.grid(row=2, column=1, columnspan=2, sticky="ew", padx=5, pady=(6, 2))
        search_entry.bind("<Return>", lambda e: self.find_model())
        tk.Button(auto, text="Find by name / ID", width=18, command=self.find_model).grid(row=2, column=3, sticky="w", padx=5, pady=(6, 2))
        tk.Label(
            auto,
            text="Accepts cell-line name, SIDM, Cellosaurus (CVCL_), DepMap (ACH-), or known alias.",
            fg="gray", anchor="w"
        ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(0, 2))

        tk.Label(auto, textvariable=self.v_status, fg="#245", anchor="w").grid(row=4, column=0, columnspan=4, sticky="ew", pady=(4, 0))
        auto.columnconfigure(1, weight=1)
        auto.columnconfigure(2, weight=1)
        auto.columnconfigure(3, weight=1)

        tk.Label(frm, text="Cell line / model:").grid(row=2, column=0, sticky="w", pady=5)
        self.cell_combo = ttk.Combobox(frm, textvariable=self.v_cell, state="readonly", width=52, values=[])
        self.cell_combo.grid(row=2, column=1, columnspan=3, sticky="ew", pady=5)
        self.cell_combo.bind("<<ComboboxSelected>>", self._on_model_selected)
        tk.Label(frm, text="Name + SIDM shown together to avoid ambiguity.", fg="gray").grid(row=2, column=4, sticky="w", padx=6)

        tk.Label(frm, text="Mode:").grid(row=3, column=0, sticky="w", pady=4)
        ttk.Combobox(frm, textvariable=self.v_mode, state="readonly", width=24, values=["conservative_union", "rna_only", "proteome_only", "strict_intersection"]).grid(row=3, column=1, sticky="w")
        tk.Label(frm, text="RNA threshold:").grid(row=3, column=2, sticky="e", padx=(20, 4))
        tk.Entry(frm, textvariable=self.v_thr, width=10).grid(row=3, column=3, sticky="w")

        def filerow(r, label, var):
            tk.Label(frm, text=label).grid(row=r, column=0, sticky="w", pady=4)
            tk.Entry(frm, textvariable=var, width=78).grid(row=r, column=1, columnspan=3, sticky="ew", pady=4)
            def browse():
                p = filedialog.askopenfilename(title=label, filetypes=[("Data files", "*.csv *.tsv *.txt *.gz *.zip"), ("All files", "*.*")])
                if p: var.set(p)
            tk.Button(frm, text="Browse...", command=browse).grid(row=r, column=4, sticky="w", padx=6)

        filerow(4, "RNA file:", self.v_rna)
        filerow(5, "Proteomics file:", self.v_prot)
        filerow(6, "Cell-line metadata:", self.v_meta)
        tk.Label(frm, text="Cache folder:").grid(row=7, column=0, sticky="nw", pady=(8, 2))
        tk.Label(frm, text=_cell_context_cache_dir(), fg="#555", justify=tk.LEFT, wraplength=760, anchor="w").grid(row=7, column=1, columnspan=4, sticky="w", pady=(8, 2))

        notes = (
            "Recommended: Conservative Union. Experimental phosphoproteomic anchors remain protected.\n"
            "Only STRING-only bridge candidates are checked. Selecting a model reuses cached SIDM files if present; otherwise it downloads that model's processed RNA/proteomics.\n"
            "Proteomics is optional. If no processed proteomics file exists, it is recorded as unavailable/UNKNOWN rather than protein absence.\n"
            "Manual Browse remains available for project-specific local files.\n"
            "Exact identity only: similarly named lines are never substituted. MDA-MB-TE is not used for MDA-MB-468, and parental HEK293/other 293 derivatives are not used for HEK293T."
        )
        tk.Label(frm, text=notes, justify=tk.LEFT, fg="#444", wraplength=920).grid(row=8, column=0, columnspan=5, sticky="w", pady=(12, 8))
        btn = tk.Frame(frm); btn.grid(row=9, column=0, columnspan=5, sticky="e", pady=8)
        tk.Button(btn, text="OK", width=12, command=self.on_ok).pack(side=tk.RIGHT, padx=4)
        tk.Button(btn, text="Cancel", width=12, command=self.destroy).pack(side=tk.RIGHT, padx=4)
        frm.columnconfigure(1, weight=1); frm.columnconfigure(2, weight=1); frm.columnconfigure(3, weight=1); frm.rowconfigure(8, weight=1)
        self.after(80, self._startup)
        self.transient(master); self.grab_set()

    def _set_models(self, models, preserve_name=True):
        current = self.v_cell.get().strip() if preserve_name else ""
        self._models = _cmp_ensure_known_models(models)
        self._display_to_model = {m["display"]: m for m in self._models}
        self.cell_combo["values"] = [m["display"] for m in self._models]
        self._selected_model = None
        if current:
            q = _norm_cell_name(current)
            for m in self._models:
                candidates = (
                    [m.get("display", ""), m.get("name", ""), m.get("id", "")]
                    + list(m.get("aliases") or [])
                    + list(m.get("identifiers") or [])
                )
                if any(_norm_cell_name(x) == q for x in candidates):
                    self.v_cell.set(m["display"]); self._selected_model = m; return

    def find_model(self):
        """Exact cross-reference search; never substitutes a merely similar name."""
        raw = self.v_model_search.get().strip()
        if not raw:
            messagebox.showinfo(
                "Find Cell Model",
                "Enter a cell-line name, SIDM, Cellosaurus ID (CVCL_), or DepMap ID (ACH-)."
            )
            return

        q = _norm_cell_name(raw)
        known_key = CMP_KNOWN_QUERY_TO_SIDM.get(q)
        matches = []

        for m in self._models:
            candidates = (
                [m.get("display", ""), m.get("name", ""), m.get("id", "")]
                + list(m.get("aliases") or [])
                + list(m.get("identifiers") or [])
            )
            normalized = {_norm_cell_name(x) for x in candidates if str(x or "").strip()}
            if q in normalized:
                matches.append(m)
                continue
            if known_key:
                xref = CMP_KNOWN_XREFS.get(known_key) or {}
                trusted = {
                    _norm_cell_name(x)
                    for x in (
                        [xref.get("name", ""), xref.get("cellosaurus", ""),
                         xref.get("depmap", ""), xref.get("cosmic", "")]
                        + list(xref.get("aliases") or [])
                    )
                    if str(x or "").strip()
                }
                if normalized.intersection(trusted):
                    matches.append(m)

        uniq = {}
        for m in matches:
            uniq[(m.get("id", ""), m.get("display", ""))] = m
        matches = list(uniq.values())

        if len(matches) == 1:
            m = matches[0]
            self.v_cell.set(m["display"])
            self._selected_model = m
            ids = ", ".join(x for x in m.get("identifiers", []) if x)
            self.v_status.set(
                f"Exact model match: {m['name']} [{m['id']}]"
                + (f" | {ids}" if ids else "")
            )
            self._on_model_selected()
            return

        if len(matches) > 1:
            messagebox.showwarning(
                "Find Cell Model",
                "More than one exact identifier/alias match was found. "
                "Please choose the intended SIDM from the list."
            )
            return

        if known_key == "__HEK293T__":
            messagebox.showwarning(
                "HEK293T exact CMP model not found",
                "No exact Cell Model Passports record matched HEK293T / HEK-293T / "
                "293T / CVCL_0063.\n\n"
                "The program will NOT substitute HEK293 or another 293-derived "
                "cell line. For this run, either leave cell-context filtering OFF "
                "(UNKNOWN retained under Conservative Union) or provide validated "
                "local HEK293T RNA/proteomics files."
            )
            self.v_status.set(
                "HEK293T not found exactly in CMP; no HEK293 substitution was made."
            )
            return

        if known_key and not str(known_key).startswith("__"):
            xref = CMP_KNOWN_XREFS[known_key]
            model = _cmp_enrich_model_record({
                "id": known_key,
                "name": xref["name"],
                "aliases": xref.get("aliases", []),
                "identifiers": [xref.get("cellosaurus", ""), xref.get("depmap", ""), xref.get("cosmic", "")],
            })
            self._models = _cmp_ensure_known_models(self._models + [model])
            self._display_to_model = {m["display"]: m for m in self._models}
            self.cell_combo["values"] = [m["display"] for m in self._models]
            self.v_cell.set(model["display"])
            self._selected_model = model
            self.v_status.set(f"Verified model identity: {model['name']} [{model['id']}]")
            self._on_model_selected()
            return

        messagebox.showwarning(
            "Cell Model Not Found",
            f"No exact Cell Model Passports match was found for '{raw}'.\n\n"
            "A similarly named model will NOT be substituted automatically. "
            "For Conservative Union, leave cell-context filtering off or supply "
            "validated local RNA/proteomics for the correct cell line."
        )

    def _startup(self):
        models = _cmp_load_cached_models(allow_stale=True)
        if models:
            self._set_models(models)
            self.v_status.set(f"Cached model list loaded: {len(models):,} models.")
            if self._selected_model: self.use_selected_cache(silent=True)
        else:
            self.v_status.set("Retrieving model list from Cell Model Passports...")
            self.refresh_model_list(force=True)

    def refresh_model_list(self, force=False):
        if self._busy: return
        if not force:
            models = _cmp_load_cached_models(allow_stale=False)
            if models:
                self._set_models(models); self.v_status.set(f"Model list loaded from cache: {len(models):,}."); return
        self._busy = True
        self.v_status.set("Retrieving Cell Model Passports models...")
        threading.Thread(target=self._refresh_worker, daemon=True).start()

    def _refresh_worker(self):
        try:
            def progress(done, total):
                msg = f"Retrieving model list: {done:,}/{total:,}" if total else f"Retrieving model list: {done:,}"
                self.after(0, lambda m=msg: self.v_status.set(m))
            models = _cmp_fetch_models(progress_cb=progress)
            self.after(0, lambda m=models: self._finish_refresh(m))
        except Exception as e:
            cached = _cmp_load_cached_models(allow_stale=True)
            self.after(0, lambda: self._finish_refresh_error(str(e), cached))

    def _finish_refresh(self, models):
        self._busy = False; self._set_models(models)
        self.v_status.set(f"Cell Model Passports model list ready: {len(models):,} models.")

    def _finish_refresh_error(self, err, cached):
        self._busy = False
        if cached:
            self._set_models(cached)
            self.v_status.set(f"Live refresh failed; using cached list ({len(cached):,}).")
            messagebox.showwarning("Cell Model Passports", f"Could not refresh live model list:\n{err}\n\nCached list loaded instead.")
        else:
            self.v_status.set("Could not retrieve Cell Model Passports model list.")
            messagebox.showerror("Cell Model Passports", f"Could not retrieve the model list:\n{err}\n\nManual local-file selection remains available.")

    def _on_model_selected(self, event=None):
        self._selected_model = self._display_to_model.get(self.v_cell.get().strip())
        if not self._selected_model: return
        self.v_rna.set(""); self.v_prot.set(""); self.v_meta.set("")
        if self.use_selected_cache(silent=True):
            self.v_enabled.set(True); return
        self.download_selected_model()

    def use_selected_cache(self, silent=False):
        model = self._selected_model or self._display_to_model.get(self.v_cell.get().strip())
        if not model:
            if not silent: messagebox.showinfo("Cell Context", "Choose a cell line/model first.")
            return False
        p = _cmp_cached_model_paths(model)
        if not p["metadata"]: p["metadata"] = _cmp_write_local_metadata(model)
        if not p["rna"]:
            if not silent: messagebox.showinfo("Cell Context Cache", f"No cached RNA file for {model['name']} ({model['id']}).")
            return False
        self.v_meta.set(p["metadata"]); self.v_rna.set(p["rna"]); self.v_prot.set(p["proteomics"])
        self.v_status.set(f"Using cache: {model['name']} [{model['id']}] | RNA=yes | Proteomics={'yes' if p['proteomics'] else 'not available/cached'}")
        return True

    def download_selected_model(self):
        if self._busy: return
        model = self._selected_model or self._display_to_model.get(self.v_cell.get().strip())
        if not model:
            messagebox.showinfo("Cell Context", "Choose a cell line/model first."); return
        if self.use_selected_cache(silent=True):
            self.v_enabled.set(True); return
        self._busy = True
        self.v_status.set(f"Resolving attached files for {model['name']} [{model['id']}]...")
        threading.Thread(target=self._download_model_worker, args=(model,), daemon=True).start()

    def _download_model_worker(self, model):
        try:
            sidm = model["id"]; folder = _cmp_model_dir(sidm); files = _cmp_model_files(sidm)
            rna_url = files["rna"].get("url") or ""; prot_url = files["proteomics"].get("url") or ""
            if not rna_url: raise RuntimeError(f"No processed RNA-seq file is attached to {model['name']} ({sidm}).")
            rna_dest = os.path.join(folder, f"{sidm}_rnaseq.csv.zip")
            prot_dest = os.path.join(folder, f"{sidm}_proteomics.csv.zip")
            if not _is_usable_auto_file(rna_dest):
                def rp(done, total):
                    msg = f"Downloading RNA-seq: {100.0*done/total:5.1f}%" if total else f"Downloading RNA-seq: {done/1024/1024:.1f} MB"
                    self.after(0, lambda m=msg: self.v_status.set(m))
                _download_stream(rna_url, rna_dest, progress_cb=rp, timeout=180)
            proteomics_available = False
            if prot_url:
                if not _is_usable_auto_file(prot_dest):
                    def pp(done, total):
                        msg = f"Downloading proteomics: {100.0*done/total:5.1f}%" if total else f"Downloading proteomics: {done/1024/1024:.1f} MB"
                        self.after(0, lambda m=msg: self.v_status.set(m))
                    _download_stream(prot_url, prot_dest, progress_cb=pp, timeout=180)
                proteomics_available = _is_usable_auto_file(prot_dest)
            result = {"rna": rna_dest, "proteomics": prot_dest if proteomics_available else "", "metadata": _cmp_write_local_metadata(model), "proteomics_available": proteomics_available}
            self.after(0, lambda r=result, m=model: self._finish_model_download(m, r))
        except Exception as e:
            self.after(0, lambda err=str(e): self._finish_model_download_error(err))

    def _finish_model_download(self, model, result):
        self._busy = False
        self.v_rna.set(result["rna"]); self.v_prot.set(result["proteomics"]); self.v_meta.set(result["metadata"]); self.v_enabled.set(True)
        if result["proteomics_available"]:
            self.v_status.set(f"Ready: {model['name']} [{model['id']}] | RNA + proteomics loaded.")
        else:
            self.v_status.set(f"Ready: {model['name']} [{model['id']}] | RNA loaded; proteomics unavailable.")
            messagebox.showinfo("Cell Context", f"{model['name']} ({model['id']}) has RNA-seq but no attached processed proteomics file was found.\n\nProteomics will be treated as unavailable/UNKNOWN, not as protein absence.")

    def _finish_model_download_error(self, err):
        self._busy = False; self.v_status.set("Selected-model data download failed.")
        messagebox.showerror("Cell Context Download", f"Could not prepare the selected model:\n{err}\n\nYou may use Browse... to supply local files manually.")

    def open_cache_folder(self):
        folder = _cell_context_cache_dir()
        try:
            if os.name == "nt": os.startfile(folder)
            elif sys.platform == "darwin": subprocess.Popen(["open", folder])
            else: subprocess.Popen(["xdg-open", folder])
        except Exception as e:
            messagebox.showerror("Open Cache Folder", str(e))

    def on_ok(self):
        model = self._selected_model or self._display_to_model.get(self.v_cell.get().strip())
        if self.v_enabled.get() and model:
            cell_name = model.get("name") or model.get("id"); sidm = model.get("id")
        else:
            cell_name = self.v_cell.get().strip(); sidm = ""
        if self.v_enabled.get() and not cell_name:
            messagebox.showerror("Cell Context", "Please choose a cell line/model."); return
        if self.v_enabled.get() and not self.v_rna.get().strip():
            messagebox.showerror("Cell Context", "No RNA file is ready. Use Load Selected Model or choose a local RNA file."); return
        for label, p in [("RNA", self.v_rna.get().strip()), ("Proteomics", self.v_prot.get().strip()), ("Metadata", self.v_meta.get().strip())]:
            if p and not os.path.isfile(p):
                messagebox.showerror("Cell Context", f"{label} file not found:\n{p}"); return
        self.result = {
            "enabled": bool(self.v_enabled.get()),
            "cell_line": cell_name,
            "cmp_model_id": sidm,
            "rna_file": self.v_rna.get().strip(),
            "proteome_file": self.v_prot.get().strip(),
            "metadata_file": self.v_meta.get().strip(),
            "mode": self.v_mode.get().strip(),
            "rna_threshold": float(self.v_thr.get()),
            "auto_cache_dir": _cell_context_cache_dir(),
            "auto_source": "Cell Model Passports API" if sidm else "manual",
        }
        self.destroy()


# ========== PathwayGUI (v12.5 upgraded) ==========

class PathwayGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Protein Pathway Simulator (STRING API v12.5 + Cell Model Passports v1.4)")
        self.root.geometry("950x800")

        # --- Top settings ---
        settings_frame = tk.LabelFrame(root, text="Settings", padx=10, pady=5)
        settings_frame.pack(fill=tk.X, padx=10, pady=5)

        # First row: Depth + Network Type
        row1 = tk.Frame(settings_frame)
        row1.pack(fill=tk.X, pady=3)

        tk.Label(row1, text="Bridge Depth D (0-10):").pack(side=tk.LEFT, padx=5)
        self.depth_var = tk.IntVar(value=1)
        ttk.Spinbox(row1, from_=0, to=10, textvariable=self.depth_var, width=6).pack(side=tk.LEFT)

        tk.Label(row1, text="    Network Type:").pack(side=tk.LEFT, padx=5)
        self.network_type_var = tk.StringVar(value="functional")
        net_combo = ttk.Combobox(row1, textvariable=self.network_type_var, width=18,
                                  values=["functional", "physical", "regulatory"], state="readonly")
        net_combo.pack(side=tk.LEFT, padx=5)
        net_combo.bind("<<ComboboxSelected>>", self._on_network_type_change)

        # Regulatory network description label
        self.reg_label = tk.Label(row1, text="", fg="blue", font=("Arial", 9))
        self.reg_label.pack(side=tk.LEFT, padx=10)

        # Second row: Species ID
        row2 = tk.Frame(settings_frame)
        row2.pack(fill=tk.X, pady=3)

        tk.Label(row2, text="Species ID:").pack(side=tk.LEFT, padx=5)
        self.species_var = tk.StringVar(value="9606")
        tk.Entry(row2, textvariable=self.species_var, width=15).pack(side=tk.LEFT, padx=5)
        tk.Label(row2, text="(Human=9606, or STRG0AXXXXX)", fg="gray").pack(side=tk.LEFT)

        tk.Label(row2, text="    Temporal Target Mode:").pack(side=tk.LEFT, padx=(12, 3))
        self.temporal_anchor_mode_var = tk.StringVar(value="measured")
        ttk.Combobox(row2, textvariable=self.temporal_anchor_mode_var, width=14,
                     values=["measured", "strict_signed"], state="readonly").pack(side=tk.LEFT, padx=3)

        # Third row: Input File
        row3 = tk.Frame(settings_frame)
        row3.pack(fill=tk.X, pady=3)

        tk.Label(row3, text="Input JSON File:").pack(side=tk.LEFT, padx=5)
        self.input_file = tk.StringVar()
        tk.Entry(row3, textvariable=self.input_file, width=80).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # --- Button row ---
        btn_frame = tk.Frame(root)
        btn_frame.pack(fill=tk.X, padx=10, pady=5)

        tk.Button(btn_frame, text="Select File", command=self.select_file).pack(side=tk.LEFT, padx=3)
        tk.Button(btn_frame, text="Sequence Input...", command=self.open_sequence_input,
                  bg="lightyellow").pack(side=tk.LEFT, padx=3)
        tk.Button(btn_frame, text="Advanced...", command=self.open_advanced).pack(side=tk.LEFT, padx=3)
        tk.Button(btn_frame, text="Cell Context...", command=self.open_cell_context, bg="lightcyan").pack(side=tk.LEFT, padx=3)

        self.run_button = tk.Button(btn_frame, text="Run Simulation (STRING v12.5)",
                                    command=self.run_simulation, bg="lightgreen")
        self.run_button.pack(side=tk.LEFT, padx=10)

        self.stop_button = tk.Button(btn_frame, text="Stop", command=self.stop_simulation,
                                      state=tk.DISABLED, bg="salmon")
        self.stop_button.pack(side=tk.LEFT, padx=3)

        # --- Log area ---
        self.log_area = scrolledtext.ScrolledText(root, height=30, width=120, wrap=tk.WORD)
        self.log_area.pack(pady=5, padx=10, fill=tk.BOTH, expand=True)

        # --- Progress bar ---
        progress_frame = tk.Frame(root)
        progress_frame.pack(fill=tk.X, padx=10, pady=(2, 5))

        self.progress_var = tk.DoubleVar(value=0.0)
        self.progress_bar = ttk.Progressbar(
            progress_frame,
            orient=tk.HORIZONTAL,
            mode="determinate",
            maximum=100.0,
            variable=self.progress_var
        )
        self.progress_bar.pack(
            side=tk.LEFT,
            fill=tk.X,
            expand=True,
            padx=(0, 8)
        )

        self.progress_label = tk.Label(
            progress_frame,
            text="0.0% | Ready",
            width=68,
            anchor=tk.W
        )
        self.progress_label.pack(side=tk.RIGHT)

        self.progress_started_at = None

        # --- Status bar ---
        self.status = tk.Label(root, text="Ready (STRING v12.5)", relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

        # --- API initialization ---
        # Initialize logging state before StringAPI emits its startup message.
        self.log_file_handle = None
        self.api = StringAPI(logger_print=self._log_threadsafe, network_type="functional")
        self.simulator = None
        self.advanced_config = {}
        self.sequence_data = None  # sequence input data
        self.cell_context_config = {"enabled": False, "mode": "conservative_union", "rna_threshold": 1.0}
        self.cell_context = None

    def _on_network_type_change(self, event=None):
        nt = self.network_type_var.get()
        if nt == "regulatory":
            self.reg_label.config(text="▶ v12.5 official: directional + sign (pos/neg/n/a)")
        elif nt == "physical":
            self.reg_label.config(text="▶ Physical interactions only")
        else:
            self.reg_label.config(text="")

    def select_file(self):
        path = filedialog.askopenfilename(
            title="Select JSON input file",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if path:
            self.input_file.set(path)

    def open_sequence_input(self):
        """Open the sequence input dialog"""
        dlg = SequenceInputDialog(self.root)
        self.root.wait_window(dlg)
        if dlg.result:
            self.sequence_data = dlg.result
            seq_count = len(dlg.result.get("sequences", {}))
            species = dlg.result.get("species_id", "") or "9606 (Human)"
            self._log_threadsafe(f"[Sequence] {seq_count} protein sequences entered (Species: {species})")

            if dlg.result.get("species_id"):
                self.species_var.set(dlg.result["species_id"])

            # Show start protein names from sequences
            names = list(dlg.result["sequences"].keys())
            self._log_threadsafe(f"[Sequence] Start proteins: {', '.join(names[:10])}")
            if len(names) > 10:
                self._log_threadsafe(f"  ... additional {len(names) - 10}")

    def open_advanced(self):
        dlg = AdvancedDialog(self.root, self.advanced_config)
        self.root.wait_window(dlg)
        self.advanced_config = dlg.result
        self._log_threadsafe(f"Advanced config updated: {self.advanced_config}")

    def open_cell_context(self):
        dlg = CellContextDialog(self.root, self.cell_context_config)
        self.root.wait_window(dlg)
        if dlg.result is not None:
            self.cell_context_config = dlg.result
            if dlg.result.get("enabled"):
                self._log_threadsafe(
                    f"[Cell Context] Configured: cell={dlg.result.get('cell_line')}, "
                    f"mode={dlg.result.get('mode')}, RNA={'yes' if dlg.result.get('rna_file') else 'no'}, "
                    f"proteome={'yes' if dlg.result.get('proteome_file') else 'no'}"
                )
            else:
                self._log_threadsafe("[Cell Context] Disabled by user.")

    def _log_threadsafe(self, msg: str):
        if getattr(self, "log_file_handle", None):
            try:
                ts = datetime.now().strftime("%H:%M:%S")
                self.log_file_handle.write(f"[{ts}] {msg}\n")
                self.log_file_handle.flush()
            except Exception:
                pass

        def insert():
            self.log_area.insert(tk.END, msg + "\n")
            self.log_area.see(tk.END)

        if threading.current_thread() is threading.main_thread():
            insert()
        else:
            self.root.after(0, insert)

    @staticmethod
    def _format_eta(seconds):
        if seconds is None or seconds < 0 or not math.isfinite(seconds):
            return "--"
        seconds = int(seconds)
        if seconds < 60:
            return f"{seconds}s"
        minutes, sec = divmod(seconds, 60)
        if minutes < 60:
            return f"{minutes}m {sec:02d}s"
        hours, minutes = divmod(minutes, 60)
        return f"{hours}h {minutes:02d}m"

    def _progress_threadsafe(self, fraction: float, phase: str = "",
                             completed: int = 0, total: int = 0):
        try:
            fraction = max(0.0, min(1.0, float(fraction)))
        except Exception:
            fraction = 0.0

        def update():
            pct = 100.0 * fraction
            self.progress_var.set(pct)

            eta_text = "--"
            if (
                self.progress_started_at is not None
                and fraction > 0.01
                and fraction < 0.999
            ):
                elapsed = max(
                    0.0, time.time() - self.progress_started_at
                )
                eta = elapsed * (1.0 - fraction) / fraction
                eta_text = self._format_eta(eta)
            elif fraction >= 0.999:
                eta_text = "0s"

            pair_text = ""
            if total and total > 0:
                pair_text = f" | {completed:,}/{total:,} pairs"

            phase_text = f" | {phase}" if phase else ""
            if fraction < 0.999:
                eta_part = f" | ETA ~ {eta_text}"
            else:
                eta_part = " | Complete"

            self.progress_label.config(
                text=(
                    f"{pct:5.1f}%"
                    f"{phase_text}"
                    f"{pair_text}"
                    f"{eta_part}"
                )
            )

        if threading.current_thread() is threading.main_thread():
            update()
        else:
            self.root.after(0, update)

    def run_simulation(self):
        fpath = self.input_file.get().strip()
        if not fpath or not os.path.isfile(fpath):
            messagebox.showerror("Error", "Please select a valid JSON file.")
            return

        # Prevent accidental semantic mismatch between a clearly named input file
        # and the GUI network selector (e.g. regulatory input run as physical).
        _lname=os.path.basename(fpath).lower(); _selected=self.network_type_var.get(); _hint=None
        if "regulatory" in _lname: _hint="regulatory"
        elif "physical" in _lname: _hint="physical"
        elif "functional" in _lname: _hint="functional"
        if _hint and _hint != _selected:
            if not messagebox.askyesno("Network type mismatch",
                f"The input filename suggests '{_hint}', but the GUI is set to '{_selected}'.\n\nContinue with '{_selected}' anyway?"):
                return

        self.run_button.config(state=tk.DISABLED)
        self.stop_button.config(state=tk.NORMAL)
        self.status.config(text=f"Running (STRING v12.5 - {self.network_type_var.get()})...")
        self.log_area.delete(1.0, tk.END)

        self.progress_started_at = time.time()
        self.progress_var.set(0.0)
        self.progress_label.config(
            text="0.0% | Initializing | ETA ~ --"
        )

        threading.Thread(target=self._run_worker, args=(fpath,), daemon=True).start()

    def _run_worker(self, fpath: str):
        d, base = os.path.split(fpath)
        base_name = os.path.splitext(base)[0]
        log_path = os.path.join(d, f"{base_name}_STRING_v12_log.txt")
        try:
            self.log_file_handle = open(log_path, 'w', encoding='utf-8')
            self.log_file_handle.write(
                f"=== Protein Pathway Simulator (STRING API v12.5) Run Log ===\n"
                f"Input file: {fpath}\n"
                f"Network type: {self.network_type_var.get()}\n"
                f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"{'=' * 60}\n\n"
            )
            self.log_file_handle.flush()
        except Exception as e:
            self._log_threadsafe(f"[Warning] Failed to open log file: {e}")
            self.log_file_handle = None

        try:
            t0 = time.time()

            # update API settings
            network_type = self.network_type_var.get()
            self.api = StringAPI(
                logger_print=self._log_threadsafe,
                network_type=network_type
            )

            # species settings
            species_id = self.species_var.get().strip()
            if species_id.startswith("STRG"):
                self.api.set_custom_species(species_id)

            # process sequence data
            if self.sequence_data:
                sequences = self.sequence_data.get("sequences", {})
                if sequences:
                    identity_thr = self.sequence_data.get("blast_identity_threshold", 30.0)
                    evalue_thr = self.sequence_data.get("blast_evalue_threshold", 1e-10)
                    self.api.submit_sequences(
                        sequences,
                        blast_identity_threshold=identity_thr,
                        blast_evalue_threshold=evalue_thr
                    )

            with open(fpath, 'r', encoding='utf-8') as f:
                inputs = json.load(f)

            # PhosphoSite manager
            site_manager = None
            if '@' in json.dumps(inputs):
                site_manager = PhosphoSiteManager()
                site_manager.load_from_json(inputs, logger_print=self._log_threadsafe)

            # Build cell-context profile once per run (bulk files are NOT queried inside BFS).
            self.cell_context = CellContextFilter(
                self.cell_context_config,
                logger_print=self._log_threadsafe
            ).build()

            _run_advanced = dict(self.advanced_config or {})
            _run_advanced["temporal_anchor_mode"] = self.temporal_anchor_mode_var.get()
            _run_advanced.setdefault("exclude_opposing_as_bridge", True)
            self.simulator = PathwaySimulator(
                self.api,
                max_depth=self.depth_var.get(),
                advanced=_run_advanced,
                logger_print=self._log_threadsafe,
                site_manager=site_manager,
                progress_callback=self._progress_threadsafe,
                cell_context=self.cell_context
            )

            start_from_json = inputs.get("start_proteins") or []
            ts = inputs.get("time_series") or {}  # legacy only; Temporal v2 uses inputs["intervals"]
            inhib = inputs.get("inhibitor_drugs") or []

            # merge logic for sequence input and JSON start_proteins
            start_from_seq = []
            if self.sequence_data:
                start_from_seq = list(self.sequence_data.get("sequences", {}).keys())

            if start_from_json and start_from_seq:
                # if both exist, use a union (remove duplicates, preserve order)
                seen = set()
                start = []
                for s in start_from_json + start_from_seq:
                    if s not in seen:
                        seen.add(s)
                        start.append(s)
                self._log_threadsafe(
                    f"[Start] JSON({len(start_from_json)}) + sequence input({len(start_from_seq)} items) "
                    f" -> merged {len(start)} start proteins"
                )
                self._log_threadsafe(f"  JSON source: {start_from_json}")
                self._log_threadsafe(f"  Sequence source: {start_from_seq}")
            elif start_from_seq:
                start = start_from_seq
                self._log_threadsafe(
                    f"[Start] Start proteins from sequence input: {len(start)} set: {start}"
                )
            else:
                start = start_from_json
                if start:
                    self._log_threadsafe(f"[Start] Start proteins loaded from JSON: {len(start)} loaded")

            final_paths, name_to_acc_mapping, intermediate_paths = self.simulator.trace_full_pathway(
                start, inputs, inhib
            )

            t1 = time.time()

            self._progress_threadsafe(
                0.95,
                "Saving pathway results"
            )

            out_path = os.path.join(d, f"{base_name}_STRING_v12_results.txt")
            all_out_path = os.path.join(d, f"{base_name}_STRING_v12_all_results.txt")
            phospho_out_path = os.path.join(d, f"{base_name}_STRING_v12_all_Phospho_results.txt")
            temporal_json_path = os.path.join(d, f"{base_name}_STRING_v12_temporal_results.json")
            interval_csv_path = os.path.join(d, f"{base_name}_STRING_v12_interval_summary.csv")
            node_roles_csv_path = os.path.join(d, f"{base_name}_STRING_v12_node_roles.csv")
            edges_csv_path = os.path.join(d, f"{base_name}_STRING_v12_edges.csv")
            cell_context_csv_path = os.path.join(d, f"{base_name}_STRING_v12_cell_context_audit.csv")

            is_regulatory = (network_type in ("regulatory",))

            # ── Build run configuration for result file headers ──
            adv = self.advanced_config or {}
            _hubs = adv.get("hub_blacklist", [])
            _start_str = ", ".join(str(s) for s in start[:10])
            if len(start) > 10:
                _start_str += f" ... (+{len(start)-10})"
            _inhib_str = ", ".join(str(s) for s in inhib) if inhib else "(none)"

            run_config_lines = (
                f"{'='*70}\n"
                f"  RUN CONFIGURATION\n"
                f"{'='*70}\n"
                f"  Run Date/Time        : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"  Input File           : {fpath}\n"
                f"  Network Type         : {network_type}\n"
                f"  STRING API Base      : {self.api.base_url}\n"
                f"  Bridge Depth D       : {self.depth_var.get()} (intermediate nodes; max edges/interval={self.depth_var.get()+1})\n"
                f"  Species ID           : {species_id}\n"
                f"  Min STRING Score     : {adv.get('min_string_score', 700)}\n"
                f"  Degree Cap H         : {adv.get('degree_cap', 50)}\n"
                f"  Segment Beam K       : {adv.get('beam_width', 150)} BFS frontier\n"
                f"  Max Paths/Target N   : {adv.get('max_paths_per_target', 300)} (legacy/non-temporal)\n"
                f"  Temporal Seg Routes  : {adv.get('temporal_segment_route_cap', 5)} per source-target pair\n"
                f"  Temporal Base K      : {adv.get('temporal_beam_width', 0) or adv.get('beam_width', 150)}\n"
                f"  Temporal Adaptive K  : max(Base K, reached target genes)\n"
                f"  Target Preservation  : >=1 route per reached target\n"
                f"  Temporal Routes/Tgt  : {adv.get('temporal_routes_per_target', 3)}\n"
                f"  Temporal Full Dedupe : {adv.get('temporal_full_path_dedupe', True)}\n"
                f"  Regulatory Wt Bonus  : {adv.get('regulatory_weight_bonus', 0.1)}\n"
                f"  Hub Blacklist        : {', '.join(_hubs) if _hubs else '(none)'}\n"
                f"  Collapse by Gene     : {adv.get('collapse_gene_level', True)}\n"
                f"  No Backtracking      : {adv.get('no_immediate_backtrack', True)}\n"
                f"  Dedupe Genes in Path : {adv.get('dedupe_gene_in_path', True)}\n"
                f"  Verbose Logging      : {adv.get('verbose', False)}\n"
            )
            if self.sequence_data:
                _seqs = self.sequence_data.get("sequences", {})
                run_config_lines += (
                    f"  Sequence Input       : Yes ({len(_seqs)} sequences)\n"
                    f"  BLAST Identity Thr   : {self.sequence_data.get('blast_identity_threshold', 30.0)}%\n"
                    f"  BLAST E-value Thr    : {self.sequence_data.get('blast_evalue_threshold', 1e-10)}\n"
                )
            _cc_summary = self.cell_context.summary() if self.cell_context else {"enabled": False}
            run_config_lines += (
                f"  Cell Context Filter  : {_cc_summary.get('enabled', False)}\n"
                f"  Cell Context Mode    : {_cc_summary.get('mode', 'off')}\n"
                f"  Cell Line            : {_cc_summary.get('cell_line', '') or '(none)'}\n"
                f"  RNA Bulk File        : {_cc_summary.get('rna_file') or '(none)'}\n"
                f"  Proteomics Bulk File : {_cc_summary.get('proteome_file') or '(none)'}\n"
                f"  RNA Threshold        : {_cc_summary.get('rna_threshold', 1.0)}\n"
            )
            run_config_lines += (
                f"  Start Proteins       : {_start_str}\n"
                f"  Input Schema         : {self.simulator.input_schema_label}\n"
                f"  Temporal Target Mode : {self.temporal_anchor_mode_var.get()}\n"
                f"  Time Points/Intervals: {len(inputs.get('intervals', [])) if inputs.get('intervals') else len(ts)}\n"
                f"  Inhibitor Drugs      : {_inhib_str}\n"
                f"  Total Elapsed Time   : {t1 - t0:.1f} sec\n"
                f"{'='*70}\n\n"
            )

            # save complete paths
            try:
                with open(out_path, "w", encoding="utf-8") as f:
                    f.write(f"*** Protein Pathway Simulator Results (STRING API v12.5 - {network_type}) ***\n\n")
                    f.write(run_config_lines)

                    f.write("*** Protein Name to Gene Mapping Table ***\n")
                    f.write("| Original Name | Gene Symbol |\n")
                    f.write("|---------------|-------------|\n")
                    for name, gene in sorted((name_to_acc_mapping or {}).items()):
                        f.write(f"| {name} | {gene} |\n")

                    f.write(f"\n*** Complete Path Results (STRING API v12.5 - weighted, {network_type} mode) ***\n")
                    f.write(f"(Total elapsed time: {t1 - t0:.1f} sec)\n\n")

                    valid = [p for p in final_paths if len(p) > 1]

                    if not valid:
                        f.write("No valid connecting path was found.\n")
                    else:
                        f.write(f"Total {len(valid)} complete paths were found.\n\n")
                        for i, path in enumerate(valid, 1):
                            genes = [self.simulator._gene(x) for x in path]
                            if is_regulatory:
                                arrow_str = self.simulator._format_pathway_arrow(genes)
                                f.write(f"Path #{i}: {arrow_str}\n")
                            else:
                                f.write(f"Path #{i}: {' -> '.join(genes)}\n")

                    # Regulatory relationship table (regulatory mode)
                    if is_regulatory and self.simulator.regulatory_info:
                        f.write(f"\n\n*** Regulatory Interaction Details (Regulatory Interactions) ***\n")
                        f.write("| Source | Target | Action | Directed | Score(A→B) | Score(B→A) |\n")
                        f.write("|--------|--------|--------|----------|------------|------------|\n")
                        for (src, tgt), info in sorted(self.simulator.regulatory_info.items()):
                            action = info.get("action", "?")
                            directed = "✓" if info.get("is_directed") else ""
                            s_ab = info.get("score_A_to_B", 0)
                            s_ba = info.get("score_B_to_A", 0)
                            if isinstance(s_ab, float):
                                s_ab = f"{s_ab:.3f}"
                            if isinstance(s_ba, float):
                                s_ba = f"{s_ba:.3f}"
                            f.write(f"| {src} | {tgt} | {action} | {directed} | {s_ab} | {s_ba} |\n")

            except Exception as e:
                self._log_threadsafe(f"Error saving results: {e}")

            # save all paths
            try:
                with open(all_out_path, "w", encoding="utf-8") as f:
                    f.write(f"*** All Path Results (Complete + Intermediate) - STRING API v12.5 ({network_type}) ***\n")
                    f.write(run_config_lines)
                    f.write(f"(Total elapsed time: {t1 - t0:.1f} sec)\n\n")

                    valid_complete = [p for p in final_paths if len(p) > 1]
                    f.write(f"=== Complete Paths ({len(valid_complete)} items) ===\n")
                    f.write("(Paths that reached the end through all time points)\n\n")

                    if valid_complete:
                        for i, path in enumerate(valid_complete, 1):
                            genes = [self.simulator._gene(x) for x in path]
                            if is_regulatory:
                                arrow_str = self.simulator._format_pathway_arrow(genes)
                                f.write(f"Complete Path #{i} (Length {len(path)}): {arrow_str}\n")
                            else:
                                f.write(f"Complete Path #{i} (Length {len(path)}): {' -> '.join(genes)}\n")
                    else:
                        f.write("No complete paths were found.\n")

                    f.write(f"\n\n=== Intermediate Paths ({len(intermediate_paths)} items) ===\n")

                    if intermediate_paths:
                        paths_by_length = {}
                        for path in intermediate_paths:
                            if len(path) > 1:
                                length = len(path)
                                if length not in paths_by_length:
                                    paths_by_length[length] = []
                                paths_by_length[length].append(path)

                        for length in sorted(paths_by_length.keys(), reverse=True):
                            paths = paths_by_length[length]
                            f.write(f"\n--- Length {length} ({len(paths)} items) ---\n")

                            for i, path in enumerate(paths[:100], 1):
                                genes = [self.simulator._gene(x) for x in path]
                                if is_regulatory:
                                    arrow_str = self.simulator._format_pathway_arrow(genes)
                                    f.write(f"  {i}. {arrow_str}\n")
                                else:
                                    f.write(f"  {i}. {' -> '.join(genes)}\n")

                            if len(paths) > 100:
                                f.write(f"  ... additional {len(paths) - 100} more\n")
                    else:
                        f.write("No intermediate paths.\n")

                    f.write("\n\n*** Path Statistics Summary ***\n")
                    f.write(f"Network type: {network_type}\n")
                    f.write(f"Complete path count: {len(valid_complete)}\n")
                    f.write(f"Intermediate path count: {len(intermediate_paths)}\n")
                    f.write(f"Total path count: {len(valid_complete) + len(intermediate_paths)}\n")

                    if is_regulatory and self.simulator.regulatory_info:
                        f.write(f"\n*** Regulatory Relationship Statistics ***\n")
                        act = sum(1 for v in self.simulator.regulatory_info.values()
                                  if v.get("action") in ("activation", "catalysis", "expression"))
                        inh = sum(1 for v in self.simulator.regulatory_info.values()
                                  if v.get("action") in ("inhibition",))
                        dr = sum(1 for v in self.simulator.regulatory_info.values() if v.get("is_directed"))
                        f.write(f"Collected regulatory relationships: {len(self.simulator.regulatory_info)}\n")
                        f.write(f"Activation edges: {act}\n")
                        f.write(f"Inhibition edges: {inh}\n")
                        f.write(f"Directional edges: {dr}\n")

            except Exception as e:
                self._log_threadsafe(f"Error saving all results: {e}")

            # ── FINAL target-preserving structured Temporal v2 outputs ──
            if self.simulator.input_schema_label.startswith("Temporal v2"):
                try:
                    temporal_payload = {
                        "schema_version": "3.4-final",
                        "model": "BFS-BS_Q_conditioned_target_preserving_adaptive_temporal_paths",
                        "source_input": fpath,
                        "network_type": network_type,
                        "temporal_target_mode": self.temporal_anchor_mode_var.get(),
                        "bridge_depth_D": self.depth_var.get(),
                        "bridge_depth_semantics": "maximum intermediate prior-network bridge nodes per interval",
                        "max_edges_per_interval": self.depth_var.get() + 1,
                        "baseline_used_as_activation_target": False,
                        "opposing_bridge_policy": "Q opposing-only excluded from active target/bridge pool; mixed genes retained",
                        "condition_functional_pool": getattr(self.simulator, "condition_functional_pool_info", {}),
                        "cell_context_filter": self.cell_context.summary() if self.cell_context else {"enabled": False},
                        "complete_path_count": len(self.simulator.temporal_path_records),
                        "paths": self.simulator.temporal_path_records,
                        "interval_summary": self.simulator.temporal_interval_summary,
                        "route_preservation_policy": {
                            "best_per_source_target_pruning": False,
                            "best_per_reached_target_pruning": False,
                            "target_stratified_interval_beam": True,
                            "target_preserving_adaptive_temporal_beam": True,
                            "adaptive_temporal_beam_rule":
                                "K_effective=max(K_temporal_base,reached_target_genes_after_exact_dedup)",
                            "reached_target_minimum_route_guarantee": 1,
                            "fixed_interval_beam_pruning": False,
                            "exact_gene_path_deduplication_only": False,
                            "full_temporal_path_gene_revisit_prevention": bool(
                                self.simulator.adv.get(
                                    "temporal_full_path_dedupe", True
                                )
                            ),
                            "temporal_segment_route_cap": int(
                                self.simulator.adv.get(
                                    "temporal_segment_route_cap", 5
                                )
                            ),
                            "temporal_base_beam_width": (
                                int(
                                    self.simulator.adv.get(
                                        "temporal_beam_width", 0
                                    ) or 0
                                )
                                or int(
                                    self.simulator.adv.get(
                                        "beam_width", 150
                                    )
                                )
                            ),
                            "temporal_effective_beam_width":
                                "reported per interval in pruning_summary and interval_summary",
                            "temporal_routes_per_target": int(
                                self.simulator.adv.get(
                                    "temporal_routes_per_target", 3
                                )
                            ),
                            "search_limits_retained": [
                                "Bridge Depth D",
                                "Min STRING Score",
                                "Degree Cap H",
                                "Segment Beam Width K",
                                "Temporal Beam K",
                                "Temporal segment route cap"
                            ],
                            "phosphosite_cartesian_expansion": False
                        },
                        "pruning_summary": self.simulator.temporal_pruning_summary,
                        "generated_at": datetime.now().isoformat(timespec="seconds")
                    }
                    with open(temporal_json_path, "w", encoding="utf-8") as tf:
                        json.dump(temporal_payload, tf, ensure_ascii=False, indent=2)

                    if self.simulator.temporal_interval_summary:
                        keys = list(self.simulator.temporal_interval_summary[0].keys())
                        with open(interval_csv_path, "w", newline="", encoding="utf-8-sig") as cf:
                            w = csv.DictWriter(cf, fieldnames=keys)
                            w.writeheader(); w.writerows(self.simulator.temporal_interval_summary)

                    role_rows = []
                    edge_rows = []
                    for rec in self.simulator.temporal_path_records:
                        pid = rec["path_id"]
                        for seg in rec.get("segments", []):
                            role_rows.append({"path_id": pid, "interval": seg["interval"], "gene": seg["source_anchor"], "role": "source_anchor"})
                            for b in seg.get("bridge_nodes", []):
                                role_rows.append({"path_id": pid, "interval": seg["interval"], "gene": b, "role": "prior_bridge"})
                            role_rows.append({"path_id": pid, "interval": seg["interval"], "gene": seg["target_anchor"], "role": "target_" + seg.get("target_role", "measured")})
                            chain = [seg["source_anchor"]] + list(seg.get("bridge_nodes", [])) + [seg["target_anchor"]]
                            for ei in range(len(chain)-1):
                                edge_rows.append({"path_id": pid, "interval": seg["interval"], "edge_index": ei+1, "source": chain[ei], "target": chain[ei+1], "segment_bridge_depth": seg.get("bridge_depth", 0), "segment_cost": seg.get("segment_cost", 0)})
                    with open(node_roles_csv_path, "w", newline="", encoding="utf-8-sig") as cf:
                        w = csv.DictWriter(cf, fieldnames=["path_id","interval","gene","role"]); w.writeheader(); w.writerows(role_rows)
                    with open(edges_csv_path, "w", newline="", encoding="utf-8-sig") as cf:
                        w = csv.DictWriter(cf, fieldnames=["path_id","interval","edge_index","source","target","segment_bridge_depth","segment_cost"]); w.writeheader(); w.writerows(edge_rows)
                    self._log_threadsafe(f"  - Temporal structured results: {temporal_json_path}")
                    self._log_threadsafe(f"  - Interval summary: {interval_csv_path}")
                    self._log_threadsafe(f"  - Node roles: {node_roles_csv_path}")
                    self._log_threadsafe(f"  - Edges: {edges_csv_path}")
                except Exception as e:
                    self._log_threadsafe(f"[Warning] Failed to save Stage 3 structured outputs: {e}")

            if self.cell_context and self.cell_context.enabled:
                try:
                    self.cell_context.write_audit_csv(cell_context_csv_path)
                    self._log_threadsafe(f"  - Cell-context bridge audit: {cell_context_csv_path}")
                    cc=self.cell_context.summary()
                    ds=cc.get("decision_stats", {})
                    self._log_threadsafe(
                        f"[Cell Context Summary] checked={ds.get('bridge_candidates_checked',0)}, "
                        f"allowed={ds.get('allowed',0)}, blocked={ds.get('blocked',0)}, "
                        f"unknown-retained={ds.get('unknown_retained',0)}"
                    )
                except Exception as e:
                    self._log_threadsafe(f"[Warning] Failed to save cell-context audit CSV: {e}")

            # Reproducibility snapshot: exact neighbors returned to BFS (cache or live API).
            try:
                snapshot_dir = os.path.join(d, f"{base_name}_STRING_v12_5_snapshot_{network_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{time.time_ns() % 1000000000:09d}")
                os.makedirs(snapshot_dir, exist_ok=True)
                input_bytes = Path(fpath).read_bytes()
                input_sha = hashlib.sha256(input_bytes).hexdigest()
                with self.api.snapshot_lock:
                    query_snapshot = {q: dict(v) for q, v in self.api.snapshot_neighbors.items()}
                snap_rows = []
                for query in sorted(query_snapshot):
                    entry = query_snapshot[query]
                    for rank, n in enumerate(entry["neighbors"], 1):
                        snap_rows.append({
                            "query_gene": query, "rank": rank, "origin": entry["origin"],
                            "partner_gene": n.get("gene", ""), "score_1000": n.get("score", ""),
                            "weight": n.get("weight", ""), "source": n.get("source", ""),
                            "target": n.get("target", ""), "is_directed": n.get("is_directed", False),
                            "action": n.get("action", ""), "score_A_to_B": n.get("score_A_to_B", ""),
                            "score_B_to_A": n.get("score_B_to_A", "")})
                columns = ["query_gene", "rank", "origin", "partner_gene", "score_1000", "weight",
                           "source", "target", "is_directed", "action", "score_A_to_B", "score_B_to_A"]
                edge_name = f"{network_type}_neighbors.csv"
                edge_path = os.path.join(snapshot_dir, edge_name)
                with open(edge_path, "w", encoding="utf-8", newline="") as sf:
                    writer = csv.DictWriter(sf, fieldnames=columns)
                    writer.writeheader(); writer.writerows(snap_rows)
                # Full precision, lossless API-transformed neighbor records and ordering.
                neighbors_path = os.path.join(snapshot_dir, "queried_neighbors.json")
                with open(neighbors_path, "w", encoding="utf-8") as sf:
                    json.dump(query_snapshot, sf, ensure_ascii=False, indent=2, allow_nan=False)
                metadata = {
                    "string_version": "12.5", "release_channel": "official", "api_base": self.api.base_url,
                    "network_type": network_type, "endpoint": "/json/interaction_partners",
                    "species": self.api._get_species(), "api_required_score": 0,
                    "created_at": datetime.now().astimezone().isoformat(),
                    "input_json_sha256": input_sha, "input_json_name": os.path.basename(fpath),
                    "queried_gene_count": len(query_snapshot), "neighbor_records": len(snap_rows),
                    "origins": dict(Counter(x["origin"] for x in query_snapshot.values())),
                    "neighbors_sha256": hashlib.sha256(Path(neighbors_path).read_bytes()).hexdigest(),
                    "neighbor_csv_sha256": hashlib.sha256(Path(edge_path).read_bytes()).hexdigest(),
                    "output_complete_paths": len(final_paths), "search_parameters": dict(_run_advanced),
                    "min_string_score": _run_advanced.get("min_string_score", 700),
                    "degree_cap_H": _run_advanced.get("degree_cap", 50),
                    "bridge_depth_D": self.depth_var.get(),
                    "cache_file": self.api.cache.cache_file,
                    "note": "Includes all returned candidates before minimum-score and beam pruning; cache entries are version-isolated."
                }
                with open(os.path.join(snapshot_dir, "metadata.json"), "w", encoding="utf-8") as sf:
                    json.dump(metadata, sf, ensure_ascii=False, indent=2, default=str)
                with open(os.path.join(snapshot_dir, "input_sha256.txt"), "w", encoding="utf-8") as sf:
                    sf.write(input_sha + "  " + os.path.basename(fpath) + "\n")
                self._log_threadsafe(f"[Snapshot] {snapshot_dir}: {len(query_snapshot)} queries, {len(snap_rows)} neighbor records")
            except Exception as e:
                self._log_threadsafe(f"[Snapshot] ERROR: {e}")
                raise RuntimeError("STRING snapshot save failed; do not use this run for publication") from e

            # Persist STRING cache for later runs.
            try:
                self.api.cache.save()
                self._log_threadsafe("[Cache] STRING cache saved.")
            except Exception as e:
                self._log_threadsafe(
                    f"[Cache] Failed to save STRING cache: {e}"
                )

            # print STRING statistics
            self.api.print_stats()

            # ── generate phosphorylation-site path file ──
            phospho_complete: List[Tuple[str, ...]] = []
            phospho_intermediate: List[Tuple[str, ...]] = []

            if site_manager and site_manager.protein_site_combos:
                self._log_threadsafe("\n[Phospho] Generating phosphorylation-site path results...")

                valid_complete = [p for p in final_paths if len(p) > 1]
                valid_intermediate = [p for p in intermediate_paths if len(p) > 1]

                if self.simulator.input_schema_label.startswith("Temporal v2"):
                    phospho_complete = generate_temporal_endpoint_phospho_paths(self.simulator.temporal_path_records)
                    # Intermediate temporal paths are intentionally not Cartesian-expanded.
                    # The gene-level intermediate file remains authoritative for partial routes.
                    phospho_intermediate = []
                    self._log_threadsafe("[Phospho FINAL] Temporal endpoint annotation mode; Cartesian site expansion disabled.")
                else:
                    phospho_complete = generate_phospho_paths(valid_complete, site_manager)
                    phospho_intermediate = generate_phospho_paths(valid_intermediate, site_manager)

                phospho_complete.sort(key=lambda t: (len(t), '->'.join(t)))
                phospho_intermediate.sort(key=lambda t: (len(t), '->'.join(t)))

                self._log_threadsafe(
                    f"[Phospho] Complete path variants: {len(phospho_complete)}, "
                    f"Intermediate path variants: {len(phospho_intermediate)} generated"
                )
            else:
                self._log_threadsafe(
                    "[Phospho] No phosphorylation-site (@XXX) data was found, so the Phospho file will be saved empty."
                )

            try:
                write_phospho_results(phospho_out_path, phospho_complete, phospho_intermediate, t1 - t0)
            except Exception as e:
                self._log_threadsafe(f"Error saving phospho results: {e}")

            self._progress_threadsafe(
                1.0,
                "Complete"
            )
            self.root.after(0, lambda: self._log_threadsafe(f"Simulation complete."))
            self.root.after(0, lambda: self._log_threadsafe(f"  - Complete path results: {out_path}"))
            self.root.after(0, lambda: self._log_threadsafe(f"  - All path results: {all_out_path}"))
            self.root.after(0, lambda: self._log_threadsafe(f"  - Phospho path results: {phospho_out_path}"))
            self.root.after(0, lambda: self._log_threadsafe(f"  - Run log: {log_path}"))
            self.root.after(0, lambda: self.status.config(
                text=f"Complete ({network_type}): {len([p for p in final_paths if len(p) > 1])} complete, "
                     f"{len(intermediate_paths)} intermediate (took {t1 - t0:.1f}s)"
            ))

        except Exception as e:
            err_msg = f"{type(e).__name__}: {e}\n" + traceback.format_exc()
            err_str = str(e)
            err_type = type(e).__name__
            self.root.after(0, lambda msg=err_msg, error=err_str, etype=err_type: (
                self._log_threadsafe(f"Error during simulation ({etype}): {error}"),
                messagebox.showerror("Error", msg)
            ))

        finally:
            if self.log_file_handle:
                try:
                    self.log_file_handle.write(
                        f"\n{'=' * 60}\n"
                        f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                        f"=== End of log ===\n"
                    )
                    self.log_file_handle.close()
                except Exception:
                    pass
                self.log_file_handle = None
            self.root.after(0, self.reset_ui)

    def stop_simulation(self):
        if self.simulator:
            self.simulator.set_stop_flag(True)
        self._log_threadsafe("Stop requested. Halting simulation...")
        self.stop_button.config(state=tk.DISABLED)
        self.status.config(text="Stopping...")
        current_pct = self.progress_var.get()
        self.progress_label.config(
            text=f"{current_pct:.1f}% | Stopping..."
        )
        # _run_worker finally block resets the UI after the worker exits.

    def reset_ui(self):
        self.run_button.config(state=tk.NORMAL)
        self.stop_button.config(state=tk.DISABLED)
        self.status.config(text="Stopped or Complete")


# ========== Main ==========

def main():
    root = tk.Tk()
    app = PathwayGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()