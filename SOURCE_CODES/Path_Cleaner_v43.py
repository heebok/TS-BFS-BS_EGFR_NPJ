#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BFS Path Cleaner — Stage 4.3 Route-Preserving
=========================================

Purpose
-------
Stage 4 is a NON-DESTRUCTIVE post-processing and validation layer for the
Stage 3.2 temporal BFS-BS output.

Design principles
-----------------
1. Prefer *_STRING_v12_temporal_results.json when available.
2. Preserve biologically distinct alternative routes.
3. Never truncate a Temporal v2 path merely because a gene appears again in
   a later interval. Temporal identity is evaluated as (gene, interval, role).
4. Remove only exact duplicate temporal records (keeping the lower-cost copy).
5. Reject/flag structural inconsistencies such as:
   - broken source -> target continuity between intervals,
   - bridge_depth inconsistent with bridge_nodes,
   - edge_count inconsistent with segment chain,
   - bridge depth greater than configured D,
   - within-segment repeated-node cycles,
   - interval order going backward,
   - mismatch between stored genes and reconstructed segment genes.
6. Condition-pool excluded genes are audited. An explicitly forced start
   protein is allowed only as the source of the first interval.
7. Legacy *_results.txt remains supported. For Temporal v2 text files the
   fallback is conservative: exact duplicate gene chains only. No gene-cycle
   truncation is performed without interval metadata.

Outputs (Temporal structured mode)
----------------------------------
  *_STAGE4_clean_temporal.json
  *_STAGE4_clean_results.txt
  *_STAGE4_topology.txt
  *_STAGE4_topology_results.txt
  *_STAGE4_audit.csv
  *_STAGE4_nodes.csv
  *_STAGE4_edges.csv

Outputs (legacy text mode)
--------------------------
  *_STAGE4_clean_results.txt
  *_STAGE4_topology.txt
  *_STAGE4_topology_results.txt
  *_STAGE4_audit.csv

All program UI text is English.
"""

from __future__ import annotations

import csv
import json
import os
import re
import threading
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

VERSION = "4.3-route-preserving-zero-edge-safe"


# =============================================================================
# Helpers
# =============================================================================

def norm_gene(x: Any) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    if not s:
        return ""
    # Stage 3 records should already be genes, but tolerate token input.
    if "@" in s:
        s = s.split("@", 1)[0]
    if s.startswith("-"):
        s = s[1:]
    return s.strip()


def uniq_preserve(seq: Iterable[Any]) -> list:
    seen = set()
    out = []
    for x in seq:
        key = json.dumps(x, sort_keys=True, ensure_ascii=False) if isinstance(x, (dict, list)) else x
        if key not in seen:
            seen.add(key)
            out.append(x)
    return out


def parse_interval_label(label: str) -> Tuple[Optional[float], Optional[float]]:
    """Parse labels such as '0min-5min', '0-5', '(0,5]' into numeric bounds."""
    if not label:
        return None, None
    nums = re.findall(r"-?\d+(?:\.\d+)?", str(label))
    if len(nums) < 2:
        return None, None
    try:
        return float(nums[0]), float(nums[1])
    except Exception:
        return None, None


def is_zero_edge_temporal_continuation(seg: Dict[str, Any]) -> bool:
    """Return True for a temporal anchor carried forward without a PPI edge.

    Example: source_anchor == target_anchor, no bridge nodes, edge_count == 0.
    This is a temporal persistence/continuation record, NOT a self-loop.
    """
    src = norm_gene(seg.get("source_anchor"))
    tgt = norm_gene(seg.get("target_anchor"))
    bridges = [norm_gene(x) for x in (seg.get("bridge_nodes") or []) if norm_gene(x)]
    if not src or not tgt or src != tgt or bridges:
        return False
    ec = seg.get("edge_count")
    if ec is None:
        return True
    try:
        return int(ec) == 0
    except Exception:
        return False


def segment_chain(seg: Dict[str, Any]) -> List[str]:
    """Canonical gene chain for one temporal segment.

    Zero-edge temporal continuations are represented as [source] rather than
    [source, source]. This prevents a legitimate carry-forward anchor from
    being misclassified as a one-edge self-loop or within-segment cycle.
    """
    src = norm_gene(seg.get("source_anchor"))
    bridges = [norm_gene(x) for x in (seg.get("bridge_nodes") or [])]
    bridges = [x for x in bridges if x]
    tgt = norm_gene(seg.get("target_anchor"))

    if src and tgt and src == tgt and not bridges:
        return [src]

    out = []
    if src:
        out.append(src)
    out.extend(bridges)
    if tgt:
        out.append(tgt)
    return out


def reconstruct_path_genes(segments: Sequence[Dict[str, Any]]) -> List[str]:
    genes: List[str] = []
    for seg in segments:
        ch = segment_chain(seg)
        if not ch:
            continue
        if not genes:
            genes.extend(ch)
        elif genes[-1] == ch[0]:
            genes.extend(ch[1:])
        else:
            # Keep both so mismatch is visible to validator rather than hidden.
            genes.extend(ch)
    return genes


def temporal_signature(rec: Dict[str, Any]) -> Tuple:
    """Exact interval-aware identity. Different bridges remain different routes."""
    sig = []
    for seg in rec.get("segments", []) or []:
        sig.append((
            str(seg.get("interval", "")),
            norm_gene(seg.get("source_anchor")),
            tuple(norm_gene(x) for x in (seg.get("bridge_nodes") or [])),
            norm_gene(seg.get("target_anchor")),
            str(seg.get("target_role", "")),
        ))
    return tuple(sig)


def topology_signature(rec: Dict[str, Any]) -> Tuple:
    """
    Stage-4 topology identity remains interval-aware and route-preserving.
    Target site tokens/costs are annotations and do not define topology.
    """
    sig = []
    for seg in rec.get("segments", []) or []:
        sig.append((
            str(seg.get("interval", "")),
            norm_gene(seg.get("source_anchor")),
            tuple(norm_gene(x) for x in (seg.get("bridge_nodes") or [])),
            norm_gene(seg.get("target_anchor")),
        ))
    return tuple(sig)


def path_cost(rec: Dict[str, Any]) -> float:
    try:
        return float(rec.get("total_cost", 0.0))
    except Exception:
        return 0.0


def drawing_topology_signature(rec: Dict[str, Any]) -> Tuple[str, ...]:
    # Flattened gene-chain identity used ONLY for topology drawing export.
    # Stage-4 validation remains interval-aware. No cycle truncation is done.
    return tuple(
        norm_gene(x)
        for x in (rec.get("genes") or [])
        if norm_gene(x)
    )


# =============================================================================
# Temporal validation
# =============================================================================

@dataclass
class AuditResult:
    path_id: Any
    valid: bool = True
    status: str = "kept"
    reasons: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    original_genes: List[str] = field(default_factory=list)
    reconstructed_genes: List[str] = field(default_factory=list)
    total_cost: float = 0.0

    def fail(self, reason: str):
        self.valid = False
        self.status = "rejected"
        self.reasons.append(reason)

    def warn(self, warning: str):
        self.warnings.append(warning)


def interval_rank_map(payload: Dict[str, Any]) -> Dict[str, int]:
    ranks: Dict[str, int] = {}
    summary = payload.get("interval_summary") or []
    for i, row in enumerate(summary):
        lab = str(row.get("interval") or row.get("id") or "").strip()
        if lab and lab not in ranks:
            ranks[lab] = i
    if ranks:
        return ranks

    # Fallback from first path segment order.
    for rec in payload.get("paths", []) or []:
        for seg in rec.get("segments", []) or []:
            lab = str(seg.get("interval", "")).strip()
            if lab and lab not in ranks:
                ranks[lab] = len(ranks)
    return ranks


def validate_temporal_record(
    rec: Dict[str, Any],
    D: Optional[int],
    ranks: Dict[str, int],
    condition_excluded: set,
    explicit_start_genes: set,
) -> AuditResult:
    ar = AuditResult(path_id=rec.get("path_id", ""), total_cost=path_cost(rec))
    ar.original_genes = [norm_gene(x) for x in (rec.get("genes") or []) if norm_gene(x)]
    segs = rec.get("segments") or []

    if not isinstance(segs, list) or not segs:
        ar.fail("No temporal segments")
        return ar

    prev_target: Optional[str] = None
    prev_rank: Optional[int] = None

    for si, seg in enumerate(segs):
        if not isinstance(seg, dict):
            ar.fail(f"Segment {si+1} is not an object")
            continue

        label = str(seg.get("interval", "")).strip()
        src = norm_gene(seg.get("source_anchor"))
        tgt = norm_gene(seg.get("target_anchor"))
        bridges = [norm_gene(x) for x in (seg.get("bridge_nodes") or []) if norm_gene(x)]
        zero_edge_continuation = bool(src and tgt and src == tgt and not bridges)
        chain = ([src] if src else []) + bridges + ([tgt] if tgt else [])
        canonical_chain = [src] if zero_edge_continuation and src else chain

        if zero_edge_continuation:
            ar.warn(f"Segment {si+1}: zero-edge temporal continuation retained ({src})")

        if not label:
            ar.fail(f"Segment {si+1}: missing interval label")
        if not src:
            ar.fail(f"Segment {si+1}: missing source_anchor")
        if not tgt:
            ar.fail(f"Segment {si+1}: missing target_anchor")

        # Interval order
        rank = ranks.get(label)
        if rank is None:
            a, b = parse_interval_label(label)
            if a is not None and b is not None and b < a:
                ar.fail(f"Segment {si+1}: backward interval {label}")
        else:
            if prev_rank is not None and rank <= prev_rank:
                ar.fail(f"Segment {si+1}: non-forward interval order ({label})")
            prev_rank = rank

        # Sequential temporal anchor continuity
        if prev_target is not None and src != prev_target:
            ar.fail(
                f"Segment {si+1}: broken temporal continuity; source {src} != previous target {prev_target}"
            )
        prev_target = tgt

        # D semantics: intermediate bridge node count.
        bd = seg.get("bridge_depth")
        if bd is not None:
            try:
                if int(bd) != len(bridges):
                    ar.fail(
                        f"Segment {si+1}: bridge_depth={bd} but bridge_nodes={len(bridges)}"
                    )
            except Exception:
                ar.fail(f"Segment {si+1}: invalid bridge_depth={bd!r}")
        if D is not None and len(bridges) > D:
            ar.fail(f"Segment {si+1}: {len(bridges)} bridge nodes exceeds D={D}")

        ec = seg.get("edge_count")
        expected_edges = 0 if zero_edge_continuation else max(0, len(chain) - 1)
        if ec is not None:
            try:
                if int(ec) != expected_edges:
                    ar.fail(
                        f"Segment {si+1}: edge_count={ec} but reconstructed edges={expected_edges}"
                    )
            except Exception:
                ar.fail(f"Segment {si+1}: invalid edge_count={ec!r}")

        # Repetition inside ONE real PPI segment is a cycle/revisit under the D model.
        # A zero-edge source==target continuation is explicitly exempt because it
        # represents the same measured anchor carried across an interval boundary.
        if len(canonical_chain) != len(set(canonical_chain)):
            ar.fail(f"Segment {si+1}: repeated gene inside a single interval PPI segment")

        # Q-derived condition pool audit. Explicit first source is the only exception.
        for ci, gene in enumerate(chain):
            if not gene or gene not in condition_excluded:
                continue
            is_explicit_first_source = (
                si == 0 and ci == 0 and gene in explicit_start_genes
            )
            if is_explicit_first_source:
                ar.warn(
                    f"Explicit experimental source {gene} overrides baseline pool exclusion only at interval 1 initiation"
                )
            else:
                ar.fail(
                    f"Condition-pool excluded gene {gene} appears in segment {si+1} ({label})"
                )

    ar.reconstructed_genes = reconstruct_path_genes(segs)

    # IMPORTANT: repeated genes across different intervals are NOT automatically a cycle.
    # Compare exact flattened representation only as a structural consistency check.
    if ar.original_genes and ar.reconstructed_genes != ar.original_genes:
        ar.fail(
            "Stored genes do not match reconstruction from temporal segments"
        )

    return ar


def clean_temporal_payload(payload: Dict[str, Any]) -> Tuple[Dict[str, Any], List[AuditResult], Dict[str, Any]]:
    paths = payload.get("paths") or []
    if not isinstance(paths, list):
        raise ValueError("Temporal JSON field 'paths' must be a list.")

    D = payload.get("bridge_depth_D")
    try:
        D = int(D) if D is not None else None
    except Exception:
        D = None

    ranks = interval_rank_map(payload)

    pool = payload.get("condition_functional_pool") or {}
    excluded = set(norm_gene(x) for x in (pool.get("excluded_low_activity_genes") or []) if norm_gene(x))

    # Explicit starts are best recovered from path first nodes because Stage 3.2 structured
    # output may not separately repeat input start_proteins.
    explicit_starts = set()
    for rec in paths:
        genes = rec.get("genes") or []
        if genes:
            g = norm_gene(genes[0])
            if g:
                explicit_starts.add(g)

    audits: List[AuditResult] = []
    valid_records: List[Dict[str, Any]] = []

    for rec in paths:
        if not isinstance(rec, dict):
            ar = AuditResult(path_id="", valid=False, status="rejected")
            ar.reasons.append("Path record is not an object")
            audits.append(ar)
            continue
        ar = validate_temporal_record(rec, D, ranks, excluded, explicit_starts)
        audits.append(ar)
        if ar.valid:
            valid_records.append(rec)

    # Exact temporal duplicate elimination only. No best-per-target or best-topology pruning.
    by_sig: Dict[Tuple, Dict[str, Any]] = {}
    duplicate_ids: Dict[Any, Any] = {}
    for rec in valid_records:
        sig = temporal_signature(rec)
        if sig not in by_sig:
            by_sig[sig] = rec
        else:
            old = by_sig[sig]
            if path_cost(rec) < path_cost(old):
                duplicate_ids[old.get("path_id")] = rec.get("path_id")
                by_sig[sig] = rec
            else:
                duplicate_ids[rec.get("path_id")] = old.get("path_id")

    kept = list(by_sig.values())
    kept.sort(key=lambda r: (path_cost(r), temporal_signature(r)))

    # Update audit records for exact duplicates.
    for ar in audits:
        if ar.valid and ar.path_id in duplicate_ids:
            ar.status = "exact_duplicate_removed"
            ar.warn(f"Exact temporal duplicate of retained path #{duplicate_ids[ar.path_id]}")

    # Renumber clean paths while preserving source id.
    cleaned_paths = []
    for new_id, rec in enumerate(kept, 1):
        nr = json.loads(json.dumps(rec))  # deep copy without extra dependency
        nr["source_path_id"] = rec.get("path_id")
        nr["path_id"] = new_id
        cleaned_paths.append(nr)

    topo_count = len({topology_signature(r) for r in kept})
    rejected = sum(1 for a in audits if not a.valid)
    duplicates = sum(1 for a in audits if a.status == "exact_duplicate_removed")

    zero_edge_segments = sum(
        1
        for rec in paths if isinstance(rec, dict)
        for seg in (rec.get("segments") or [])
        if isinstance(seg, dict) and is_zero_edge_temporal_continuation(seg)
    )
    reject_reason_counts = Counter(
        reason
        for a in audits if not a.valid
        for reason in a.reasons
    )

    stats = {
        "input_path_count": len(paths),
        "structurally_valid_before_dedup": len(valid_records),
        "rejected_path_count": rejected,
        "exact_duplicate_removed_count": duplicates,
        "clean_path_count": len(cleaned_paths),
        "unique_interval_aware_topologies": topo_count,
        "zero_edge_temporal_continuation_count": zero_edge_segments,
        "condition_pool_excluded_gene_count": len(excluded),
        "interval_count": len(ranks),
        "rejection_reason_count": len(reject_reason_counts),
    }

    out = dict(payload)
    out["schema_version"] = "4.0"
    out["model"] = "BFS-BS_temporal_paths_stage4_validated"
    out["stage4"] = {
        "version": VERSION,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "policy": {
            "temporal_identity": "(gene, interval, role)",
            "cross_interval_gene_repetition_auto_removed": False,
            "within_segment_cycles_allowed": False,
            "exact_temporal_duplicates_removed": True,
            "best_per_target_pruning": False,
            "cost_top_k_pruning": False,
            "alternative_routes_preserved": True,
            "zero_edge_source_equals_target_is_temporal_continuation": True,
            "zero_edge_continuation_is_not_cycle": True,
            "condition_pool_validation": True,
        },
        "statistics": stats,
    }
    out["complete_path_count"] = len(cleaned_paths)
    out["paths"] = cleaned_paths
    return out, audits, stats


# =============================================================================
# Temporal outputs
# =============================================================================

def temporal_base(path: str) -> str:
    stem, _ = os.path.splitext(path)
    for suffix in ("_STRING_v12_temporal_results", "_temporal_results"):
        if stem.endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    return stem


def temporal_output_paths(input_path: str) -> Dict[str, str]:
    b = temporal_base(input_path)
    return {
        "json": b + "_STAGE4_clean_temporal.json",
        "clean": b + "_STAGE4_clean_results.txt",
        "topology": b + "_STAGE4_topology.txt",
        "topology_results": b + "_STAGE4_topology_results.txt",
        "audit": b + "_STAGE4_audit.csv",
        "nodes": b + "_STAGE4_nodes.csv",
        "edges": b + "_STAGE4_edges.csv",
        "validation_input": b + "_STAGE4_validation_input.txt",
        "validation_manifest": b + "_STAGE4_validation_manifest.json",
    }


def write_temporal_text(path: str, payload: Dict[str, Any], title: str):
    paths = payload.get("paths") or []
    s4 = payload.get("stage4") or {}
    stats = s4.get("statistics") or {}
    out = []
    out.append(f"*** {title} ***")
    out.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    out.append(f"Stage 4 version: {VERSION}")
    out.append("")
    out.append("Stage 4 policy:")
    out.append("  - Preserve biologically distinct alternative routes")
    out.append("  - Do NOT truncate a path because a gene reappears in a later interval")
    out.append("  - Temporal identity uses (gene, interval, role)")
    out.append("  - Remove exact temporal duplicates only")
    out.append("  - Preserve zero-edge source==target temporal continuations as anchor carry-forward")
    out.append("  - Reject genuine within-segment PPI cycles and broken interval continuity")
    out.append("  - No best-per-target or top-k post-processing pruning")
    out.append("")
    for k, v in stats.items():
        out.append(f"{k}: {v}")
    out.append("")
    out.append("*** Clean Temporal Path Results ***")
    for rec in paths:
        pid = rec.get("path_id")
        source_id = rec.get("source_path_id")
        genes = [norm_gene(x) for x in (rec.get("genes") or [])]
        cost = path_cost(rec)
        suffix = f" [source #{source_id}]" if source_id not in (None, "") else ""
        out.append(f"Path #{pid}{suffix}: {' -> '.join(genes)}  | total_cost={cost:.6f}")
        for seg in rec.get("segments", []) or []:
            ch = segment_chain(seg)
            role = seg.get("target_role", "")
            toks = seg.get("target_anchor_tokens") or []
            token_text = ", ".join(str(x) for x in toks) if toks else "-"
            out.append(
                f"    {seg.get('interval','?')}: {' -> '.join(ch)}"
                f" | D_used={len(seg.get('bridge_nodes') or [])}"
                f" | target_role={role} | target_sites={token_text}"
            )
        out.append("")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def write_temporal_topology(path: str, payload: Dict[str, Any]):
    groups: Dict[Tuple, List[Dict[str, Any]]] = defaultdict(list)
    for rec in payload.get("paths", []) or []:
        groups[topology_signature(rec)].append(rec)

    out = []
    out.append("*** Stage 4 Unique Interval-Aware Topologies ***")
    out.append("Distinct bridge routes remain distinct topologies; site tokens/costs are annotations.")
    out.append(f"Unique topologies: {len(groups)}")
    out.append("")

    for i, (sig, recs) in enumerate(sorted(groups.items(), key=lambda kv: min(path_cost(r) for r in kv[1])), 1):
        exemplar = min(recs, key=path_cost)
        out.append(
            f"Topology #{i}: {' -> '.join(exemplar.get('genes') or [])}"
            f" | supporting records={len(recs)} | best_cost={path_cost(exemplar):.6f}"
        )
        for seg in exemplar.get("segments", []) or []:
            out.append(f"    {seg.get('interval','?')}: {' -> '.join(segment_chain(seg))}")
        out.append("")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def write_drawing_topology_results(path: str, payload: Dict[str, Any]):
    # Separate flattened topology export for drawing/Stage-5 use.
    # Exact flattened gene-chain duplicates are grouped, but paths are never
    # truncated or cycle-cleaned here.
    records = list(payload.get("paths", []) or [])

    groups: Dict[Tuple[str, ...], List[Dict[str, Any]]] = defaultdict(list)
    for rec in records:
        sig = drawing_topology_signature(rec)
        if sig:
            groups[sig].append(rec)

    rows = []
    for sig, recs in groups.items():
        exemplar = min(
            recs,
            key=lambda r: (
                path_cost(r),
                str(r.get("source_path_id", r.get("path_id", ""))),
                str(r.get("path_id", ""))
            )
        )
        rows.append((sig, recs, exemplar))

    rows.sort(
        key=lambda x: (
            path_cost(x[2]),
            "->".join(x[0])
        )
    )

    protein_freq = Counter()
    first_hop_freq = Counter()

    for rec in records:
        genes = list(drawing_topology_signature(rec))
        # Frequency = number of Stage-4-clean paths containing the gene.
        for gene in set(genes):
            protein_freq[gene] += 1
        if len(genes) > 1:
            first_hop_freq[genes[1]] += 1

    path_den = max(1, len(records))
    hop_den = max(1, sum(first_hop_freq.values()))

    out = []
    out.append("*** Stage 4 Drawing Topology Results ***")
    out.append("")
    out.append(
        "*** Unique Topology Paths "
        "(Stage-4 validated; flattened gene-chain view) ***"
    )
    out.append(
        f"(Validated clean paths {len(records)} | "
        f"Unique topologies {len(rows)})"
    )
    out.append(
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )
    out.append("")
    out.append(
        "NOTE: Drawing export only. Stage-4 validation remains "
        "interval-aware; cross-interval recurrence is not cycle-truncated."
    )
    out.append("")

    out.append("Top 10 protein frequencies:")
    if protein_freq:
        for gene, count in protein_freq.most_common(10):
            out.append(
                f"  {gene}: {count} "
                f"({100.0 * count / path_den:.1f}%)"
            )
    else:
        out.append("  (none)")
    out.append("")

    out.append("Top 1st-hop frequencies:")
    if first_hop_freq:
        for gene, count in first_hop_freq.most_common(10):
            out.append(
                f"  {gene}: {count} "
                f"({100.0 * count / hop_den:.1f}%)"
            )
    else:
        out.append("  (none)")
    out.append("")

    for i, (sig, recs, exemplar) in enumerate(rows, 1):
        original_id = exemplar.get(
            "source_path_id",
            exemplar.get("path_id", "")
        )
        out.append(
            f"Unique Path #{i} (Original #{original_id}): "
            f"{' -> '.join(sig)}"
        )

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def write_audit_csv(path: str, audits: Sequence[AuditResult]):
    fields = [
        "path_id", "status", "valid", "total_cost", "reasons", "warnings",
        "original_genes", "reconstructed_genes"
    ]
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for a in audits:
            w.writerow({
                "path_id": a.path_id,
                "status": a.status,
                "valid": a.valid,
                "total_cost": a.total_cost,
                "reasons": " | ".join(a.reasons),
                "warnings": " | ".join(a.warnings),
                "original_genes": " -> ".join(a.original_genes),
                "reconstructed_genes": " -> ".join(a.reconstructed_genes),
            })


def write_graph_csvs(node_path: str, edge_path: str, payload: Dict[str, Any]):
    node_rows = []
    edge_rows = []
    for rec in payload.get("paths", []) or []:
        pid = rec.get("path_id")
        for seg in rec.get("segments", []) or []:
            interval = str(seg.get("interval", ""))
            src = norm_gene(seg.get("source_anchor"))
            tgt = norm_gene(seg.get("target_anchor"))
            bridges = [norm_gene(x) for x in (seg.get("bridge_nodes") or []) if norm_gene(x)]
            tgt_role = str(seg.get("target_role", "measured"))
            target_tokens = ";".join(str(x) for x in (seg.get("target_anchor_tokens") or []))

            if src:
                node_rows.append({"path_id": pid, "interval": interval, "gene": src,
                                  "role": "source_anchor", "target_role": "", "site_tokens": ""})
            for b in bridges:
                node_rows.append({"path_id": pid, "interval": interval, "gene": b,
                                  "role": "prior_bridge", "target_role": "", "site_tokens": ""})
            if tgt:
                node_rows.append({"path_id": pid, "interval": interval, "gene": tgt,
                                  "role": "target_anchor", "target_role": tgt_role,
                                  "site_tokens": target_tokens})

            chain = ([src] if src else []) + bridges + ([tgt] if tgt else [])
            for ei in range(len(chain) - 1):
                edge_rows.append({
                    "path_id": pid,
                    "interval": interval,
                    "edge_index": ei + 1,
                    "source": chain[ei],
                    "target": chain[ei+1],
                    "bridge_depth": len(bridges),
                    "target_role": tgt_role,
                })

    with open(node_path, "w", newline="", encoding="utf-8-sig") as f:
        fields = ["path_id", "interval", "gene", "role", "target_role", "site_tokens"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(node_rows)
    with open(edge_path, "w", newline="", encoding="utf-8-sig") as f:
        fields = ["path_id", "interval", "edge_index", "source", "target", "bridge_depth", "target_role"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(edge_rows)


# =============================================================================
# Text / legacy parser
# =============================================================================
SECTION_MARKERS = (
    "Complete Path Results", "Complete Paths", "Full Path Results", "완전 경로 결과"
)
PATH_PATTERNS = [
    re.compile(r"^(?:Path|경로)\s*#\s*(\d+)\s*:\s*(.+)$", re.I),
    re.compile(r"^(?:Path|경로)\s*(\d+)\s*:\s*(.+)$", re.I),
    re.compile(r"^#\s*(\d+)\s*:\s*(.+)$", re.I),
    re.compile(r"^(\d+)\s*[:.)-]\s*(.+)$", re.I),
]


def parse_path_line(line: str):
    s = line.strip()
    for pat in PATH_PATTERNS:
        m = pat.match(s)
        if m:
            nodes = [x.strip() for x in m.group(2).split("->") if x.strip()]
            # Remove optional annotations after | on final node.
            if nodes and "|" in nodes[-1]:
                nodes[-1] = nodes[-1].split("|", 1)[0].strip()
            return int(m.group(1)), nodes
    return None


def parse_text_file(path: str):
    with open(path, encoding="utf-8-sig") as f:
        lines = f.readlines()
    temporal = any("Input Schema" in x and "Temporal v2" in x for x in lines)
    paths = []
    header = []
    started = False
    for line in lines:
        if any(marker.lower() in line.lower() for marker in SECTION_MARKERS):
            started = True
        if not started:
            header.append(line.rstrip("\r\n"))
        p = parse_path_line(line)
        if p:
            paths.append(p)
    return temporal, header, paths


def find_companion_temporal_json(txt_path: str) -> Optional[str]:
    d = os.path.dirname(txt_path)
    name = os.path.basename(txt_path)
    candidates = []
    if "_STRING_v12_results" in name:
        candidates.append(name.replace("_STRING_v12_results", "_STRING_v12_temporal_results").rsplit(".",1)[0] + ".json")
    stem = os.path.splitext(name)[0]
    # Tolerate suffixes like (1) in copied files by searching same prefix.
    prefix = stem.split("_STRING_v12_results", 1)[0]
    try:
        for fn in os.listdir(d):
            if fn.startswith(prefix) and fn.endswith("_STRING_v12_temporal_results.json"):
                candidates.append(fn)
    except Exception:
        pass
    for fn in uniq_preserve(candidates):
        p = os.path.join(d, fn)
        if os.path.isfile(p):
            return p
    return None


def remove_consecutive_duplicates(nodes: List[str]) -> List[str]:
    out = []
    for x in nodes:
        if not out or x != out[-1]:
            out.append(x)
    return out


def legacy_remove_cycle(nodes: List[str]) -> List[str]:
    out, seen = [], set()
    for x in nodes:
        if x in seen:
            break
        out.append(x); seen.add(x)
    return out


def clean_text_mode(path: str, legacy_cycle_cleanup: bool):
    temporal, header, raw = parse_text_file(path)
    audits = []
    cleaned = []
    seen = set()

    for idx, nodes in raw:
        orig = list(nodes)
        # Temporal flat text lacks interval identity. Never infer cycles from gene recurrence.
        if temporal:
            new = orig
            mode = "temporal_text_exact_dedup_only"
        else:
            new = remove_consecutive_duplicates(orig)
            if legacy_cycle_cleanup:
                new = legacy_remove_cycle(new)
            mode = "legacy_cleanup"

        key = tuple(new)
        status = "kept"
        warnings = []
        if key in seen:
            status = "exact_duplicate_removed"
        else:
            seen.add(key)
            cleaned.append((idx, new))
        audits.append({
            "path_id": idx,
            "status": status,
            "valid": True,
            "total_cost": "",
            "reasons": "",
            "warnings": ("Temporal text has no interval metadata; no recurrence/cycle truncation was attempted"
                         if temporal else mode),
            "original_genes": " -> ".join(orig),
            "reconstructed_genes": " -> ".join(new),
        })

    stem, _ = os.path.splitext(path)
    if stem.endswith("_results"):
        stem = stem[:-len("_results")]
    outs = {
        "clean": stem + "_STAGE4_clean_results.txt",
        "topology": stem + "_STAGE4_topology.txt",
        "topology_results": stem + "_STAGE4_topology_results.txt",
        "audit": stem + "_STAGE4_audit.csv",
    }

    out_lines = []
    out_lines.extend(header)
    out_lines.append("")
    out_lines.append("*** Stage 4 Clean Path Results ***")
    out_lines.append(f"Input mode: {'Temporal v2 flat-text fallback' if temporal else 'Legacy'}")
    out_lines.append(f"Input paths: {len(raw)} | Retained unique paths: {len(cleaned)}")
    if temporal:
        out_lines.append("Policy: exact duplicate chains only; repeated genes are NOT interpreted as cycles without interval metadata.")
    else:
        out_lines.append(f"Legacy cycle cleanup: {legacy_cycle_cleanup}")
    out_lines.append("")
    for n, (orig_id, genes) in enumerate(cleaned, 1):
        out_lines.append(f"Path #{n} [source #{orig_id}]: {' -> '.join(genes)}")
    with open(outs["clean"], "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines))

    with open(outs["topology"], "w", encoding="utf-8") as f:
        f.write("*** Stage 4 Unique Gene-Chain Topologies ***\n")
        f.write(f"Unique paths: {len(cleaned)}\n\n")
        for n, (orig_id, genes) in enumerate(cleaned, 1):
            f.write(f"Topology #{n} [source #{orig_id}]: {' -> '.join(genes)}\n")

    protein_freq = Counter()
    first_hop_freq = Counter()
    for _, genes in cleaned:
        for gene in set(genes):
            protein_freq[gene] += 1
        if len(genes) > 1:
            first_hop_freq[genes[1]] += 1

    with open(outs["topology_results"], "w", encoding="utf-8") as f:
        f.write("*** Stage 4 Drawing Topology Results ***\n\n")
        f.write("*** Unique Topology Paths ***\n")
        f.write(
            f"(Retained clean paths {len(cleaned)} | "
            f"Unique topologies {len(cleaned)})\n"
        )
        f.write(
            f"Generated: "
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        )

        f.write("Top 10 protein frequencies:\n")
        denom = max(1, len(cleaned))
        for gene, count in protein_freq.most_common(10):
            f.write(
                f"  {gene}: {count} "
                f"({100.0 * count / denom:.1f}%)\n"
            )

        f.write("\nTop 1st-hop frequencies:\n")
        hop_denom = max(1, sum(first_hop_freq.values()))
        for gene, count in first_hop_freq.most_common(10):
            f.write(
                f"  {gene}: {count} "
                f"({100.0 * count / hop_denom:.1f}%)\n"
            )
        f.write("\n")

        for n, (orig_id, genes) in enumerate(cleaned, 1):
            f.write(
                f"Unique Path #{n} (Original #{orig_id}): "
                f"{' -> '.join(genes)}\n"
            )

    with open(outs["audit"], "w", newline="", encoding="utf-8-sig") as f:
        fields = ["path_id", "status", "valid", "total_cost", "reasons", "warnings", "original_genes", "reconstructed_genes"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(audits)

    stats = {"input_path_count": len(raw), "clean_path_count": len(cleaned), "temporal_text_fallback": temporal}
    return outs, stats



def _collect_validation_site_tokens(payload: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Collect observed temporal-anchor phosphosites as (token, gene).

    Only experimentally observed target-anchor tokens are exported. Prior-network
    bridge proteins are intentionally not assigned synthetic phosphosites.
    """
    seen = set()
    rows = []
    for rec in payload.get("paths", []) or []:
        for seg in rec.get("segments", []) or []:
            for tok in seg.get("target_anchor_tokens", []) or []:
                s = str(tok).strip().lstrip("-")
                if "@" not in s:
                    continue
                gene = norm_gene(s.split("@", 1)[0])
                if not gene:
                    continue
                key = (s, gene)
                if key not in seen:
                    seen.add(key)
                    rows.append(key)
    rows.sort(key=lambda x: (x[1], x[0]))
    return rows


def write_validation_input(path: str, payload: Dict[str, Any]):
    """Write one flat-text file accepted by both supplied validators.

    It contains a Gene@Site mapping table plus Stage-4-clean pathway lines.
    Temporal metadata itself stays authoritative in the companion manifest/JSON.
    """
    site_rows = _collect_validation_site_tokens(payload)
    paths = payload.get("paths", []) or []
    out = []
    out.append("*** Stage 4 Validation Input ***")
    out.append("Purpose: shared input for PPI Cell-Type Validator and PPI Phosphosite Validator")
    out.append("Only Stage-4-retained temporal routes are included.")
    out.append("Observed phosphosites are temporal target-anchor measurements; bridge proteins are not assigned synthetic sites.")
    out.append("")
    out.append("*** Protein Name to Gene Mapping Table ***")
    out.append("| Original Name | Gene Symbol |")
    out.append("|---------------|-------------|")
    for tok, gene in site_rows:
        out.append(f"| {tok} | {gene} |")
    out.append("")
    out.append("*** Complete Path Results ***")
    out.append(f"Total {len(paths)} Stage-4-clean paths.")
    out.append("")
    for i, rec in enumerate(paths, 1):
        pid = rec.get("path_id", i)
        genes = [norm_gene(x) for x in (rec.get("genes") or []) if norm_gene(x)]
        out.append(f"Path #{pid}: {' -> '.join(genes)}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def write_validation_manifest(path: str, payload: Dict[str, Any]):
    """Machine-readable bridge between Stage 4 validators and Stage 5 graphics."""
    pool = payload.get("condition_functional_pool") or {}
    manifest_paths = []
    for rec in payload.get("paths", []) or []:
        item = {
            "path_id": rec.get("path_id"),
            "source_path_id": rec.get("source_path_id"),
            "genes": [norm_gene(x) for x in (rec.get("genes") or [])],
            "total_cost": path_cost(rec),
            "segments": [],
        }
        for seg in rec.get("segments", []) or []:
            item["segments"].append({
                "interval": seg.get("interval"),
                "source_anchor": norm_gene(seg.get("source_anchor")),
                "bridge_nodes": [norm_gene(x) for x in (seg.get("bridge_nodes") or [])],
                "target_anchor": norm_gene(seg.get("target_anchor")),
                "target_role": seg.get("target_role", ""),
                "target_anchor_tokens": list(seg.get("target_anchor_tokens") or []),
                "bridge_depth": len(seg.get("bridge_nodes") or []),
            })
        manifest_paths.append(item)
    doc = {
        "schema_version": "4.1",
        "purpose": "validation_and_stage5_bridge",
        "stage4_version": VERSION,
        "validation_policy": {
            "cleaning_does_not_use_celltype_grade": True,
            "cleaning_does_not_use_phosphosite_grade": True,
            "validators_annotate_but_do_not_delete_routes": True,
            "recommended_join_key_edge": ["source", "target"],
            "recommended_join_key_path": ["path_id"],
        },
        "condition_functional_pool": pool,
        "paths": manifest_paths,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=2)


# =============================================================================
# GUI
# =============================================================================
class Stage4CleanerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("BFS Path Cleaner — Stage 4.3 Validation-Ready")
        self.geometry("940x720")
        self.minsize(840, 620)
        self.configure(bg="#f3f6f9")
        self.input_path: Optional[str] = None
        self._build_ui()

    def _build_ui(self):
        head = tk.Frame(self, bg="#102a43", padx=16, pady=13)
        head.pack(fill="x")
        tk.Label(head, text="BFS Path Cleaner — Stage 4.3",
                 bg="#102a43", fg="white", font=("Arial", 17, "bold")).pack(anchor="w")
        tk.Label(head,
                 text="Temporal-aware cleaning • preserve alternative routes • validator-ready + topology-drawing export",
                 bg="#102a43", fg="#bcccdc", font=("Arial", 10)).pack(anchor="w")

        body = tk.Frame(self, bg="#f3f6f9", padx=18, pady=14)
        body.pack(fill="both", expand=True)

        f1 = ttk.LabelFrame(body, text="1) Input", padding=10)
        f1.pack(fill="x", pady=5)
        self.file_var = tk.StringVar(value="Select Stage 3.2 temporal JSON or *_results.txt...")
        ttk.Entry(f1, textvariable=self.file_var).pack(side="left", fill="x", expand=True, padx=(0,8))
        ttk.Button(f1, text="Browse", command=self.choose_file).pack(side="left")

        info = ttk.LabelFrame(body, text="Stage 4 safety policy", padding=10)
        info.pack(fill="x", pady=5)
        policy = (
            "Temporal v2: uses structured interval segments when available. A gene may legitimately reappear in a later interval; "
            "such recurrence is NOT removed. Zero-edge source==target segments are treated as temporal anchor continuation, not PPI self-loops. "
            "Only exact temporal duplicates and genuinely structurally invalid records are removed."
        )
        ttk.Label(info, text=policy, wraplength=850, justify="left").pack(anchor="w")
        ttk.Label(info, text="Preferred input: *_STRING_v12_temporal_results.json",
                  foreground="#0b7285").pack(anchor="w", pady=(4,0))

        f2 = ttk.LabelFrame(body, text="2) Legacy fallback options", padding=10)
        f2.pack(fill="x", pady=5)
        self.legacy_cycle_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            f2,
            text="For Legacy v1 text only: remove repeated-gene cycles using the old rule",
            variable=self.legacy_cycle_var,
        ).pack(anchor="w")
        ttk.Label(
            f2,
            text="This option is ignored for Temporal v2 files.",
            foreground="#6b7280",
        ).pack(anchor="w")

        f3 = ttk.LabelFrame(body, text="3) Output", padding=10)
        f3.pack(fill="x", pady=5)
        self.output_var = tk.StringVar(value="Outputs will be generated beside the selected file.")
        ttk.Label(f3, textvariable=self.output_var, wraplength=850, justify="left").pack(anchor="w")

        controls = tk.Frame(body, bg="#f3f6f9")
        controls.pack(fill="x", pady=8)
        self.run_btn = ttk.Button(controls, text="Run Stage 4", command=self.run)
        self.run_btn.pack(side="left")
        ttk.Button(controls, text="Open Output Folder", command=self.open_folder).pack(side="left", padx=8)
        ttk.Button(controls, text="Reset", command=self.reset).pack(side="left")

        self.progress = ttk.Progressbar(body, mode="determinate", maximum=100)
        self.progress.pack(fill="x", pady=(2,4))
        self.status_var = tk.StringVar(value="Waiting...")
        ttk.Label(body, textvariable=self.status_var).pack(anchor="w")

        lf = ttk.LabelFrame(body, text="Execution Log", padding=6)
        lf.pack(fill="both", expand=True, pady=5)
        self.log_box = scrolledtext.ScrolledText(lf, height=16, font=("Consolas", 9), wrap="word")
        self.log_box.pack(fill="both", expand=True)
        self.log_box.configure(state="disabled")

    def log(self, msg: str):
        def _do():
            self.log_box.configure(state="normal")
            self.log_box.insert("end", msg.rstrip() + "\n")
            self.log_box.see("end")
            self.log_box.configure(state="disabled")
        self.after(0, _do)

    def set_status(self, text: str, progress: Optional[float] = None):
        self.after(0, lambda: self.status_var.set(text))
        if progress is not None:
            self.after(0, lambda: self.progress.configure(value=progress))

    def choose_file(self):
        p = filedialog.askopenfilename(
            title="Select Stage 3.2 temporal JSON or results text",
            filetypes=[
                ("Stage 3.2 outputs", "*.json *.txt"),
                ("JSON files", "*.json"),
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
        )
        if not p:
            return
        self.input_path = p
        self.file_var.set(p)
        self.output_var.set("Outputs will be generated beside: " + os.path.dirname(p))
        self.log(f"Selected: {p}")

    def run(self):
        if not self.input_path or not os.path.isfile(self.input_path):
            messagebox.showerror("Input Error", "Please select a valid Stage 3.2 output file.")
            return
        self.run_btn.configure(state="disabled")
        self.progress.configure(value=0)
        threading.Thread(target=self._process, daemon=True).start()

    def _process(self):
        try:
            p = self.input_path
            assert p is not None
            self.log("=" * 70)
            self.log(f"Stage 4 started: {os.path.basename(p)}")
            self.set_status("Detecting input type...", 5)

            actual = p
            structured = p.lower().endswith(".json")
            if p.lower().endswith(".txt"):
                companion = find_companion_temporal_json(p)
                if companion:
                    self.log(f"Companion temporal JSON detected: {os.path.basename(companion)}")
                    self.log("Using structured JSON to avoid unsafe cycle inference from flat text.")
                    actual = companion
                    structured = True

            if structured:
                self.set_status("Reading structured temporal JSON...", 15)
                with open(actual, encoding="utf-8-sig") as f:
                    payload = json.load(f)
                if not isinstance(payload, dict) or not isinstance(payload.get("paths"), list):
                    raise ValueError("Selected JSON is not a Stage 3 temporal-results JSON (missing 'paths').")

                self.log(f"Input temporal paths: {len(payload.get('paths') or [])}")
                self.log(f"Bridge Depth D: {payload.get('bridge_depth_D', 'unknown')}")
                pool = payload.get("condition_functional_pool") or {}
                self.log(f"Condition-pool excluded genes: {len(pool.get('excluded_low_activity_genes') or [])}")
                self.set_status("Validating temporal structure...", 40)

                cleaned, audits, stats = clean_temporal_payload(payload)
                outs = temporal_output_paths(actual)
                self.log("Validation summary:")
                for k, v in stats.items():
                    self.log(f"  {k}: {v}")

                self.set_status("Writing Stage 4 outputs...", 75)
                with open(outs["json"], "w", encoding="utf-8") as f:
                    json.dump(cleaned, f, ensure_ascii=False, indent=2)
                write_temporal_text(outs["clean"], cleaned, "Stage 4 Temporal-Aware Clean Results")
                write_temporal_topology(outs["topology"], cleaned)
                write_drawing_topology_results(
                    outs["topology_results"], cleaned
                )
                write_audit_csv(outs["audit"], audits)
                write_graph_csvs(outs["nodes"], outs["edges"], cleaned)
                write_validation_input(outs["validation_input"], cleaned)
                write_validation_manifest(outs["validation_manifest"], cleaned)

                self.log("Validation-ready exports generated:")
                self.log(
                    "  - Drawing topology paths: "
                    + os.path.basename(outs["topology_results"])
                )
                self.log("  - Shared validator input: " + os.path.basename(outs["validation_input"]))
                self.log("  - Temporal/path join manifest: " + os.path.basename(outs["validation_manifest"]))
                self.log("  Cell-type and phosphosite grades are annotations only; Stage 4 does not delete routes based on them.")

                self.log("Outputs:")
                for k, fp in outs.items():
                    self.log(f"  {k}: {fp}")
                self.output_var.set(" | ".join(os.path.basename(x) for x in outs.values()))
            else:
                self.set_status("Processing text fallback mode...", 35)
                temporal, _, _ = parse_text_file(actual)
                if temporal:
                    self.log("Temporal v2 text detected but no companion structured JSON was found.")
                    self.log("SAFE FALLBACK: exact duplicate gene chains only; no cycle truncation.")
                else:
                    self.log("Legacy text mode detected.")
                outs, stats = clean_text_mode(actual, self.legacy_cycle_var.get())
                for k, v in stats.items():
                    self.log(f"  {k}: {v}")
                self.log("Outputs:")
                for k, fp in outs.items():
                    self.log(f"  {k}: {fp}")
                self.output_var.set(" | ".join(os.path.basename(x) for x in outs.values()))

            self.set_status("Completed.", 100)
            self.log("Stage 4 completed successfully.")
            self.log("=" * 70)
            self.after(0, lambda: messagebox.showinfo(
                "Stage 4 Completed",
                "Temporal-aware path cleaning and validation completed successfully.\n\n"
                "Alternative routes and zero-edge temporal continuations were preserved; no top-k or best-target pruning was applied."
            ))
        except Exception as e:
            tb = traceback.format_exc()
            self.log(f"ERROR: {e}")
            self.log(tb)
            self.set_status("Error.", 0)
            msg = str(e)
            self.after(0, lambda m=msg: messagebox.showerror("Stage 4 Error", m))
        finally:
            self.after(0, lambda: self.run_btn.configure(state="normal"))

    def open_folder(self):
        if not self.input_path:
            messagebox.showinfo("Info", "Please select an input file first.")
            return
        folder = os.path.dirname(self.input_path)
        try:
            if os.name == "nt":
                os.startfile(folder)  # type: ignore[attr-defined]
            elif os.uname().sysname == "Darwin":
                import subprocess; subprocess.Popen(["open", folder])
            else:
                import subprocess; subprocess.Popen(["xdg-open", folder])
        except Exception as e:
            messagebox.showerror("Open Folder Error", str(e))

    def reset(self):
        self.input_path = None
        self.file_var.set("Select Stage 3.2 temporal JSON or *_results.txt...")
        self.output_var.set("Outputs will be generated beside the selected file.")
        self.progress.configure(value=0)
        self.status_var.set("Waiting...")
        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", "end")
        self.log_box.configure(state="disabled")


if __name__ == "__main__":
    app = Stage4CleanerApp()
    app.mainloop()
