#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reference transcript anchoring with adjacency-aware matching.

Match levels:
  3 = FULL             (SE: up-skip-dn adjacent in order; RI: up-dn adjacent with a cleanly retained intron)
  2 = PARTIAL          (regions match exactly but adjacency is not satisfied, etc.)
  1 = ANCHOR_ONLY      (only a single anchor region matches)
  0 = NONE             (no match)
"""

import logging
from collections import defaultdict

from .events import _extract_event_regions
from .nmd import transcript_self_nmd

logger = logging.getLogger(__name__)

# Transcript biotypes that the annotation itself marks as NMD substrates.
# BioSpliceView assumes the anchor (reference) transcript is normally expressed,
# so such transcripts must not be preferred as anchors.
NMD_BIOTYPES = {'nonsense_mediated_decay'}


def _check_boundary_match(query_region, tx_exons, mode='exact'):
    qs, qe = query_region
    for es, ee in tx_exons:
        if mode == 'exact' and es == qs and ee == qe:
            return True, 'EXACT'
        if mode == 'start' and es == qs:
            return True, 'START_MATCH'
        if mode == 'end' and ee == qe:
            return True, 'END_MATCH'
        if mode == 'contains' and es <= qs and ee >= qe:
            return True, 'CONTAINED'
    return False, 'NONE'


def _exon_in_ref(exon_region, ref_exons):
    """Return True if the genomic exon interval appears (contained) in the reference transcript's exon list."""
    es_, ee_ = exon_region
    for res, ree in ref_exons:
        if res <= es_ and ree >= ee_:
            return True
    return False


# =============================================================================
# ✅ Adjacency-check helpers
# =============================================================================

def _find_exon_index(region, tx_exons, tol=0):
    """Find the index of the exon in tx_exons exactly (or within tolerance) matching region."""
    rs, re_ = region
    for i, (es, ee) in enumerate(tx_exons):
        if abs(es - rs) <= tol and abs(ee - re_) <= tol:
            return i
    return -1


def _are_adjacent_in_tx(idx_a, idx_b, tx_exons):
    """
    Return True if the exons at idx_a and idx_b in tx_exons are adjacent
    (no exon in between), i.e. the absolute index difference is 1.
    """
    if idx_a < 0 or idx_b < 0:
        return False
    return abs(idx_a - idx_b) == 1


def _check_se_adjacency(tx_exons, up_region, skip_region, dn_region):
    """
    SE adjacency check: the upstream, skipped and downstream exons must be
    adjacent in order within the transcript, i.e.
    upstream_exon -> skipped_exon -> downstream_exon are consecutive.

    Returns (is_adjacent: bool, detail: str)
    """
    up_idx = _find_exon_index(up_region, tx_exons)
    skip_idx = _find_exon_index(skip_region, tx_exons)
    dn_idx = _find_exon_index(dn_region, tx_exons)

    if up_idx < 0 or skip_idx < 0 or dn_idx < 0:
        return False, f"missing_exon(up={up_idx},skip={skip_idx},dn={dn_idx})"

    # Must be strictly adjacent in order: skip between up and dn, indices consecutive
    if skip_idx == up_idx + 1 and dn_idx == skip_idx + 1:
        return True, "up-skip-dn_adjacent"

    return False, f"not_adjacent(up={up_idx},skip={skip_idx},dn={dn_idx})"


def _check_ri_adjacency(tx_exons, up_region, dn_region):
    """
    RI adjacency check: the upstream and downstream exons must be adjacent
    within the transcript.

    Returns (is_adjacent: bool, detail: str)
    """
    up_idx = _find_exon_index(up_region, tx_exons)
    dn_idx = _find_exon_index(dn_region, tx_exons)

    if up_idx < 0 or dn_idx < 0:
        return False, f"missing_exon(up={up_idx},dn={dn_idx})"

    if _are_adjacent_in_tx(up_idx, dn_idx, tx_exons):
        return True, "up-dn_adjacent"

    return False, f"not_adjacent(up={up_idx},dn={dn_idx})"


def _check_mxe_adjacency(tx_exons, up_region, dn_region, e1_region, e2_region):
    """
    MXE adjacency check:
      1. both upstream and downstream exons must be present
      2. a single transcript must not contain both mutually exclusive exons
      3. adjacent layouts are preferred (up-dn adjacent, or up-exon-dn consecutive)

    Returns (priority: int, detail: str); higher priority wins:
        4 = perfect: up-e1-dn or up-e2-dn consecutive, only one exclusive exon present
        3 = good: up and dn adjacent, only one exclusive exon present
        2 = fair: up and dn both present, only one exclusive exon, but not adjacent
        1 = poor: up and dn present, both or neither exclusive exon present
        0 = disqualified
    """
    up_idx = _find_exon_index(up_region, tx_exons)
    dn_idx = _find_exon_index(dn_region, tx_exons)
    e1_idx = _find_exon_index(e1_region, tx_exons)
    e2_idx = _find_exon_index(e2_region, tx_exons)

    # Both flanks must be present
    if up_idx < 0 or dn_idx < 0:
        return 0, f"missing_anchor(up={up_idx},dn={dn_idx})"

    has_e1 = e1_idx >= 0
    has_e2 = e2_idx >= 0

    # Mutual exclusivity: a transcript must not contain both
    if has_e1 and has_e2:
        return 1, "both_exons_present"

    # Exactly one exclusive exon present (or none)
    only_one = has_e1 != has_e2

    # Adjacency check
    up_dn_adjacent = _are_adjacent_in_tx(up_idx, dn_idx, tx_exons)

    # up - exon - dn consecutive check
    if has_e1:
        up_e1_dn_adjacent = (e1_idx == up_idx + 1 and dn_idx == e1_idx + 1)
    elif has_e2:
        up_e2_dn_adjacent = (e2_idx == up_idx + 1 and dn_idx == e2_idx + 1)
        up_e1_dn_adjacent = False
    else:
        up_e1_dn_adjacent = False

    if only_one and up_e1_dn_adjacent:
        which = "e1" if has_e1 else "e2"
        return 4, f"perfect(up-{which}-dn_adjacent)"

    if only_one and up_dn_adjacent:
        which = "e1" if has_e1 else "e2"
        return 3, f"good(up-dn_adjacent;{which}_only)"

    if only_one:
        which = "e1" if has_e1 else "e2"
        return 2, f"ok(up,dn_present;{which}_only;not_adjacent)"

    # Neither exclusive exon present
    return 1, "neither_exon_present"


def _classify_match_v324(tx_exons, anchors, diff_regions, etype, row=None):
    """
    Transcript match classification (adjacency-aware version).

    Core rules:
      - SE: up - skipped - dn must be adjacent in order
      - RI: up - dn must be adjacent
      - MXE: both flanks present, only one exclusive exon, adjacency preferred
      - A5SS/A3SS: original logic (short/long exon + flanking match)
    """
    tx_set = set(tx_exons)
    tx_exons_sorted = sorted(tx_exons, key=lambda x: x[0])

    # ===================== SE =====================
    if etype == "SE":
        if len(anchors) < 2 or not diff_regions:
            return 0, "NONE", "missing_regions"

        up_region = anchors[0]
        dn_region = anchors[1]
        skip_region = diff_regions[0]

        up_exact = up_region in tx_set
        dn_exact = dn_region in tx_set
        skip_exact = skip_region in tx_set

        # Exact-match counting
        exact_count = sum([up_exact, dn_exact, skip_exact])

        if exact_count == 3:
            # All three match exactly; check adjacency
            is_adj, adj_detail = _check_se_adjacency(tx_exons_sorted, up_region, skip_region, dn_region)
            if is_adj:
                return 3, "FULL", f"all_3_exact+{adj_detail}"
            else:
                # All present but not consecutive -> downgrade
                return 2, "PARTIAL", f"all_3_exact_but_{adj_detail}"

        if exact_count == 2:
            # Check for up-skip or skip-dn adjacency (partial match)
            up_idx = _find_exon_index(up_region, tx_exons_sorted)
            skip_idx = _find_exon_index(skip_region, tx_exons_sorted)
            dn_idx = _find_exon_index(dn_region, tx_exons_sorted)

            if up_exact and skip_exact and _are_adjacent_in_tx(up_idx, skip_idx, tx_exons_sorted):
                # up-skip adjacent, dn missing
                if dn_idx >= 0 and dn_idx == skip_idx + 1:
                    return 3, "FULL", "up_skip_dn_adjacent(partial_exact)"
                return 2, "PARTIAL", "up_skip_adjacent;dn_missing"

            if skip_exact and dn_exact and _are_adjacent_in_tx(skip_idx, dn_idx, tx_exons_sorted):
                # skip-dn adjacent, up missing
                if up_idx >= 0 and skip_idx == up_idx + 1:
                    return 3, "FULL", "up_skip_dn_adjacent(partial_exact)"
                return 2, "PARTIAL", "skip_dn_adjacent;up_missing"

            if up_exact and dn_exact and _are_adjacent_in_tx(up_idx, dn_idx, tx_exons_sorted):
                # up-dn adjacent with the exon skipped (alt form of an SE event)
                return 2, "PARTIAL", "up_dn_adjacent(skip_skipped_form)"

            # Fallback: 2 exact matches but not adjacent
            matched = []
            if up_exact: matched.append("up")
            if skip_exact: matched.append("skip")
            if dn_exact: matched.append("dn")
            return 2, "PARTIAL", f"2_exact({'+'.join(matched)});not_adjacent"

        if exact_count == 1:
            return 1, "ANCHOR_ONLY", f"1_exact_only"

        return 0, "NONE", "no_exact_match"

    # ===================== RI =====================
    elif etype == "RI":
        if len(anchors) < 2:
            return 0, "NONE", "missing_anchors"

        up_region = anchors[0]
        dn_region = anchors[1]

        up_exact = up_region in tx_set
        dn_exact = dn_region in tx_set

        if up_exact and dn_exact:
            is_adj, adj_detail = _check_ri_adjacency(tx_exons_sorted, up_region, dn_region)
            if is_adj:
                # Check whether the intron is cleanly retained (not covered by another exon)
                ri_clean = True
                if diff_regions:
                    rs, re_ = diff_regions[0]
                    for es, ee in tx_exons_sorted:
                        if es <= rs and ee >= re_:
                            ri_clean = False
                            break
                if ri_clean:
                    return 3, "FULL", f"up_dn_adjacent+ri_clean"
                else:
                    return 2, "PARTIAL", f"up_dn_adjacent;ri_overlaps_exon"
            else:
                return 2, "PARTIAL", f"up_dn_exact_but_{adj_detail}"

        if up_exact or dn_exact:
            matched = "up" if up_exact else "dn"
            return 1, "ANCHOR_ONLY", f"{matched}_exact_only"

        return 0, "NONE", "no_anchor_match"

    # ===================== MXE =====================
    elif etype == "MXE":
        if len(diff_regions) < 2:
            return 0, "NONE", "no_diff_regions"

        e1_region = diff_regions[0]
        e2_region = diff_regions[1]

        if len(anchors) < 2:
            return 0, "NONE", "missing_anchors"

        up_region = anchors[0]
        dn_region = anchors[1]

        priority, detail = _check_mxe_adjacency(
            tx_exons_sorted, up_region, dn_region, e1_region, e2_region
        )

        if priority >= 4:
            return 3, "FULL", detail
        elif priority == 3:
            return 2, "PARTIAL", detail
        elif priority == 2:
            return 2, "PARTIAL", detail
        elif priority == 1:
            return 1, "ANCHOR_ONLY", detail
        else:
            return 0, "NONE", detail

    # ===================== A5SS / A3SS =====================
    else:
        if len(anchors) < 2:
            return 0, "NONE", "missing_anchors"

        short_exon = anchors[0]
        flank_exon = anchors[1]

        short_exact = short_exon in tx_set
        long_exact = False
        if row is not None:
            try:
                long_exon = (int(float(row['longExonStart_0base'])), int(float(row['longExonEnd'])))
                long_exact = long_exon in tx_set
            except (KeyError, ValueError, TypeError):
                long_exact = False

        flank_exact = flank_exon in tx_set
        flank_boundary_ok, flank_boundary_type = _check_boundary_match(flank_exon, tx_exons_sorted, 'contains')

        if (short_exact or long_exact) and flank_exact:
            iso = "short" if short_exact else "long"
            return 3, "FULL", f"{iso}_exon_exact+flank_exact"

        if (short_exact or long_exact) and flank_boundary_ok:
            iso = "short" if short_exact else "long"
            return 2, "PARTIAL", f"{iso}_exact+flank_boundary({flank_boundary_type})"

        if short_exact or long_exact:
            iso = "short" if short_exact else "long"
            return 1, "ANCHOR_ONLY", f"{iso}_exact_only"

        if flank_exact:
            return 1, "ANCHOR_ONLY", "flank_exact_only"

        if flank_boundary_ok:
            return 1, "ANCHOR_ONLY", f"flank_boundary({flank_boundary_type})"

        return 0, "NONE", "no_match"


def _anchor_safety(tid, cds_list, tx_exons, strand, tid_biotype):
    """Anchor safety tier of a coding candidate transcript.
    0 = safe (coding biotype, stop codon not NMD-sensitive by our own rules)
    1 = coding biotype but its annotated stop looks NMD-sensitive (e.g. long 3'UTR)
    2 = annotation itself marks the transcript as an NMD substrate
    """
    bt = (tid_biotype or {}).get(tid, '')
    if bt in NMD_BIOTYPES:
        return 2
    try:
        if transcript_self_nmd(cds_list, tx_exons, strand):
            return 1
    except Exception:  # never let the safety check break anchoring
        pass
    return 0


def anchor_ref_transcript(row, etype, gene_sym, gene_idx, gene_exons, tid_strand, alias_map,
                          longest_tx_per_gene, tid_biotype=None):
    """Anchor the event's reference transcript in the annotation.

    Candidates within the same match level are preferred in the order:
    (a) not NMD-annotated and not self-NMD, (b) not NMD-annotated but self-NMD,
    (c) NMD-annotated (kept only as a fallback and flagged in the status).
    Within each tier the longest total CDS wins, as before.

    Returns (transcript_id, cds_segments, exon_segments, strand, status_message).
    """
    sym = str(gene_sym).strip()
    canonical = alias_map.get(sym, sym)
    candidates = gene_idx.get(canonical, {})

    if not candidates:
        for alias, canon in alias_map.items():
            if alias == sym or alias == canonical:
                candidates = gene_idx.get(canon, {})
                if candidates:
                    canonical = canon
                    break

    if not candidates:
        for g in gene_idx.keys():
            if sym in g or g in sym:
                candidates = gene_idx.get(g, {})
                if candidates:
                    canonical = g
                    logger.debug(f"Found gene by partial match: {sym} -> {canonical}")
                    break

    if not candidates:
        has_exon = canonical in gene_exons and len(gene_exons[canonical]) > 0
        if has_exon:
            for tid, exons in gene_exons[canonical].items():
                if exons:
                    return tid, [], exons, '+', f"NON_CODING_RNA(gene={canonical})"
        return None, [], [], '+', f"UNANCHORED(GENE_NOT_FOUND;symbol={sym};canonical={canonical})"

    try:
        anchors, diff_regions = _extract_event_regions(row, etype)
    except KeyError as e:
        return None, [], [], '+', f"UNANCHORED(MISSING_COLUMN;{e})"
    if not anchors and not diff_regions:
        return None, [], [], '+', "UNANCHORED(NO_EVENT_REGIONS_PARSED)"

    buckets = defaultdict(list)
    noncoding_buckets = defaultdict(list)

    for tid in candidates:
        tx_exons = gene_exons.get(canonical, {}).get(tid, [])
        if not tx_exons:
            continue
        level_code, level_label, detail = _classify_match_v324(tx_exons, anchors, diff_regions, etype, row=row)
        if level_code == 0:
            continue
        cds_list = candidates[tid]
        is_coding = len(cds_list) > 0
        entry = (tid, cds_list, tx_exons, detail) if is_coding else (tid, detail)
        if is_coding:
            buckets[level_code].append(entry)
        else:
            noncoding_buckets[level_code].append(entry)

    for lv in [3, 2, 1]:
        if buckets[lv]:
            scored = [
                (_anchor_safety(tid, cds, exons, tid_strand.get(tid, '+'), tid_biotype), entry)
                for entry in buckets[lv]
                for tid, cds, exons, detail in [entry]
            ]
            best_tier = min(t for t, _ in scored)
            pool = [entry for t, entry in scored if t == best_tier]
            best = max(pool, key=lambda x: sum(e - s for s, e, _ in x[1]))
            tid, cds, exons, detail = best
            if lv == 3:
                status = "FULL"
            elif lv == 2:
                status = f"PARTIAL({detail})"
            else:
                status = f"ANCHOR_ONLY({detail})"
            if best_tier == 1:
                status += ";ANCHOR_SELF_NMD"
            elif best_tier == 2:
                status += ";ANCHOR_NMD_BIOTYPE"
            return tid, cds, exons, tid_strand.get(tid, '+'), status

    for lv in [3, 2, 1]:
        if noncoding_buckets[lv]:
            entries = noncoding_buckets[lv]
            tid = entries[0][0]
            tx_exons = gene_exons.get(canonical, {}).get(tid, [])
            detail = entries[0][1] if len(entries[0]) > 1 else ""
            lv_label = {3: "NON_CODING_FULL", 2: "NON_CODING_PARTIAL", 1: "NON_CODING_ANCHOR"}[lv]
            return tid, [], tx_exons, tid_strand.get(tid, '+'), f"{lv_label}({detail})"

    return None, [], [], '+', f"UNANCHORED(PRECISE_MATCH_FAIL;gene={canonical};n_tx={len(candidates)})"
