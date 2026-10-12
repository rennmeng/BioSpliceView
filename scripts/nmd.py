#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PTC genomic mapping and NMD sensitivity classification.

All exon/CDS intervals are 0-based half-open genomic coordinates.
"Last exon" and "penultimate exon" always refer to TRANSCRIPT order,
i.e. for '-' strand transcripts the last exon is the genomically lowest one.
"""


def _map_ptc_to_genomic(cds_ph_list, strand, ptc_nt_pos):
    if not cds_ph_list or ptc_nt_pos <= 0:
        return None

    if strand == '+':
        ordered_cds = sorted(cds_ph_list, key=lambda x: x[0])
    else:
        ordered_cds = sorted(cds_ph_list, key=lambda x: x[0], reverse=True)

    cum = 0
    for cs, ce, ph in ordered_cds:
        seg_len = ce - cs
        if cum + seg_len >= ptc_nt_pos:
            offset = ptc_nt_pos - cum - 1
            if strand == '+':
                return cs + offset
            else:
                return ce - offset - 1
        cum += seg_len

    return None


def _tx_order_last_exon_indices(strand, n_exons):
    """Indices (into genomically ascending exon list) of the transcript-order
    last and penultimate exons. For '+' strand the last exon is the highest
    coordinate; for '-' strand it is the lowest."""
    if strand == '+':
        return n_exons - 1, n_exons - 2
    return 0, 1


def _tx_offset(genomic_pos, sorted_exons, strand):
    """Transcript-space offset (nt from the 5' end) of a genomic position.
    Returns None if the position is not inside any exon."""
    offset = 0
    if strand == '+':
        for es, ee in sorted_exons:
            if es <= genomic_pos < ee:
                return offset + (genomic_pos - es)
            offset += ee - es
    else:
        for es, ee in reversed(sorted_exons):
            if es <= genomic_pos < ee:
                return offset + (ee - 1 - genomic_pos)
            offset += ee - es
    return None


def _classify_nmd_by_ptc_position(ptc_genomic_pos, full_tx_exons, strand, start_codon_pos=None):
    if not full_tx_exons:
        return "NMD_Target"

    sorted_exons = sorted(full_tx_exons, key=lambda x: x[0])
    n_exons = len(sorted_exons)

    if n_exons == 1:
        return "NMDER3"

    ptc_exon_idx = -1
    for i, (es, ee) in enumerate(sorted_exons):
        if es <= ptc_genomic_pos < ee:
            ptc_exon_idx = i
            break

    if ptc_exon_idx == -1:
        return "NMD_Target"

    last_idx, penult_idx = _tx_order_last_exon_indices(strand, n_exons)

    # NMDER1: PTC < 150 nt downstream of the start codon (transcript distance)
    if start_codon_pos is not None:
        ptc_off = _tx_offset(ptc_genomic_pos, sorted_exons, strand)
        start_off = _tx_offset(start_codon_pos, sorted_exons, strand)
        if ptc_off is not None and start_off is not None:
            if ptc_off - start_off <= 150:
                return "NMDER1"

    # NMDER3: PTC within the transcript-order last exon
    if ptc_exon_idx == last_idx:
        return "NMDER3"

    # NMDER3: PTC <= 50 nt upstream of the terminal exon-exon junction
    if ptc_exon_idx == penult_idx:
        if strand == '+':
            junction_pos = sorted_exons[penult_idx][1]
            dist_to_terminal_ej = junction_pos - ptc_genomic_pos
        else:
            junction_pos = sorted_exons[penult_idx][0]
            dist_to_terminal_ej = ptc_genomic_pos - junction_pos
        if 0 <= dist_to_terminal_ej <= 50:
            return "NMDER3"

    # NMDER2: PTC > 400 nt from the next downstream exon-exon junction
    if strand == '+':
        downstream_junction_pos = sorted_exons[ptc_exon_idx][1]
        dist_to_downstream_ej = downstream_junction_pos - ptc_genomic_pos
    else:
        downstream_junction_pos = sorted_exons[ptc_exon_idx][0]
        dist_to_downstream_ej = ptc_genomic_pos - downstream_junction_pos

    if dist_to_downstream_ej > 400:
        return "NMDER2"

    return "NMD_Target"


def absolute_nmd_status(prot, cds_ph_truncated, full_tx_exons, strand,
                        start_codon_pos=None, premature_ptc=False):
    """Absolute NMD assessment of an isoform's OWN stop codon, independent of
    whether a splicing event introduced it. Returns:
      Normal        stop codon not predicted to trigger NMD
      NMD_Target    stop codon predicted NMD-sensitive (canonical EJC rule)
      NMDER1/2/3    PREMATURE stop codon in an NMD-escape zone (only reported
                    when premature_ptc=True, i.e. the stop is event-induced;
                    a full-length stop codon sitting in an escape zone — the
                    normal case, e.g. stop in the last exon — is "Normal")
      NO_PROTEIN / NO_STOP_CODON / UNMAPPED when the call cannot be made.
    """
    if not prot or not cds_ph_truncated:
        return "NO_PROTEIN"

    star_pos = prot.find('*')
    if star_pos == -1:
        return "NO_STOP_CODON"

    stop_nt_pos = (star_pos + 1) * 3
    stop_genomic_pos = _map_ptc_to_genomic(cds_ph_truncated, strand, stop_nt_pos)
    if stop_genomic_pos is None:
        return "UNMAPPED"

    cls = _classify_nmd_by_ptc_position(stop_genomic_pos, full_tx_exons, strand,
                                        start_codon_pos)
    if cls in ("NMDER1", "NMDER2", "NMDER3"):
        return cls if premature_ptc else "Normal"
    if cls == "NMD_Target":
        return "NMD_Target"
    return "Normal"


def judge_nmd_status(wt_prot, alt_prot, is_frameshift,
                     alt_cds_ph_truncated, full_tx_exons, strand, start_codon_pos=None):
    """NMD consequence of the EVENT: whether the alternative isoform's ORF
    acquires a premature termination codon relative to the reference.
    Frame-preserving events -> Express; otherwise the alt isoform's premature
    stop codon is classified (NMD_Target / NMDER1-3)."""
    if not wt_prot or not alt_prot:
        return "NO_PROTEIN"

    if not is_frameshift:
        return "Express"

    cls = absolute_nmd_status(alt_prot, alt_cds_ph_truncated, full_tx_exons,
                              strand, start_codon_pos, premature_ptc=True)
    if cls == "UNMAPPED":
        return "NMD_Target"
    if cls == "Normal":
        # a frameshift that ends on the reference stop codon cannot happen;
        # treat defensively as Express
        return "Express"
    return cls


def _stop_key(cds_ph_trunc, strand):
    """Genomic span of the termination codon (the last 3 nt of the CDS in
    transcript order). Two isoforms whose stop codons sit at the same genomic
    position terminate identically regardless of read-through length."""
    if not cds_ph_trunc:
        return None
    if strand == '+':
        s, e, _ = sorted(cds_ph_trunc, key=lambda x: x[0])[-1]
        return (e - 3, e)
    s, e, _ = sorted(cds_ph_trunc, key=lambda x: x[0])[0]
    return (s, s + 3)


def _get_start_codon_position(ref_cds_ph_trunc, strand):
    if not ref_cds_ph_trunc:
        return None
    if strand == '+':
        return ref_cds_ph_trunc[0][0]
    else:
        return ref_cds_ph_trunc[-1][1]


def transcript_self_nmd(cds_ph_list, tx_exons, strand):
    """Whether a CDS-annotated transcript is ITSELF NMD-sensitive according to
    the same EJC rules (its annotated stop codon acting as the PTC). Used to
    avoid anchoring events on transcripts that are already NMD substrates —
    BioSpliceView assumes the reference transcript is normally expressed.

    Returns True if the annotated stop is judged NMD-sensitive; False otherwise
    (including when a judgement cannot be made).
    """
    if not cds_ph_list or not tx_exons:
        return False

    # accept both (start, end) and (start, end, frame) segments
    cds_ph_list = [(seg[0], seg[1], seg[2] if len(seg) > 2 else 0) for seg in cds_ph_list]

    sorted_exons = sorted(tx_exons, key=lambda x: x[0])
    if len(sorted_exons) <= 1:
        return False

    total_cds = sum(ce - cs for cs, ce, _ in cds_ph_list)
    if total_cds < 3:
        return False

    stop_nt_pos = total_cds  # last codon of the annotated CDS (inclusive of stop)
    stop_genomic_pos = _map_ptc_to_genomic(cds_ph_list, strand, stop_nt_pos)
    if stop_genomic_pos is None:
        return False

    return _classify_nmd_by_ptc_position(stop_genomic_pos, sorted_exons, strand) == "NMD_Target"
