#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GenBank record generation.

Two GenBank records are written per event:
  - EXON-*.gb : mRNA sequence spliced from the reference transcript exons
                (exon/CDS features plus event-landmark misc_features)
  - Tran-*.gb : genomic sequence spanned by the reference transcript
                (exon/CDS features plus event-landmark misc_features)

Landmark labels follow the manuscript vocabulary (anchor exons, skipped
exon, retained intron, mutually exclusive exons, alternative splice site
region) and every record carries a structured comment reporting the NMD
consequence of the event and the absolute NMD status of each isoform.
"""

import logging
import os
import re
from datetime import datetime

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO

from .events import _extract_event_regions
from .matching import _exon_in_ref
from .sequences import get_exon_sequence, safe_get_seq

logger = logging.getLogger(__name__)


def sanitize_filename(s: str) -> str:
    return re.sub(r'[^\w\-\.]', '_', str(s).strip())


def _gb1(s, e):
    """Convert a 0-based half-open interval to a 1-based inclusive string."""
    return f"{s + 1}-{e}"


def _nmd_comment(context, gene_symbol, etype, event_id, ref_tid):
    """Structured record comment carrying the dual NMD report."""
    if not context:
        return (f"BioSpliceView {etype} event {event_id} of {gene_symbol}; "
                f"reference transcript {ref_tid}")
    lines = [
        f"BioSpliceView {etype} event {event_id} of {gene_symbol}",
        f"Reference transcript: {ref_tid}"
        + (f" [{context.get('anchor_biotype')}]" if context.get('anchor_biotype') else ""),
        f"NMD consequence of event: {context.get('nmd_status', 'NA')} "
        "(Express = no PTC introduced; NMD_Target = predicted NMD-sensitive; "
        "NMDER1/2/3 = PTC with predicted NMD escape)",
        f"Reference isoform NMD status (absolute): {context.get('ref_abs_nmd', 'NA')}",
        f"Alternative isoform NMD status (absolute): {context.get('alt_abs_nmd', 'NA')}",
    ]
    return " | ".join(lines)


def _event_landmarks(etype, row, anchors, diff_regions, ref_exons):
    """Build (label, note, interval) landmark triples using the manuscript
    vocabulary. `ref_exons` decides which event form the reference realises."""
    labels = []

    def form_note(in_ref):
        return ("present in reference isoform" if in_ref
                else "specific to alternative isoform")

    if etype == "SE":
        if len(anchors) >= 2:
            labels.append(("anchor exon (upstream)", "constitutive exon upstream of the event",
                           anchors[0]))
            labels.append(("anchor exon (downstream)", "constitutive exon downstream of the event",
                           anchors[1]))
        for d in diff_regions:
            labels.append(("skipped exon",
                           f"cassette exon; {form_note(_exon_in_ref(d, ref_exons))}", d))

    elif etype == "RI":
        if len(anchors) >= 2:
            labels.append(("anchor exon (upstream)", "exon upstream of the retained intron",
                           anchors[0]))
            labels.append(("anchor exon (downstream)", "exon downstream of the retained intron",
                           anchors[1]))
        for d in diff_regions:
            labels.append(("retained intron",
                           "intron retained in the alternative isoform", d))

    elif etype == "MXE":
        if len(anchors) >= 2:
            labels.append(("anchor exon (upstream)", "constitutive exon upstream of the event",
                           anchors[0]))
            labels.append(("anchor exon (downstream)", "constitutive exon downstream of the event",
                           anchors[1]))
        for i, d in enumerate(diff_regions):
            labels.append((f"mutually exclusive exon {i + 1}",
                           f"{form_note(_exon_in_ref(d, ref_exons))}", d))

    elif etype in ("A5SS", "A3SS"):
        short = (int(float(row['shortES'])), int(float(row['shortEE'])))
        long_ = (int(float(row['longExonStart_0base'])), int(float(row['longExonEnd'])))
        flank = (int(float(row['flankingES'])), int(float(row['flankingEE'])))
        labels.append(("short isoform exon",
                       f"exon using the short-form splice site; "
                       f"{form_note(_exon_in_ref(short, ref_exons))}", short))
        labels.append(("long isoform exon",
                       f"exon using the long-form splice site; "
                       f"{form_note(_exon_in_ref(long_, ref_exons))}", long_))
        labels.append(("flanking exon", "constitutive exon flanking the event", flank))
        for d in diff_regions:
            labels.append(("alternative splice site region",
                           f"sequence differing between the short and long forms of this "
                           f"{etype} event", d))

    return labels


def generate_exon_gb(ref_tid, exon_bounds, cds_ph_list, genome, chrom, strand,
                     gene_symbol, etype, row, output_dir, job_id, event_id="", context=None):
    if not exon_bounds or not ref_tid:
        return None
    try:
        is_minus = (strand == '-')
        n_exons = len(exon_bounds)

        mrna_seq = get_exon_sequence(exon_bounds, genome, chrom, strand)
        if not mrna_seq:
            return None
        seq_len = len(mrna_seq)
        feat_strand = 1

        def _mirror(pos):
            return seq_len - pos if is_minus else pos

        safe_gene = sanitize_filename(gene_symbol)
        id_prefix = f"{job_id}_" if job_id else ""
        locus_name = f"{id_prefix}{event_id}_EXON-{etype}-{safe_gene}"

        record = SeqRecord(Seq(mrna_seq), id=ref_tid, name=locus_name[:16],
                           description=f"{gene_symbol} {etype} event {event_id} spliced mRNA "
                                       f"(reference isoform)")
        record.annotations["molecule_type"] = "mRNA"
        record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
        record.annotations["keywords"] = ["alternative splicing", etype, "BioSpliceView"]
        record.annotations["comment"] = _nmd_comment(context, gene_symbol, etype, event_id, ref_tid)

        cum = 0
        for idx, (s, e) in enumerate(exon_bounds):
            elen = e - s
            rel_s, rel_e = cum, cum + elen
            if is_minus:
                rel_s, rel_e = _mirror(cum + elen), _mirror(cum)
            tx_num = (n_exons - idx) if is_minus else (idx + 1)
            record.features.append(SeqFeature(
                FeatureLocation(rel_s, rel_e, strand=feat_strand),
                type="exon", qualifiers={"number": str(tx_num)}))
            cum += elen

        if cds_ph_list:
            cds_sorted = sorted(cds_ph_list, key=lambda x: x[0])
            cds_rel_parts = []
            for cs, ce, _ in cds_sorted:
                c = 0
                for es, ee in exon_bounds:
                    elen = ee - es
                    if cs >= es and ce <= ee:
                        off = cs - es
                        rs, re_ = c + off, c + off + (ce - cs)
                        if is_minus:
                            rs, re_ = _mirror(c + off + (ce - cs)), _mirror(c + off)
                        cds_rel_parts.append((rs, re_))
                        break
                    elif cs < ee and ce > es:
                        os_ = max(cs, es)
                        oe = min(ce, ee)
                        off = os_ - es
                        rs, re_ = c + off, c + off + (oe - os_)
                        if is_minus:
                            rs, re_ = _mirror(c + off + (oe - os_)), _mirror(c + off)
                        cds_rel_parts.append((rs, re_))
                    c += elen

            if cds_rel_parts:
                cds_rel_parts.sort(key=lambda x: x[0])
                locs = [FeatureLocation(s, e, strand=feat_strand) for s, e in cds_rel_parts]
                cds_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)

                cds_seq = "".join(mrna_seq[s:e] for s, e in cds_rel_parts)
                translated = str(Seq(cds_seq).translate(to_stop=False)) if cds_seq else ""
                sp = translated.find('*')
                if sp >= 0:
                    translated = translated[:sp + 1]

                quals = {
                    "gene": gene_symbol,
                    "transcript_id": ref_tid,
                    "codon_start": "1",
                    "transl_table": "1",
                    "product": f"{gene_symbol} protein (reference isoform)"
                }
                if translated:
                    quals["translation"] = translated

                record.features.append(SeqFeature(cds_loc, type="CDS", qualifiers=quals))
        else:
            base_comment = record.annotations.get("comment", "")
            record.annotations["comment"] = (
                f"{base_comment} | Non-coding RNA transcript for {gene_symbol}")

        try:
            anchors, diff_regions = _extract_event_regions(row, etype)
            event_labels = _event_landmarks(etype, row, anchors, diff_regions, exon_bounds)

            for label, note_text, (gs, ge) in event_labels:
                if gs >= ge:  # zero-length region (e.g. A3SS/A5SS with identical sites)
                    continue
                rel_parts = []
                c = 0
                found = False
                for es, ee in exon_bounds:
                    elen = ee - es
                    os_ = max(gs, es)
                    oe = min(ge, ee)
                    if os_ < oe:
                        off = os_ - es
                        rs, re_ = c + off, c + off + (oe - os_)
                        if is_minus:
                            rs, re_ = _mirror(re_), _mirror(rs)
                        rel_parts.append((rs, re_))
                        found = True
                    c += elen
                if found and rel_parts:
                    rel_parts.sort(key=lambda x: x[0])
                    merged = [rel_parts[0]]
                    for s, e in rel_parts[1:]:
                        if s <= merged[-1][1]:
                            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                        else:
                            merged.append((s, e))
                    locs = [FeatureLocation(s, e, strand=feat_strand) for s, e in merged]
                    feat_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
                    record.features.append(SeqFeature(
                        feat_loc, type="misc_feature",
                        qualifiers={
                            "label": label,
                            "note": f"{etype} landmark; {note_text}; "
                                    f"genomic {chrom}:{_gb1(gs, ge)} (1-based)"}))
                else:
                    existing_note = record.annotations.get("comment", "")
                    add_note = f"{label}: {chrom}:{_gb1(gs, ge)} (not present in this reference transcript)"
                    record.annotations["comment"] = (f"{existing_note}; {add_note}"
                                                     if existing_note else add_note)
        except Exception as e:
            logger.debug(f"EXON-GB misc_feature annotation failed: {e}")

        gb_path = os.path.join(output_dir, f"{locus_name}.gb")
        SeqIO.write(record, gb_path, "genbank")
        return gb_path

    except Exception as e:
        logger.warning(f"EXON-GB failed {ref_tid}: {e}")
        return None


def generate_tran_gb(ref_tid, full_tx_exons, cds_ph_list, genome, chrom, strand,
                     gene_symbol, etype, row, output_dir, job_id, event_id="", context=None):
    if not ref_tid:
        return None
    try:
        all_coords = []
        if full_tx_exons:
            for s, e in full_tx_exons:
                all_coords.extend([s, e])
        if cds_ph_list:
            for cs, ce, _ in cds_ph_list:
                all_coords.extend([cs, ce])
        if not all_coords:
            return None

        gene_start = min(all_coords)
        gene_end = max(all_coords)
        genomic_seq = safe_get_seq(genome, chrom, gene_start, gene_end)
        if not genomic_seq:
            return None

        is_minus = (strand == '-')
        display_seq = str(Seq(genomic_seq).reverse_complement()) if is_minus else genomic_seq
        seq_len = len(display_seq)
        feat_strand = 1

        def _mirror(pos):
            return seq_len - pos

        safe_gene = sanitize_filename(gene_symbol)
        id_prefix = f"{job_id}_" if job_id else ""
        locus_name = f"{id_prefix}{event_id}_Tran-{etype}-{safe_gene}"

        record = SeqRecord(Seq(display_seq), id=ref_tid, name=locus_name[:16],
                           description=f"{gene_symbol} {etype} event {event_id} genomic context "
                                       f"(reference transcript)")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
        record.annotations["keywords"] = ["alternative splicing", etype, "BioSpliceView"]
        record.annotations["comment"] = _nmd_comment(context, gene_symbol, etype, event_id, ref_tid)

        if full_tx_exons:
            n_exons = len(full_tx_exons)
            for idx, (s, e) in enumerate(full_tx_exons):
                cs = max(gene_start, min(gene_end, s))
                ce = max(gene_start, min(gene_end, e))
                if cs >= ce:
                    continue
                rel_s, rel_e = cs - gene_start, ce - gene_start
                if is_minus:
                    rel_s, rel_e = _mirror(rel_e), _mirror(rel_s)
                tx_num = (n_exons - idx) if is_minus else (idx + 1)
                record.features.append(SeqFeature(
                    FeatureLocation(rel_s, rel_e, strand=feat_strand),
                    type="exon", qualifiers={"number": str(tx_num)}))

        if cds_ph_list:
            raw_parts = []
            for gs, ge, _ in sorted(cds_ph_list, key=lambda x: x[0]):
                cs = max(gene_start, min(gene_end, gs))
                ce = max(gene_start, min(gene_end, ge))
                if cs >= ce:
                    continue
                rel_s, rel_e = cs - gene_start, ce - gene_start
                if is_minus:
                    rel_s, rel_e = _mirror(rel_e), _mirror(rel_s)
                raw_parts.append((rel_s, rel_e))
            raw_parts.sort(key=lambda x: x[0])
            merged = []
            for s, e in raw_parts:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            if merged:
                locs = [FeatureLocation(s, e, strand=feat_strand) for s, e in merged]
                cds_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
                record.features.append(SeqFeature(
                    cds_loc, type="CDS",
                    qualifiers={"gene": gene_symbol, "transcript_id": ref_tid,
                                "codon_start": "1", "transl_table": "1",
                                "product": f"{gene_symbol} protein (reference isoform)"}))

        try:
            anchors, diff_regions = _extract_event_regions(row, etype)
            event_labels = _event_landmarks(etype, row, anchors, diff_regions, full_tx_exons)

            for label, note_text, (gs, ge) in event_labels:
                if gs >= ge:
                    continue
                cs = max(gene_start, min(gene_end, gs))
                ce = max(gene_start, min(gene_end, ge))
                if cs < ce:
                    rel_s, rel_e = cs - gene_start, ce - gene_start
                    if is_minus:
                        rel_s, rel_e = _mirror(rel_e), _mirror(rel_s)
                    record.features.append(SeqFeature(
                        FeatureLocation(rel_s, rel_e, strand=feat_strand),
                        type="misc_feature",
                        qualifiers={
                            "label": label,
                            "note": f"{etype} landmark; {note_text}; "
                                    f"genomic {chrom}:{_gb1(gs, ge)} (1-based)"}))
        except Exception as e:
            logger.debug(f"Tran-GB misc_feature annotation failed: {e}")

        gb_path = os.path.join(output_dir, f"{locus_name}.gb")
        SeqIO.write(record, gb_path, "genbank")
        return gb_path

    except Exception as e:
        logger.warning(f"Tran-GB failed {ref_tid}: {e}")
        return None
