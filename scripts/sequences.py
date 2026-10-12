#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genome sequence retrieval and CDS/mRNA sequence operations."""

import logging

from pyfaidx import Fasta
from Bio.Seq import Seq

from .settings import ALT_DOWNSTREAM_EXTENSION_BP

logger = logging.getLogger(__name__)


def normalize_chr(chr_str: str) -> str:
    s = str(chr_str).strip()
    return s[3:] if s.lower().startswith('chr') else s


def safe_get_seq(genome: Fasta, chrom: str, start_0: int, end_0: int) -> str:
    key = normalize_chr(chrom)
    actual_key = None
    for candidate in [key, f'chr{key}', f'Chr{key}']:
        if candidate in genome:
            actual_key = candidate
            break
    if actual_key is None:
        for k in genome.keys():
            if normalize_chr(k) == key:
                actual_key = k
                break
    if actual_key is None:
        return ""
    cl = len(genome[actual_key])
    s, e = max(0, int(start_0)), min(int(end_0), cl)
    return str(genome[actual_key][s:e]) if s < e else ""


def _starts_with_atg(cds_seq: str) -> bool:
    if not cds_seq or len(cds_seq) < 3:
        return False
    return cds_seq[:3].upper() == "ATG"


def get_coding_sequence(cds_segments, genome, chrom, strand):
    if not cds_segments:
        return ""
    if strand == '+':
        parts = [safe_get_seq(genome, chrom, s, e) for s, e, _ in cds_segments]
        return "".join(parts)
    else:
        parts = []
        for s, e, _ in reversed(cds_segments):
            seq = safe_get_seq(genome, chrom, s, e)
            parts.append(str(Seq(seq).reverse_complement()))
        return "".join(parts)


def get_exon_sequence(exon_segments, genome, chrom, strand):
    if not exon_segments:
        return ""
    if strand == '+':
        parts = [safe_get_seq(genome, chrom, s, e) for s, e in exon_segments]
        return "".join(parts)
    else:
        parts = []
        for s, e in reversed(exon_segments):
            seq = safe_get_seq(genome, chrom, s, e)
            parts.append(str(Seq(seq).reverse_complement()))
        return "".join(parts)


def extend_cds_with_3utr(cds_ph_list, genome, chrom, strand, extension_bp=None, tx_exons=None):
    """Append 3'UTR sequence to the CDS chain so the ORF can run through to a
    stop codon. When tx_exons (the isoform's exon chain, 0-based half-open) is
    provided the extension follows the SPLICED chain — jumping introns — which
    is required for alternative isoforms; otherwise it falls back to
    genomically contiguous extension (safe for reference CDS whose annotated
    stop codon immediately follows the last CDS block)."""
    if extension_bp is None:
        extension_bp = ALT_DOWNSTREAM_EXTENSION_BP
    if not cds_ph_list or extension_bp <= 0:
        return list(cds_ph_list)
    extended = list(cds_ph_list)
    norm_key = normalize_chr(chrom)
    cl = 0
    for candidate in [norm_key, f'chr{norm_key}', f'Chr{norm_key}']:
        if candidate in genome:
            cl = len(genome[candidate])
            break
    if cl == 0:
        for k in genome.keys():
            if normalize_chr(k) == norm_key:
                cl = len(genome[k])
                break
    cum_len = sum(e-s for s, e, _ in cds_ph_list)

    if tx_exons:
        exons = sorted(tx_exons)
        ext_segs = []
        remaining = extension_bp
        if strand == '+':
            last_e = cds_ph_list[-1][1]
            # 当前 CDS 末端所在外显子及其后继外显子（转录本序 = 升序）
            for i, (es, ee) in enumerate(exons):
                if es <= last_e <= ee:
                    if last_e < ee and remaining > 0:
                        take = min(ee - last_e, remaining)
                        ext_segs.append((last_e, last_e + take, cum_len % 3))
                        remaining -= take
                    for nes, nee in exons[i+1:]:
                        if remaining <= 0:
                            break
                        take = min(nee - nes, remaining)
                        ext_segs.append((nes, nes + take, cum_len % 3))
                        remaining -= take
                    break
        else:
            first_s = cds_ph_list[0][0]
            for i, (es, ee) in enumerate(exons):
                if es <= first_s <= ee:
                    if es < first_s and remaining > 0:
                        take = min(first_s - es, remaining)
                        ext_segs.append((first_s - take, first_s, cum_len % 3))
                        remaining -= take
                    for nes, nee in reversed(exons[:i]):
                        if remaining <= 0:
                            break
                        take = min(nee - nes, remaining)
                        ext_segs.append((nee - take, nee, cum_len % 3))
                        remaining -= take
                    break
        if ext_segs:
            extended.extend(ext_segs)
        return sorted(extended, key=lambda x: x[0])

    # 无外显子链时的回退：基因组连续延伸（参考转录本注释 stop 紧邻 CDS 末端时安全）
    if strand == '+':
        last_e = cds_ph_list[-1][1]
        ext_end = min(last_e + extension_bp, cl) if cl > 0 else last_e + extension_bp
        if last_e < ext_end:
            extended.append((last_e, ext_end, cum_len % 3))
    else:
        first_s = cds_ph_list[0][0]
        ext_start = max(0, first_s - extension_bp)
        if ext_start < first_s:
            extended.insert(0, (ext_start, first_s, cum_len % 3))
    return sorted(extended, key=lambda x: x[0])


def truncate_cds_at_stop(cds_ph_list, protein, strand):
    if not cds_ph_list or not protein:
        return cds_ph_list, protein
    star_pos = protein.find('*')
    if star_pos < 0:
        return cds_ph_list, protein
    stop_nt_end = (star_pos + 1) * 3
    if strand == '+':
        ordered = sorted(cds_ph_list, key=lambda x: x[0])
    else:
        ordered = sorted(cds_ph_list, key=lambda x: x[0], reverse=True)
    truncated_ordered = []
    cum = 0
    for cs, ce, ph in ordered:
        seg_len = ce - cs
        if cum >= stop_nt_end:
            break
        if cum + seg_len <= stop_nt_end:
            truncated_ordered.append((cs, ce, ph))
        else:
            keep_nt = stop_nt_end - cum
            if keep_nt > 0:
                if strand == '+':
                    truncated_ordered.append((cs, cs + keep_nt, ph))
                else:
                    truncated_ordered.append((ce - keep_nt, ce, ph))
            break
        cum += seg_len
    return sorted(truncated_ordered, key=lambda x: x[0]), protein[:star_pos + 1]


def _compute_trimmed_diff(ref_seq, alt_seq):
    if not ref_seq or not alt_seq or ref_seq == alt_seq:
        return ""
    left = 0
    ml = min(len(ref_seq), len(alt_seq))
    while left < ml and ref_seq[left] == alt_seq[left]:
        left += 1
    right = 0
    while right < ml - left and ref_seq[-(right+1)] == alt_seq[-(right+1)]:
        right += 1
    rd = ref_seq[left:-right] if right else ref_seq[left:]
    ad = alt_seq[left:-right] if right else alt_seq[left:]
    if not rd and not ad:
        return ""
    return f"{rd}|{ad}"
