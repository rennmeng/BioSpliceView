#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rMATS event coordinate parsing, frameshift detection and table schema constants."""

import logging
from typing import List, Tuple

logger = logging.getLogger(__name__)

# Annotation columns appended after the raw rMATS columns in the output Excel
ANNOTATION_COLUMNS = [
    'Ref_Transcript_ID', 'Anchor_Biotype', 'Ref_Full_CDS', 'Ref_Full_Protein', 'Ref_CDS_Length',
    'Alt_Full_CDS', 'Alt_Full_Protein', 'Alt_CDS_Length',
    'Diff_CDS', 'Diff_Protein', 'Diff_Length',
    'NMD_Status', 'Ref_Absolute_NMD', 'Alt_Absolute_NMD', 'Match_Status'
]

# Coordinate columns of the rMATS output (selected per event type)
RMATS_COORD_COLS = [
    'exonStart_0base', 'exonEnd', 'riExonStart_0base', 'riExonEnd',
    '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd',
    'shortES', 'shortEE', 'longExonStart_0base', 'longExonEnd',
    'flankingES', 'flankingEE', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
]

# Metadata columns of the rMATS output
RMATS_META_COLS = [
    'chr', 'strand', 'GeneID', 'geneSymbol',
    'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2',
    'IncFormLen', 'SkipFormLen', 'PValue', 'FDR',
    'IncLevel1', 'IncLevel2', 'IncLevelDifference'
]


def _extract_event_regions(row, etype) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """
    Return (anchors, diff_regions) as genomic 0-based half-open intervals.
    anchors: the "invariant" regions used to match the reference transcript
    diff_regions: the event difference regions
    """
    anchors, diff_regions = [], []
    if etype == "SE":
        anchors.append((int(float(row['upstreamES'])), int(float(row['upstreamEE']))))
        anchors.append((int(float(row['downstreamES'])), int(float(row['downstreamEE']))))
        diff_regions.append((int(float(row['exonStart_0base'])), int(float(row['exonEnd']))))
    elif etype == "RI":
        up_es, up_ee = int(float(row['upstreamES'])), int(float(row['upstreamEE']))
        dn_es, dn_ee = int(float(row['downstreamES'])), int(float(row['downstreamEE']))
        anchors.extend([(up_es, up_ee), (dn_es, dn_ee)])
        if up_ee < dn_es:
            diff_regions.append((up_ee, dn_es))
    elif etype == "MXE":
        diff_regions.append((int(float(row['1stExonStart_0base'])), int(float(row['1stExonEnd']))))
        diff_regions.append((int(float(row['2ndExonStart_0base'])), int(float(row['2ndExonEnd']))))
        if 'upstreamES' in row and 'upstreamEE' in row:
            anchors.append((int(float(row['upstreamES'])), int(float(row['upstreamEE']))))
        if 'downstreamES' in row and 'downstreamEE' in row:
            anchors.append((int(float(row['downstreamES'])), int(float(row['downstreamEE']))))
    elif etype in ("A5SS", "A3SS"):
        short_es = int(float(row['shortES']))
        short_ee = int(float(row['shortEE']))
        long_es = int(float(row['longExonStart_0base']))
        long_ee = int(float(row['longExonEnd']))
        flank_es = int(float(row['flankingES']))
        flank_ee = int(float(row['flankingEE']))

        # 可变边界随链变化（A3SS 在 + 链共享终点、在 - 链共享起点，A5SS 相反），
        # 按几何关系判定，而不是按事件类型硬编码
        if long_es == short_es:          # 起点共享 -> 终点可变
            diff_regions.append((min(short_ee, long_ee), max(short_ee, long_ee)))
        elif long_ee == short_ee:        # 终点共享 -> 起点可变
            diff_regions.append((min(short_es, long_es), max(short_es, long_es)))
        elif etype == "A5SS":
            diff_regions.append((min(short_es, long_es), max(short_es, long_es)))
        else:
            diff_regions.append((min(short_ee, long_ee), max(short_ee, long_ee)))

        anchors.append((short_es, short_ee))
        anchors.append((flank_es, flank_ee))
    return anchors, diff_regions


def get_event_genomic_length(row, etype):
    """Genomic length of the event difference region (MXE: difference between the two exclusive exons)."""
    if etype == "SE":
        return abs(int(float(row['exonEnd'])) - int(float(row['exonStart_0base'])))
    elif etype == "RI":
        return abs(int(float(row['downstreamES'])) - int(float(row['upstreamEE'])))
    elif etype == "MXE":
        e1 = int(float(row['1stExonEnd'])) - int(float(row['1stExonStart_0base']))
        e2 = int(float(row['2ndExonEnd'])) - int(float(row['2ndExonStart_0base']))
        return abs(e2 - e1)
    elif etype in ("A5SS", "A3SS"):
        long_len = int(float(row['longExonEnd'])) - int(float(row['longExonStart_0base']))
        short_len = int(float(row['shortEE'])) - int(float(row['shortES']))
        return abs(long_len - short_len)
    return 0


def _region_overlaps_cds(region, cds_sorted):
    rs, re_ = region
    if rs >= re_:
        return False
    for cs, ce, _ in cds_sorted:
        if max(cs, rs) < min(ce, re_):
            return True
    return False


def _is_frameshift_by_event(etype, row, ref_cds_ph=None, strand='+'):
    try:
        if not ref_cds_ph:
            return False

        cds_sorted = sorted(ref_cds_ph, key=lambda x: x[0])

        change_region = None
        diff_len = 0

        if etype == "SE":
            s = int(float(row['exonStart_0base']))
            e = int(float(row['exonEnd']))
            if e <= s:
                return False
            change_region = (s, e)
            diff_len = e - s

        elif etype == "RI":
            up_ee = int(float(row['upstreamEE']))
            dn_es = int(float(row['downstreamES']))
            lo, hi = min(up_ee, dn_es), max(up_ee, dn_es)
            if hi <= lo:
                return False
            change_region = (lo, hi)
            diff_len = hi - lo

        elif etype == "MXE":
            e1_s = int(float(row['1stExonStart_0base']))
            e1_e = int(float(row['1stExonEnd']))
            e2_s = int(float(row['2ndExonStart_0base']))
            e2_e = int(float(row['2ndExonEnd']))
            e1_len = e1_e - e1_s
            e2_len = e2_e - e2_s
            if e1_len <= 0 or e2_len <= 0:
                return False
            change_region = (min(e1_s, e2_s), max(e1_e, e2_e))
            diff_len = abs(e1_len - e2_len)

        elif etype == "A5SS":
            short_s = int(float(row['shortES']))
            short_e = int(float(row['shortEE']))
            long_s = int(float(row['longExonStart_0base']))
            long_e = int(float(row['longExonEnd']))
            short_len = short_e - short_s
            long_len = long_e - long_s
            if short_len <= 0 or long_len <= 0:
                return False
            change_region = (min(short_s, long_s), max(short_e, long_e))
            diff_len = abs(short_len - long_len)

        elif etype == "A3SS":
            short_s = int(float(row['shortES']))
            short_e = int(float(row['shortEE']))
            long_s = int(float(row['longExonStart_0base']))
            long_e = int(float(row['longExonEnd']))
            short_len = short_e - short_s
            long_len = long_e - long_s
            if short_len <= 0 or long_len <= 0:
                return False
            change_region = (min(short_s, long_s), max(short_e, long_e))
            diff_len = abs(short_len - long_len)

        else:
            return True

        if change_region is None:
            return False

        if not _region_overlaps_cds(change_region, cds_sorted):
            return False

        if diff_len <= 0:
            return False
        return diff_len % 3 != 0

    except (KeyError, ValueError, TypeError) as e:
        logger.warning(f"Frameshift check failed for {etype}: {e}")
        return True
