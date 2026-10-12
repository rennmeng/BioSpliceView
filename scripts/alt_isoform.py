#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Derive the alternative-isoform CDS segments from the reference CDS.

Both input and output are lists of half-open segments
[(cds_start_0, cds_end_0, phase), ...].

Coordinate convention (both sides are 0-based half-open and compared with
EXACT equality, no tolerance):
  - GTF index (gtf_index.build_gtf_index): (gtf_start - 1, gtf_end)
  - rMATS event columns: an rMATS interval (columnStart_0base, columnEnd) maps
    directly to 0-based half-open, because rMATS starts are 0-based and rMATS
    ends are 1-based inclusive (= 0-based half-open end). No ±1 shift is
    applied anywhere in this module; a 1-based GTF would silently fail the
    exact matches and fall back to the documented heuristics.
"""


def build_alt_exon_chain(anchor_exons, row, etype):
    """Apply the event to the anchor transcript's exon chain and return the
    alternative isoform's exons (0-based half-open, genomically sorted).
    Exon identity requires EXACT coordinate equality with the rMATS interval.
    Returns None when the transformation is ambiguous (caller should then
    keep the anchor chain as a fallback)."""
    if not anchor_exons:
        return None
    exons = sorted(anchor_exons)

    if etype == "SE":
        cas = (int(float(row['exonStart_0base'])), int(float(row['exonEnd'])))
        if cas in exons:
            i = exons.index(cas)
            return exons[:i] + exons[i+1:]                     # 锚定含盒 -> Alt 去盒
        return sorted(exons + [cas])                           # 锚定为跳过形 -> Alt 加盒

    if etype == "RI":
        iS, iE = int(float(row['upstreamEE'])), int(float(row['downstreamES']))
        if iS >= iE:
            return None
        up_i = next((k for k, (es, ee) in enumerate(exons) if ee == iS), None)
        dn_i = next((k for k, (es, ee) in enumerate(exons) if es == iE), None)
        if up_i is not None and dn_i is not None and dn_i == up_i + 1:
            merged = (exons[up_i][0], exons[dn_i][1])
            return exons[:up_i] + [merged] + exons[dn_i+1:]    # 锚定为剪接形 -> Alt 融合
        big = next((k for k, (es, ee) in enumerate(exons) if es <= iS and ee >= iE), None)
        if big is not None:
            es, ee = exons[big]
            return (exons[:big] + [(es, iS), (iE, ee)] + exons[big+1:])  # 锚定为滞留形 -> 拆分
        return None

    if etype == "MXE":
        e1 = (int(float(row['1stExonStart_0base'])), int(float(row['1stExonEnd'])))
        e2 = (int(float(row['2ndExonStart_0base'])), int(float(row['2ndExonEnd'])))
        has_e1 = e1 in exons
        has_e2 = e2 in exons
        if has_e1 and not has_e2:
            i1 = exons.index(e1)
            return sorted(exons[:i1] + exons[i1+1:] + [e2])
        if has_e2 and not has_e1:
            i2 = exons.index(e2)
            return sorted(exons[:i2] + exons[i2+1:] + [e1])
        return None                                            # 锚定两者皆含/皆无 -> 歧义

    if etype in ("A5SS", "A3SS"):
        short = (int(float(row['shortES'])), int(float(row['shortEE'])))
        long_ = (int(float(row['longExonStart_0base'])), int(float(row['longExonEnd'])))
        flank_es = int(float(row['flankingES']))
        flank_ee = int(float(row['flankingEE']))
        short_s, short_e = short
        long_s, long_e = long_

        # 可变区 V = 长形相对短形多出的区段（长形可能把上游小外显子+内含子一并融合，
        # 因此换形时必须删除 V 内的全部锚定外显子，而不只是替换等长的那一个）
        if long_s == short_s:                 # 共享起点 → 终点可变
            v_s, v_e = short_e, long_e
        elif long_e == short_e:               # 共享终点 → 起点可变
            v_s, v_e = long_s, short_s
        else:
            return None

        # 形式判定：剪接连接证据优先（锚定是否存在 flanking→长/短 形剪接位点），
        # 其次外显子精确匹配。连接 = 相邻外显子 (prev_ee, next_es)。
        J = {(exons[i][1], exons[i+1][0]) for i in range(len(exons) - 1)}
        if flank_ee <= min(long_s, short_s):      # flanking 在低坐标侧 → 连接指向起点
            uses_long = (flank_ee, long_s) in J
            uses_short = (flank_ee, short_s) in J
        elif flank_es >= max(long_e, short_e):    # flanking 在高坐标侧 → 连接来自终点
            uses_long = (long_e, flank_es) in J
            uses_short = (short_e, flank_es) in J
        else:
            uses_long = uses_short = None

        if uses_long is not None and uses_long != uses_short:
            ref_is_long = uses_long
        elif long_ in exons and short not in exons:
            ref_is_long = True
        elif short in exons and long_ not in exons:
            ref_is_long = False
        else:
            return None

        ref_exon = long_ if ref_is_long else short
        alt_exon = short if ref_is_long else long_
        # Alt 链 = 锚定外显子 − ref 形式外显子 − V 区内的外显子(被长形融合/被短形跳过)
        new_chain = [ex for ex in exons
                     if ex != ref_exon and not (min(ex[1], v_e) > max(ex[0], v_s))]
        if alt_exon not in new_chain:
            new_chain.append(alt_exon)
        return sorted(new_chain)

    return None


def compute_alt_from_ref(ref_cds_ph, row, etype, strand, tx_exons=None):
    """tx_exons: anchor transcript exons [(start_0, end_0), ...] — used to decide
    which event form (short/long exon) the reference realises; falls back to
    CDS-overlap heuristics when unavailable."""
    ref_sorted = sorted(ref_cds_ph, key=lambda x: x[0])

    if etype == "RI":
        up_end = int(float(row['upstreamEE']))
        dn_start = int(float(row['downstreamES']))
        if up_end >= dn_start:
            return ref_sorted
        new_cds, inserted = [], False
        for cs, ce, ph in ref_sorted:
            if not inserted and cs >= dn_start:
                new_cds.append((up_end, dn_start, 0))
                inserted = True
            if not inserted and cs < up_end < ce:
                new_cds.append((cs, up_end, ph))
                new_cds.append((up_end, dn_start, 0))
                inserted = True
                new_cds.append((dn_start, ce, 0))
                continue
            new_cds.append((cs, ce, ph))
        if not inserted:
            new_cds.append((up_end, dn_start, 0))
        return sorted(new_cds, key=lambda x: x[0])

    elif etype == "SE":
        skip_s = int(float(row['exonStart_0base']))
        skip_e = int(float(row['exonEnd']))
        new_cds = []
        for cs, ce, ph in ref_sorted:
            os_ = max(cs, skip_s)
            oe = min(ce, skip_e)
            if os_ >= oe:
                new_cds.append((cs, ce, ph))
            else:
                if cs < os_:
                    new_cds.append((cs, os_, ph))
                if ce > oe:
                    new_cds.append((oe, ce, ph))
        return sorted(new_cds, key=lambda x: x[0]) if new_cds else ref_sorted

    elif etype == "MXE":
        # ==================== MXE handling v2 (exclusive-exon fix) ====================
        e1_s = int(float(row['1stExonStart_0base']))
        e1_e = int(float(row['1stExonEnd']))
        e2_s = int(float(row['2ndExonStart_0base']))
        e2_e = int(float(row['2ndExonEnd']))

        def _ov_len(rs, re_):
            tot = 0
            for cs, ce, _ in ref_sorted:
                o_s, o_e = max(cs, rs), min(ce, re_)
                if o_s < o_e:
                    tot += (o_e - o_s)
            return tot

        ov1 = _ov_len(e1_s, e1_e)
        ov2 = _ov_len(e2_s, e2_e)

        if ov1 == 0 and ov2 == 0:
            if strand == '+':
                before_len = sum(ce - cs for cs, ce, _ in ref_sorted if ce <= e2_s)
            else:
                before_len = sum(ce - cs for cs, ce, _ in ref_sorted if cs >= e2_e)
            new_cds = list(ref_sorted)
            new_cds.append((e2_s, e2_e, before_len % 3))
            return sorted(new_cds, key=lambda x: x[0])

        if ov1 >= ov2:
            ref_s, ref_e = e1_s, e1_e
            alt_s, alt_e = e2_s, e2_e
        else:
            ref_s, ref_e = e2_s, e2_e
            alt_s, alt_e = e1_s, e1_e

        kept = []
        for cs, ce, ph in ref_sorted:
            o_s, o_e = max(cs, ref_s), min(ce, ref_e)
            if o_s >= o_e:
                kept.append((cs, ce, ph))
                continue
            if cs < o_s:
                kept.append((cs, o_s, ph))
            if ce > o_e:
                kept.append((o_e, ce, ph))

        has_alt_in_ref = any(max(cs, alt_s) < min(ce, alt_e) for cs, ce, _ in kept)

        if not has_alt_in_ref:
            if strand == '+':
                before_len = sum(ce - cs for cs, ce, _ in kept if ce <= alt_s)
            else:
                before_len = sum(ce - cs for cs, ce, _ in kept if cs >= alt_e)
            kept.append((alt_s, alt_e, before_len % 3))

        kept.sort(key=lambda x: x[0])

        merged = []
        for s, e, ph in kept:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e), merged[-1][2])
            else:
                merged.append((s, e, ph))
        return merged

    elif etype in ("A5SS", "A3SS"):
        short_s = int(float(row['shortES']))
        short_e = int(float(row['shortEE']))
        long_s = int(float(row['longExonStart_0base']))
        long_e = int(float(row['longExonEnd']))
        flank_es = int(float(row['flankingES']))
        flank_ee = int(float(row['flankingEE']))

        # ---- 1) 可变区 V = 长形相对短形多出的区段。长形可能把上游小外显子
        #      与内含子一并融合（含终止密码子），因此换形 = 增/删整个 V 区 ----
        if long_s == short_s:                 # 共享起点 → 终点可变
            v_s, v_e = short_e, long_e
        elif long_e == short_e:               # 共享终点 → 起点可变
            v_s, v_e = long_s, short_s
        else:
            return ref_sorted

        # ---- 2) Ref 形式判定：剪接连接证据优先（锚定是否存在 flanking→
        #      长/短 形剪接位点），其次外显子精确匹配，最后 CDS 重叠(平票→长形)。
        #      不做包含判定（长形必然包含短形区域）----
        junctions = set()
        if tx_exons:
            exons_sorted = sorted(tx_exons)
            junctions = {(exons_sorted[i][1], exons_sorted[i+1][0])
                         for i in range(len(exons_sorted) - 1)}
        if junctions:
            if flank_ee <= min(long_s, short_s):
                uses_long = (flank_ee, long_s) in junctions
                uses_short = (flank_ee, short_s) in junctions
            elif flank_es >= max(long_e, short_e):
                uses_long = (long_e, flank_es) in junctions
                uses_short = (short_e, flank_es) in junctions
            else:
                uses_long = uses_short = None
        else:
            uses_long = uses_short = None

        def _exon_equal(lo, hi):
            return bool(tx_exons) and (lo, hi) in set(tx_exons)

        def _overlap_len(rs, re_):
            tot = 0
            for cs, ce, _ in ref_sorted:
                o_s, o_e = max(cs, rs), min(ce, re_)
                if o_s < o_e:
                    tot += (o_e - o_s)
            return tot

        if uses_long is not None and uses_long != uses_short:
            ref_is_long = uses_long
        elif _exon_equal(long_s, long_e) and not _exon_equal(short_s, short_e):
            ref_is_long = True
        elif _exon_equal(short_s, short_e) and not _exon_equal(long_s, long_e):
            ref_is_long = False
        else:
            ref_is_long = not (_overlap_len(short_s, short_e) > _overlap_len(long_s, long_e))

        # ---- 3) CDS 层面的长短替换 ----
        #  ref=长 → Alt=短: 删除/修剪 CDS 中与 V 重叠的段（长形独有的编码区被
        #                   短形跳过）；若 CDS 从 flanking 剪接位点连续读入而
        #                   修剪后不再有窗口，则补短形外显子窗口（读通至终止子）。
        #  ref=短 → Alt=长: 读通守卫成立时在 flanking 剪接位点一侧插入 V 段
        #                   （短形跳过的编码区在长形中恢复为编码）。
        new_cds = []
        for cs, ce, ph in ref_sorted:
            if ref_is_long and min(ce, v_e) > max(cs, v_s):
                if cs < v_s:
                    new_cds.append((cs, v_s, ph))
                if ce > v_e:
                    new_cds.append((v_e, ce, ph))
            else:
                new_cds.append((cs, ce, ph))

        # 读通守卫: CDS 是否触及 flanking 剪接连接点（编码框跨越该连接点）
        if flank_ee <= min(long_s, short_s):
            reads_through = any(ce == flank_ee for _, ce, _ in ref_sorted)
        else:
            reads_through = any(cs == flank_es for cs, _, _ in ref_sorted)

        alt_es, alt_ee = (short_s, short_e) if ref_is_long else (long_s, long_e)
        if reads_through:
            if strand == '+':
                before = sum(e - s for s, e, _ in new_cds if e <= alt_es)
            else:
                before = sum(e - s for s, e, _ in new_cds if s >= alt_ee)
            if ref_is_long:
                # 修剪后若无从共享边界读入的窗口，补短形外显子窗口
                already = any((cs == alt_es) if strand == '+' else (ce == alt_ee)
                              for cs, ce, _ in new_cds)
                if not already:
                    new_cds.append((alt_es, alt_ee, before % 3))
            else:
                # ref=短 → Alt=长: 插入 V 段（长形多出的编码区）
                if v_s != v_e:
                    new_cds.append((v_s, v_e, before % 3))

        if not new_cds:
            return ref_sorted

        new_cds = sorted(new_cds, key=lambda x: x[0])
        merged = []
        for s, e, ph in new_cds:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e), merged[-1][2])
            else:
                merged.append((s, e, ph))
        return merged

    return ref_sorted
