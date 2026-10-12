#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GTF index construction and disk caching.

Index structure:
  gene_idx   {canonical_gene: {transcript_id: [(cds_start_0, cds_end_0, frame), ...]}}
  gene_exons {canonical_gene: {transcript_id: [(exon_start_0, exon_end_0), ...]}}
  tid_strand {transcript_id: strand}
  alias_map  {alias(gene_name/gene_id/gene): canonical_gene}
  longest_tx_per_gene {canonical_gene: (transcript_id, cds_length)}
  tid_biotype {transcript_id: transcript_type}   # e.g. protein_coding / nonsense_mediated_decay
"""

# Bump when the cached tuple layout changes so stale caches are never loaded
_INDEX_FORMAT_VERSION = 2

import gzip
import hashlib
import logging
import os
import pickle
import re
from collections import defaultdict
from datetime import datetime

from . import settings

logger = logging.getLogger(__name__)


def get_gtf_cache_key(gtf_file: str) -> str:
    stat = os.stat(gtf_file)
    return hashlib.md5(f"{gtf_file}_{stat.st_mtime}_{stat.st_size}".encode()).hexdigest()


def build_gtf_index_with_cache(gtf_file: str, cache_dir: str = None):
    if cache_dir is None:
        cache_dir = settings.get_cache_dir()

    os.makedirs(cache_dir, exist_ok=True)
    cache_key = get_gtf_cache_key(gtf_file)
    cache_file = os.path.join(cache_dir, f'gtf_index_{cache_key}_v{_INDEX_FORMAT_VERSION}.pkl')
    temp_cache_file = cache_file + '.tmp'

    logger.debug(f"Cache directory: {cache_dir}")
    logger.debug(f"Cache file: {cache_file}")

    if os.path.exists(cache_file):
        file_size = os.path.getsize(cache_file)
        if file_size < 10240:
            logger.warning(f"Cache file too small ({file_size} bytes), removing ...")
            try:
                os.remove(cache_file)
                logger.info(f"Removed small cache file: {cache_file}")
            except Exception as e:
                logger.warning(f"Failed to remove small cache file: {e}")
        else:
            logger.info(f"Loading cached GTF index ({file_size/1024/1024:.1f} MB) ...")
            try:
                with open(cache_file, 'rb') as f:
                    gene_idx, gene_exons, tid_strand, alias_map, longest_tx_per_gene, tid_biotype = pickle.load(f)

                if not gene_idx or len(gene_idx) == 0:
                    raise ValueError("Loaded empty gene index")

                total_tx = sum(len(tx) for tx in gene_idx.values())
                total_exons = sum(sum(len(exons) for exons in tx.values()) for tx in gene_exons.values())
                logger.info(f"✓ GTF index: {len(gene_idx):,} genes · {total_tx:,} transcripts · "
                            f"{total_exons:,} exons (from cache)")

                return gene_idx, gene_exons, tid_strand, alias_map, longest_tx_per_gene, tid_biotype
            except (pickle.UnpicklingError, EOFError, ValueError, KeyError, AttributeError) as e:
                logger.warning(f"❌ Cache load failed: {e}, rebuilding...")
                try:
                    backup_file = cache_file + '.corrupted'
                    if os.path.exists(backup_file):
                        os.remove(backup_file)
                    os.rename(cache_file, backup_file)
                    logger.info(f"Corrupted cache backed up to: {backup_file}")
                except Exception as backup_err:
                    logger.warning(f"Failed to backup corrupted cache: {backup_err}")
                    try:
                        os.remove(cache_file)
                        logger.info(f"Removed corrupted cache file: {cache_file}")
                    except Exception as remove_err:
                        logger.warning(f"Failed to remove corrupted cache: {remove_err}")

    logger.info("Building GTF index (this may take a while) ...")
    start_time = datetime.now()
    gene_idx, gene_exons, tid_strand, alias_map, longest_tx_per_gene, tid_biotype = build_gtf_index(gtf_file)
    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"GTF index built in {elapsed:.1f}s")

    if not gene_idx or len(gene_idx) == 0:
        logger.error("Failed to build GTF index: empty result")
        return gene_idx, gene_exons, tid_strand, alias_map, longest_tx_per_gene, tid_biotype

    try:
        logger.debug(f"Saving GTF index to cache: {cache_file}")
        with open(temp_cache_file, 'wb') as f:
            pickle.dump((gene_idx, gene_exons, tid_strand, alias_map, longest_tx_per_gene, tid_biotype), f,
                        protocol=pickle.HIGHEST_PROTOCOL)

        test_size = os.path.getsize(temp_cache_file)
        if test_size < 10240:
            raise ValueError(f"Cache file too small: {test_size} bytes")

        if os.path.exists(cache_file):
            os.remove(cache_file)
        os.rename(temp_cache_file, cache_file)

        cache_size = os.path.getsize(cache_file) / (1024 * 1024)
        logger.info(f"✓ GTF index cached ({cache_size:.1f} MB): {cache_file}")
    except Exception as e:
        logger.warning(f"❌ Cache save failed: {e}")
        if os.path.exists(temp_cache_file):
            try:
                os.remove(temp_cache_file)
                logger.info(f"Removed temporary cache file: {temp_cache_file}")
            except Exception as cleanup_err:
                logger.warning(f"Failed to remove temporary file: {cleanup_err}")

    return gene_idx, gene_exons, tid_strand, alias_map, longest_tx_per_gene, tid_biotype


def _default_gene_dict():
    return defaultdict(list)


def build_gtf_index(gtf_file: str):
    gene_idx = defaultdict(_default_gene_dict)
    gene_exons = defaultdict(_default_gene_dict)
    tid_strand = {}
    tid_biotype = {}
    alias_map = {}
    tx_cds_len = defaultdict(int)

    if not gtf_file or not os.path.exists(gtf_file):
        logger.error(f"GTF file not found: {gtf_file}")
        return {}, {}, {}, {}, {}, {}

    logger.debug(f"Parsing GTF: {gtf_file}")
    opener = gzip.open if str(gtf_file).endswith('.gz') else open

    line_count = 0
    exon_count = 0
    cds_count = 0

    with opener(gtf_file, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            if len(p) < 9:
                continue
            ft = p[2]
            attrs = parse_gtf_attrs(p[8])
            tid = attrs.get('transcript_id', '')
            gene_name = attrs.get('gene_name', '')
            gene_id = attrs.get('gene_id', '')
            gene_alt = attrs.get('gene', '')

            canonical = gene_name or gene_alt or gene_id
            if not tid or not canonical:
                continue

            if ft in ('transcript', 'exon', 'CDS', 'stop_codon', 'start_codon', 'UTR'):
                bt = attrs.get('transcript_type', attrs.get('transcript_biotype', ''))
                if bt and tid not in tid_biotype:
                    tid_biotype[tid] = bt

            start_0 = int(p[3]) - 1
            end_0 = int(p[4])
            strand = p[6]
            tid_strand[tid] = strand

            for alias in [gene_name, gene_id, gene_alt]:
                if alias and alias != canonical:
                    alias_map[alias] = canonical

            if gene_name and gene_name != canonical:
                alias_map[gene_name] = canonical

            if ft == 'CDS':
                frame = int(p[7]) if p[7] in ('0', '1', '2') else 0
                gene_idx[canonical][tid].append((start_0, end_0, frame))
                tx_cds_len[tid] += (end_0 - start_0)
                cds_count += 1
            elif ft == 'exon':
                gene_exons[canonical][tid].append((start_0, end_0))
                exon_count += 1

            line_count += 1
            if line_count % 100000 == 0:
                logger.info(f"  · {line_count:,} GTF lines processed")

    logger.debug("Sorting CDS and exon coordinates ...")
    for canonical in gene_idx:
        for tid in gene_idx[canonical]:
            gene_idx[canonical][tid].sort(key=lambda x: x[0])
    for canonical in gene_exons:
        for tid in gene_exons[canonical]:
            gene_exons[canonical][tid].sort(key=lambda x: x[0])

    logger.debug("Identifying the longest transcript per gene ...")
    longest_tx_per_gene = {}
    for g in gene_idx:
        best_tid, best_len = None, 0
        for t in gene_idx[g]:
            clen = tx_cds_len.get(t, 0)
            if clen > best_len:
                best_len = clen
                best_tid = t
        if best_tid:
            longest_tx_per_gene[g] = (best_tid, best_len)

    gene_count = len(gene_idx)
    transcript_count = sum(len(tx) for tx in gene_idx.values())

    logger.info(f"✓ GTF index built: {gene_count:,} genes · {transcript_count:,} transcripts · "
                f"{exon_count:,} exons · {cds_count:,} CDS · {len(alias_map):,} aliases")

    gene_idx_dict = {}
    for gene, transcripts in gene_idx.items():
        gene_idx_dict[gene] = dict(transcripts)

    gene_exons_dict = {}
    for gene, transcripts in gene_exons.items():
        gene_exons_dict[gene] = dict(transcripts)

    return gene_idx_dict, gene_exons_dict, tid_strand, alias_map, longest_tx_per_gene, tid_biotype


def parse_gtf_attrs(attr_string: str) -> dict:
    d = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_string):
        d[m.group(1)] = m.group(2)
    for m in re.finditer(r'(\w+)\s+([^\s;"\']+)', attr_string):
        k, v = m.group(1), m.group(2).strip('"').strip("'")
        if k not in d and v:
            d[k] = v
    return d
