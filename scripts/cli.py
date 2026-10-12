#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command-line entry point."""

import argparse
import logging
import sys

from . import __version__, console, settings
from .annotate import filter_and_summarize

logger = logging.getLogger(__name__)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog='biospliceview',
        description='BioSpliceView — annotate rMATS alternative splicing events '
                    '(adjacency-aware transcript matching)')
    p.add_argument('-r', '--rmats-dir', required=True, help='rMATS output directory')
    p.add_argument('-g', '--genome', required=True, help='genome FASTA file')
    p.add_argument('-t', '--gtf', required=True, help='GTF annotation file')
    p.add_argument('-o', '--output-dir', default='./results', help='output directory')
    p.add_argument('--fdr', type=float, default=0.05, help='FDR threshold')
    p.add_argument('--psi', type=float, default=0.2, help='delta-PSI threshold')
    p.add_argument('--minjc', type=float, default=15, help='minimum junction count')
    p.add_argument('--job-id', default=None, help='job ID (optional)')
    p.add_argument('--no-gb', action='store_true', help='do not generate GenBank files')
    p.add_argument('--no-cache', action='store_true', help='disable the GTF index cache (enabled by default)')
    p.add_argument('--cache-dir', default=None, help='GTF cache directory (default: see config)')
    return p.parse_args(argv)


def _setup_logging():
    handler = logging.StreamHandler()
    handler.setFormatter(console.ColorFormatter(use_color=console.ansi_enabled(sys.stderr)))
    logging.basicConfig(level=logging.INFO, handlers=[handler])


def main(argv=None):
    _setup_logging()
    args = parse_args(argv)

    cache_note = ('disabled' if args.no_cache
                  else f"on ({args.cache_dir or settings.GTF_CACHE_DIR})")
    print(console.banner(
        title=f"BioSpliceView v{__version__}",
        subtitle="rMATS alternative splicing annotation · adjacency-aware matching",
        rows=[
            ("rMATS dir", args.rmats_dir),
            ("Genome", args.genome),
            ("GTF", args.gtf),
            ("Output", f"{args.output_dir}  (job: {args.job_id or 'auto-generated'})"),
            ("Thresholds", f"FDR ≤ {args.fdr:g} · |ΔPSI| > {args.psi:g} · MinJC ≥ {args.minjc:g}"),
            ("Features", f"GenBank {'off' if args.no_gb else 'on'} · GTF cache {cache_note}"),
        ],
    ), flush=True)

    filter_and_summarize(
        args.rmats_dir,
        args.output_dir,
        args.genome,
        args.gtf,
        args.fdr,
        args.psi,
        args.minjc,
        args.job_id,
        not args.no_gb,
        not args.no_cache,
        args.cache_dir
    )


if __name__ == '__main__':
    main()
