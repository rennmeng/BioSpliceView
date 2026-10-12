#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BioSpliceView — rMATS alternative splicing event annotation toolkit.

Split from the original single-file ``spliceproduce.py``:

  - settings.py        global settings & constants
  - gtf_index.py       GTF index construction & disk cache
  - sequences.py       genome sequence retrieval / CDS assembly / truncation
  - nmd.py             PTC genomic mapping and NMD classification
  - events.py          rMATS event coordinate parsing & schema constants
  - matching.py        reference transcript anchoring & adjacency matching
  - alt_isoform.py     alt-isoform CDS derivation from the reference CDS
  - genbank_output.py  GenBank (EXON/Tran) record generation
  - annotate.py        per-event annotation and the filter_and_summarize pipeline
  - cli.py             command-line entry point
"""

import warnings

# Translating isoform CDS whose length is not a multiple of three (frameshift
# events, stop-codon truncation, phase carry-over) triggers Biopython's
# "Partial codon" warning. This is expected behaviour in BioSpliceView - the
# partial codon is intentionally kept - so the warning is silenced here.
warnings.filterwarnings("ignore", message=".*Partial codon.*")

__version__ = "1.1.0"

from .annotate import filter_and_summarize, process_event  # noqa: F401

__all__ = ["filter_and_summarize", "process_event", "__version__"]
