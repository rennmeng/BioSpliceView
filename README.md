# BioSpliceView

<p align="center">
  <b>Annotation of rMATS alternative splicing events at the isoform level</b><br>
  Reference/alternative CDS &amp; protein reconstruction · frameshift detection · NMD classification · GenBank export
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/python-3.8%2B-blue">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-green">
</p>

**BioSpliceView** annotates significant alternative splicing events (SE / RI / A3SS / A5SS / MXE) detected by
[rMATS](http://rnaseq-mats.sourceforge.net/): it anchors every event onto the reference transcriptome,
reconstructs the reference- and alternative-isoform CDS/protein sequences, classifies NMD sensitivity, and
writes a summary Excel workbook plus per-event GenBank records.

## Quick start

```bash
git clone https://github.com/rennmeng/BioSpliceView.git
cd BioSpliceView
pip install -r requirements.txt    # pandas, openpyxl, pyfaidx, biopython (Python 3.8+)
python example.py                  # runs the bundled mini dataset end-to-end in seconds
```

Results are written to `results_example/Job_mini/`: `Job_mini_results.xlsx` plus two GenBank records per event
under `GB_Files/`. Optionally, `pip install .` provides the `biospliceview` console command.

## Usage

```bash
python main.py \
    -r  <rmats-output-dir> \
    -g  <genome.fa> \
    -t  <annotation.gtf> \
    -o  ./results \
    --job-id myjob \
    --fdr 0.05 --psi 0.2 --minjc 15
```

| Parameter | Description |
|-----------|-------------|
| `-r` | Directory containing the rMATS analysis output files: `SE.MATS.JC.txt`, `A3SS.MATS.JC.txt`, `A5SS.MATS.JC.txt`, `RI.MATS.JC.txt`, `MXE.MATS.JC.txt` (missing event types are skipped) |
| `-g` | Genome FASTA file (a `.fai` index is created on first use). **Must be the same reference genome as used in the rMATS analysis** |
| `-t` | GTF annotation file (requires `gene_name` / `transcript_id` attributes). **Must be the same reference annotation as used in the rMATS analysis** |
| `-o` | Output directory |
| `--job-id` | Job ID used in the output folder/file names; auto-generated if omitted |
| `--fdr` | FDR threshold (default `0.05`) — keep events with `FDR <= fdr` |
| `--psi` | ΔPSI threshold (default `0.2`) — keep events with `abs(IncLevelDifference) > psi` |
| `--minjc` | Minimum junction count (default `15`) — exclude events whose IJC **and** SJC averages are too low in both samples |

Optional flags: `--no-gb` (skip GenBank output) and `--no-cache` / `--cache-dir` (GTF index cache control).
The default cache directory can also be set via the environment variable `SPLICE_GTF_CACHE_DIR`.

## Output

`<output>/Job_<id>/Job_<id>_results.xlsx` contains one sheet per event type plus an `All_Significant` summary.
The raw rMATS columns are kept and 14 annotation columns are appended:

| Column | Meaning |
|--------|---------|
| `Ref_Transcript_ID` | anchored reference transcript |
| `Anchor_Biotype` | annotation biotype of the anchor (e.g. `protein_coding`, `nonsense_mediated_decay`) |
| `Ref_Full_CDS` / `Ref_Full_Protein` / `Ref_CDS_Length` | reference isoform, truncated at the first stop codon found after extending into the 3'UTR |
| `Alt_Full_CDS` / `Alt_Full_Protein` / `Alt_CDS_Length` | alternative isoform; empty and flagged `ALT_NO_START_CODON` if it does not start on ATG |
| `Diff_CDS` / `Diff_Protein` / `Diff_Length` | difference segments after trimming common prefixes/suffixes (`ref\|alt` format) |
| `NMD_Status` | **event consequence** — whether the event itself introduces a PTC: `Express` (frame-preserving / protein identical) · `NMD_Target` · `NMDER1/2/3` (PTC predicted to escape NMD) · `NO_START_CODON` · `NON_CODING` |
| `Ref_Absolute_NMD` / `Alt_Absolute_NMD` | **absolute NMD status of each isoform's own stop codon** (`Normal` / `NMD_Target` / `NMDER1/2/3`), independent of the event — reveals isoforms that are already NMD substrates in the annotation |
| `Match_Status` | `FULL` / `PARTIAL(...)` / `ANCHOR_ONLY(...)` / `NON_CODING*` / `UNANCHORED(...)`; anchor warnings are appended as `;ANCHOR_SELF_NMD` / `;ANCHOR_NMD_BIOTYPE` |

### Dual NMD reporting

`NMD_Status` answers *"does this event introduce a premature stop codon relative to the
reference isoform?"* — the question an experimenter usually asks. `Ref_Absolute_NMD` /
`Alt_Absolute_NMD` answer *"is this isoform itself an NMD substrate according to its own
stop-codon position?"*. The two can legitimately disagree: e.g. an in-frame event on an
NMD-annotated reference yields `NMD_Status = Express` (no new PTC) while
`Alt_Absolute_NMD = NMD_Target` (the isoform, like its reference, remains NMD-sensitive).

NMD escape subclasses follow the manuscript definitions:
`NMDER1` PTC < 150 nt downstream of the start codon · `NMDER2` PTC > 400 nt from the next
downstream exon-junction complex · `NMDER3` PTC within the last exon or ≤ 50 nt upstream of
the terminal exon–exon junction. All distances are measured in transcript orientation
(strand-aware).

### Anchor preference

Within the same match level, candidate anchors are preferred in the order
(1) not NMD-annotated and stop codon not NMD-sensitive, (2) not NMD-annotated but
self-NMD, (3) NMD-annotated (kept only when no safer candidate realises the event and
flagged `;ANCHOR_NMD_BIOTYPE`). This keeps BioSpliceView's assumption — the reference
transcript is normally expressed — explicit and auditable.

## How events are matched

- **SE**: upstream → skipped → downstream exons must be adjacent in order
- **RI**: upstream and downstream exons adjacent, retained intron cleanly preserved
- **MXE**: both flanks present, only one mutually exclusive exon, adjacency preferred
- **A5SS/A3SS**: short/long exon + flanking exon boundary matching

Every event is graded `FULL` / `PARTIAL` / `ANCHOR_ONLY` / `NON_CODING` / `UNANCHORED`, and frameshift events are
further classified for NMD sensitivity by the canonical EJC rules (`NMD_Target` / `NMDER1/2/3`).

## GenBank output

Per event, two records are written under `GB_Files/`: `EXON-*.gb` (spliced mRNA view) and
`Tran-*.gb` (genomic view with introns). Records carry exon/CDS features (`transl_table=1`),
`misc_feature` landmarks using the manuscript vocabulary — anchor exons, skipped exon,
retained intron, mutually exclusive exons, alternative splice site region — and a structured
`COMMENT` block with the dual NMD report, directly openable in SnapGene or any GenBank viewer.

## Bundled mini dataset

`examples/mini_dataset/` contains only what the 262 events of a real annotation job involve — 228 transcripts on
20 chromosomes, shrinking the 2.7 GB mm10 genome to 15 MB — with all coordinates consistently shifted so the
dataset is fully self-contained:

```
examples/mini_dataset/
├── mm10_mini.fa (+ .fai)   # 216 genomic windows (transcript spans / event clusters ± 2000 bp)
├── refGene_mini.gtf        # transcript/exon/CDS lines of those transcripts
├── regions.tsv             # window mapping: new = orig - shift (0-based half-open)
└── rmats/                  # rebuilt rMATS inputs (coordinates shifted too)
```

Rerunning it reproduces the original full-genome annotation exactly: every matched event is identical on all
annotation columns; the only difference affects 4 events that were already unmatched (`UNANCHORED`) in the
original job, which report `GENE_NOT_FOUND` instead of `PRECISE_MATCH_FAIL`.

## Project layout

```
BioSpliceView/
├── main.py                # CLI entry point
├── example.py             # one-command demo on the bundled mini dataset
├── scripts/               # core package
│   ├── settings.py        # defaults & constants
│   ├── gtf_index.py       # GTF index construction + disk cache
│   ├── sequences.py       # sequence retrieval / CDS assembly / stop-codon truncation
│   ├── events.py          # event coordinate parsing, frameshift detection
│   ├── matching.py        # transcript anchoring & adjacency matching
│   ├── alt_isoform.py     # alt-isoform CDS derivation
│   ├── nmd.py             # PTC genomic mapping & NMD classification
│   ├── genbank_output.py  # EXON/Tran GenBank records
│   ├── annotate.py        # per-event annotation + main pipeline
│   └── cli.py             # argparse CLI
└── examples/mini_dataset/
```

Event coordinates follow the rMATS convention (`*Start_0base`/`*ES` columns are 0-based, `*End`/`*EE` columns
are 1-based); the GTF index works with 0-based half-open intervals internally.
