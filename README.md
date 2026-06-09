# bait2contig

`bait2contig` is a lightweight command-line tool for finding contigs that match bait or reference sequences under user-defined identity and bait-coverage thresholds. It maps bait sequences to a contig FASTA with `minimap2`, parses PAF output, writes a hit-level TSV, and can summarize contigs anchored by each bait.

Typical bait sequences include full-length 16S rRNA sequences, ITS sequences, marker genes, MAG markers, viral markers, and custom reference sequences.

## Installation

Python 3.9 or newer is required. The package has no required Python dependencies. If `psutil` is installed, `bait2contig` uses it for more detailed CPU and memory monitoring; otherwise it falls back to the Python standard library.

Install the latest published package from PyPI:

```bash
python -m pip install bait2contig
```

Install with optional detailed resource monitoring support:

```bash
python -m pip install "bait2contig[monitor]"
```

Upgrade an existing installation:

```bash
python -m pip install --upgrade bait2contig
```

Force reinstall and upgrade, which is useful when the local environment may be using an old cached or editable copy:

```bash
python -m pip install --upgrade --force-reinstall --no-cache-dir bait2contig
```

Install from the current project directory:

```bash
python -m pip install .
```

Force reinstall from the current project directory:

```bash
python -m pip install --force-reinstall --no-cache-dir .
```

For editable development from the current project directory:

```bash
python -m pip install -e .
```

If the `bait2contig` command still resolves to an older installation, inspect the command path and package location:

```bash
which bait2contig
python -m pip show bait2contig
```

## minimap2 Dependency

`bait2contig search` requires `minimap2` on `PATH`, or a path supplied with `--minimap2`.

```bash
bait2contig search --minimap2 /path/to/minimap2 --contigs contigs.fa --bait bait.fa --out bait2contig.hits.tsv
```

The internal mapping command is:

```bash
minimap2 -x {preset} -t {threads} {contigs} {bait} > {tmp_paf}
```

The order is contigs first and bait second, so PAF query IDs are bait IDs and PAF target IDs are contig IDs.

## Quick Start

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32

bait2contig summarize \
  --hits bait2contig.hits.tsv \
  --out bait2contig.summary.tsv \
  --best-hit \
  --include-contigs
```

## Command Overview

`bait2contig` has two subcommands:

| Command | Purpose | Required inputs | Main output |
| --- | --- | --- | --- |
| `bait2contig search` | Map bait/reference sequences to contigs and write hit-level TSV output. | `--contigs`, `--bait`, `--out` | Hit TSV, optional matched-contig FASTA, optional kept PAF |
| `bait2contig summarize` | Summarize hit-level TSV output by bait. | `--hits`, `--out` | Bait-level summary TSV |

Top-level help:

```bash
bait2contig --help
bait2contig search --help
bait2contig summarize --help
```

## External Invocation

External programs should call `bait2contig` as a command-line program. The Python functions inside `bait2contig` are implementation details and are not a stable public API.

Recommended Python subprocess usage:

```python
import subprocess

cmd = [
    "bait2contig",
    "search",
    "--contigs", "contigs.fa",
    "--bait", "bait.fa",
    "--out", "bait2contig.hits.tsv",
    "--identity", "0.97",
    "--coverage", "0.80",
    "--threads", "8",
    "--quiet",
]

result = subprocess.run(cmd, text=True, capture_output=True)
if result.returncode != 0:
    raise RuntimeError(result.stderr)
```

Recommended summarize call:

```python
import subprocess

cmd = [
    "bait2contig",
    "summarize",
    "--hits", "bait2contig.hits.tsv",
    "--out", "bait2contig.summary.tsv",
    "--best-hit",
    "--include-contigs",
    "--quiet",
]

result = subprocess.run(cmd, text=True, capture_output=True)
if result.returncode != 0:
    raise RuntimeError(result.stderr)
```

External caller contract:

| Item | Behavior |
| --- | --- |
| Success exit code | `0` |
| Runtime or validation failure | `1` |
| Command-line parsing failure | `2` |
| Interrupted process | `130` |
| Machine-readable data | TSV and FASTA output files |
| Progress and error text | stderr and plain-text log file |
| stdout | Reserved for help output; do not rely on stdout during normal workflows |
| Logs | Plain text, never gzip-compressed |
| Help color | Disable with `--no-color` or `NO_COLOR=1` when capturing help text |

Use argument lists in external programs instead of shell strings. This avoids quoting errors with paths that contain spaces.

For deterministic external integration:

- Use absolute paths for inputs and outputs.
- Pass `--quiet` to reduce terminal noise while retaining warnings and errors.
- Pass `--no-color` if capturing help output.
- Use `--gzip` only when callers expect `.gz` output paths.
- Check the actual output path after applying gzip rules.
- Check exit code before reading output files.
- Parse TSV output, not human-readable log messages.
- Parse log marker blocks only if resume status needs to be audited.

## Search Parameters

Required arguments:

| Option | Type | Required | Description |
| --- | --- | --- | --- |
| `--contigs FILE` | path | yes | Input contig FASTA. Plain and `.gz` files are supported. |
| `--bait FILE` | path | yes | Input bait/reference FASTA. Plain and `.gz` files are supported. |
| `--out FILE` | path | yes | Output hit TSV path. With `--gzip`, `.gz` is appended if not already present. |

Filtering arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--identity FLOAT` | 0 to 1 | `0.97` | Minimum alignment identity. |
| `--coverage FLOAT` | 0 to 1 | `0.80` | Minimum bait coverage. |
| `--min-aln-length INT` | integer >= 0 | `0` | Minimum alignment block length. |
| `--terminal-tolerance INT` | integer >= 0 | `5` | Allowed unaligned bases at sequence ends for partial bait alignments. |
| `--no-terminal-filter` | flag | off | Disable terminal-placement filtering for partial bait hits. |
| `--best-only` | flag | off | Keep only one best contig per bait in the hit TSV. |

Annotation arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--lineage FILE` | path | not set | Optional bait lineage TSV with `bait_id` and `lineage` columns. |
| `--circular-list FILE` | path | not set | Optional list of circular contig IDs, one ID per line. Overrides FASTA header circularity inference. |

Mapping arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--preset STR` | string | `asm10` | minimap2 preset passed with `-x`. |
| `--threads INT` | integer >= 1 | `8` | Number of minimap2 threads. |
| `--minimap2 PATH` | path or executable name | `minimap2` | minimap2 executable. |
| `--keep-paf` | flag | off | Keep intermediate PAF output next to the hit TSV. |
| `--tmp-dir DIR` | path | output directory | Directory for the temporary plain-text PAF file. |

Contig extraction arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--extract-contigs FILE` | path | not set | Write matched contig sequences to FASTA. Extraction runs only when this is provided. |
| `--extract-mode STR` | choice | `all` | One of `all`, `best`, `circular`, or `non-circular`. |
| `--extract-min-identity FLOAT` | 0 to 1 | `--identity` | Minimum identity for extracted contigs. |
| `--extract-min-coverage FLOAT` | 0 to 1 | `--coverage` | Minimum bait coverage for extracted contigs. |
| `--extract-min-aln-length INT` | integer >= 0 | `--min-aln-length` | Minimum alignment length for extracted contigs. |
| `--extract-rename` | flag | off | Rename extracted FASTA headers to include bait and hit metrics. |
| `--extract-include-lineage` | flag | off | Include lineage in renamed FASTA headers. Requires `--extract-rename`. |
| `--extract-dedup` | flag | on | Deduplicate extracted contigs by `ctg_id`. |
| `--no-extract-dedup` | flag | off | Allow repeated contig sequences when a contig is matched by multiple bait sequences. |

Resume and output arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--resume` | flag | off | Skip work only if the log contains a matching successful DONE marker and the output is valid. |
| `--rerun` | flag | off | Force recomputation and overwrite output. Cannot be used with `--resume`. |
| `--force` | flag | off | Allow overwriting existing output without resume checks. |
| `--gzip` | flag | off | Compress TSV, kept PAF, and extracted FASTA outputs when applicable. |
| `--log FILE` | path | `<actual_out>.log` | Plain-text log path. |

Runtime and logging arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--monitor-interval INT` | integer >= 1 | `30` | Seconds between resource log lines. |
| `--no-color` | flag | off | Disable colored terminal output and help. |
| `--quiet` | flag | off | Show only warnings and errors on screen. |
| `--verbose` | flag | off | Show detailed logs on screen. |

## Summarize Parameters

Required arguments:

| Option | Type | Required | Description |
| --- | --- | --- | --- |
| `--hits FILE` | path | yes | Hit TSV produced by `bait2contig search`. Plain and `.gz` files are supported. |
| `--out FILE` | path | yes | Output summary TSV. With `--gzip`, `.gz` is appended if not already present. |

Filtering arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--min-identity FLOAT` | 0 to 1 | not set | Additional identity filter applied to input hits. |
| `--min-coverage FLOAT` | 0 to 1 | not set | Additional bait coverage filter applied to input hits. |
| `--min-aln-length INT` | integer >= 0 | not set | Additional alignment length filter applied to input hits. |

Summary arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--best-hit` | flag | off | Add best-contig columns for each bait. |
| `--include-contigs` | flag | off | Add ordered contig ID list for each bait. |
| `--contig-sep STR` | string | `,` | Separator used for the contig list. |

Resume and output arguments:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--resume` | flag | off | Skip work only if the log contains a matching successful DONE marker and output is valid. |
| `--rerun` | flag | off | Force recomputation and overwrite output. Cannot be used with `--resume`. |
| `--force` | flag | off | Allow overwriting existing output without resume checks. |
| `--gzip` | flag | off | Compress summary TSV output. |
| `--log FILE` | path | `<actual_out>.log` | Plain-text log path. |

Runtime and logging arguments are the same as `search`: `--monitor-interval`, `--no-color`, `--quiet`, and `--verbose`.

## Output Path Rules

`bait2contig` distinguishes requested output paths from actual output paths.

| Request | Actual output |
| --- | --- |
| `--out bait2contig.hits.tsv` | `bait2contig.hits.tsv` |
| `--out bait2contig.hits.tsv --gzip` | `bait2contig.hits.tsv.gz` |
| `--out bait2contig.hits.tsv.gz --gzip` | `bait2contig.hits.tsv.gz` |
| `--extract-contigs matched.fa --gzip` | `matched.fa.gz` |
| `--extract-contigs matched.fa.gz --gzip` | `matched.fa.gz` |

Default log path is always based on the actual output path:

```text
<actual_out>.log
```

Examples:

```text
bait2contig.hits.tsv.log
bait2contig.hits.tsv.gz.log
bait2contig.summary.tsv.log
bait2contig.summary.tsv.gz.log
```

Logs are plain text even when `--gzip` is used.

## Search

Basic search:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32
```

Search with lineage:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --lineage bait.lineage.tsv \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32
```

Search with gzip output:

```bash
bait2contig search \
  --contigs contigs.fa.gz \
  --bait bait.fa.gz \
  --out bait2contig.hits.tsv \
  --gzip \
  --threads 32
```

Actual output:

```text
bait2contig.hits.tsv.gz
```

If minimap2 is not found, `bait2contig` exits with:

```text
ERROR: minimap2 was not found. Please install minimap2 or provide its path with --minimap2.
```

## Summarize

```bash
bait2contig summarize \
  --hits bait2contig.hits.tsv \
  --out bait2contig.summary.tsv \
  --best-hit \
  --include-contigs
```

With gzip output:

```bash
bait2contig summarize \
  --hits bait2contig.hits.tsv.gz \
  --out bait2contig.summary.tsv \
  --gzip \
  --best-hit \
  --include-contigs
```

Actual output:

```text
bait2contig.summary.tsv.gz
```

## Input FASTA Format

Plain and gzip-compressed FASTA are supported:

```text
>contig_001 len=15320 circular=true
ATGC...
>contig_002 len=8401
ATGC...
```

The sequence ID is the first token before whitespace. For `>contig_001 len=15320 circular=true`, the contig ID is `contig_001`. Contig lengths are calculated from sequence length.

## Lineage TSV Format

The optional lineage file is a plain or gzip TSV with two columns:

```text
bait_id	lineage
Pace_16S_001	d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__;g__;s__
Pace_16S_002	d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__;g__;s__
```

Headerless input is also accepted. Duplicate `bait_id` values are rejected.

## Circular List Format

The optional circular list is a plain or gzip text file with one contig ID per line:

```text
contig_001
contig_008
contig_109
```

`is_circular` is inferred only from FASTA headers or `--circular-list`. It is not de novo circularity detection. Header keywords are case-insensitive and include `circular=true`, `is_circular=true`, `circular`, and `circ=true`.

## Search Output TSV

Without lineage:

```text
ctg_id	bait_id	identity	aln_length	cov_bait	ctg_len	is_circular
```

With lineage:

```text
ctg_id	bait_id	identity	aln_length	cov_bait	ctg_len	is_circular	lineage
contig_00012	Pace_16S_001	0.991234	1450	0.982394	32781	False	d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__;g__;s__
```

Column definitions:

| Column | Type | Description |
| --- | --- | --- |
| `ctg_id` | string | Contig ID, taken from the first FASTA header token. |
| `bait_id` | string | Bait/reference ID, taken from the first FASTA header token. |
| `identity` | float | `residue_matches / alignment_block_length`, formatted with six decimals. |
| `aln_length` | integer | PAF alignment block length. |
| `cov_bait` | float | Aligned bait length divided by bait length, formatted with six decimals. |
| `ctg_len` | integer | Contig sequence length from FASTA when available. |
| `is_circular` | boolean text | `True` or `False`. Inferred from `--circular-list` or contig FASTA headers. |
| `lineage` | string | Optional lineage copied from `--lineage`; present only when lineage input is provided. |

If no alignments pass filters, the TSV is still written with only the header.

## Summary Output TSV

Default columns:

```text
bait_id	contig_count	circular_contig_count	total_ctg_len	mean_ctg_len	max_ctg_len	mean_identity	max_identity	mean_cov_bait	max_cov_bait	mean_aln_length	max_aln_length
```

If lineage is present in the hits TSV, `lineage` is preserved after `bait_id`. `--best-hit` adds best-hit columns, and `--include-contigs` adds a contig list.

Default summary column definitions:

| Column | Type | Description |
| --- | --- | --- |
| `bait_id` | string | Bait/reference ID. |
| `lineage` | string | Optional lineage, present only when the input hits TSV contains lineage. |
| `contig_count` | integer | Number of unique contigs anchored by the bait. |
| `circular_contig_count` | integer | Number of unique contigs marked circular. |
| `total_ctg_len` | integer | Sum of unique contig lengths. |
| `mean_ctg_len` | float | Mean unique contig length, formatted with two decimals. |
| `max_ctg_len` | integer | Maximum unique contig length. |
| `mean_identity` | float | Mean identity across retained best bait-contig hits, six decimals. |
| `max_identity` | float | Maximum identity across retained best bait-contig hits, six decimals. |
| `mean_cov_bait` | float | Mean bait coverage across retained best bait-contig hits, six decimals. |
| `max_cov_bait` | float | Maximum bait coverage across retained best bait-contig hits, six decimals. |
| `mean_aln_length` | float | Mean alignment length, formatted with two decimals. |
| `max_aln_length` | integer | Maximum alignment length. |

Additional `--best-hit` columns:

| Column | Type | Description |
| --- | --- | --- |
| `best_ctg_id` | string | Best contig for the bait by ranking rules. |
| `best_identity` | float | Best-hit identity, six decimals. |
| `best_cov_bait` | float | Best-hit bait coverage, six decimals. |
| `best_aln_length` | integer | Best-hit alignment length. |
| `best_ctg_len` | integer | Best-hit contig length. |
| `best_is_circular` | boolean text | `True` or `False` for the best hit. |

Additional `--include-contigs` column:

| Column | Type | Description |
| --- | --- | --- |
| `contigs` | string | Unique contig IDs ordered by best-hit ranking and joined with `--contig-sep`. |

## Identity and Coverage

`bait2contig` parses minimap2 PAF fields and calculates:

```text
identity = residue_matches / alignment_block_length
cov_bait = aligned bait length / bait length
```

More explicitly:

```text
cov_bait = (query_end - query_start) / query_len
```

Hits are filtered with `identity >= --identity`, `cov_bait >= --coverage`, and `aln_length >= --min-aln-length`.

By default, partial bait alignments are also filtered by end placement. If `cov_bait < 1.0`, the alignment must touch at least one bait end and one contig end within `--terminal-tolerance` bases. The default tolerance is 5 bp, which allows a few terminal bases to remain unaligned. Use `--no-terminal-filter` to disable this check.

Best-hit ranking is:

1. identity descending
2. cov_bait descending
3. aln_length descending
4. ctg_len descending

## Contig Extraction

Extract all matched contigs:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --extract-contigs matched_contigs.fa
```

Extract the best contig per bait:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --extract-contigs best_contigs.fa \
  --extract-mode best
```

Extract circular contigs only:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --extract-contigs circular_contigs.fa \
  --extract-mode circular
```

Use stricter extraction thresholds:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.95 \
  --coverage 0.70 \
  --extract-contigs high_confidence_contigs.fa \
  --extract-min-identity 0.99 \
  --extract-min-coverage 0.90
```

Use renamed extracted headers:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --lineage bait.lineage.tsv \
  --out bait2contig.hits.tsv \
  --extract-contigs matched_contigs.fa \
  --extract-rename \
  --extract-include-lineage
```

Renamed headers use:

```text
>{ctg_id}|bait={bait_id}|identity={identity:.6f}|cov_bait={cov_bait:.6f}|aln_length={aln_length}|ctg_len={ctg_len}|circular={is_circular}
```

## Gzip Examples

When `--gzip` is used, `.gz` is appended only if needed:

```bash
bait2contig search \
  --contigs contigs.fa.gz \
  --bait bait.fa.gz \
  --lineage bait.lineage.tsv.gz \
  --out bait2contig.hits.tsv \
  --gzip \
  --extract-contigs matched_contigs.fa \
  --extract-mode best \
  --threads 32
```

Actual outputs:

```text
bait2contig.hits.tsv.gz
matched_contigs.fa.gz
bait2contig.hits.tsv.gz.log
```

Log files are not gzip-compressed even when `--gzip` is used.

## Log-Based Resume

`--resume` relies on the plain-text log file and standardized `[BAIT2CONTIG_DONE]` markers. It does not only check whether output files exist.

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32 \
  --resume
```

If a successful matching run is found, `bait2contig` verifies the output path, output size, command, and key parameters before skipping. If parameters changed, the run is recomputed.

Force recomputation:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32 \
  --rerun
```

## CPU and Memory Monitoring

The log records periodic resource lines:

```text
[RESOURCE] elapsed=30.0s rss_mb=842.5 cpu_percent=1250.4
```

CPU percent may exceed 100 because minimap2 can use multiple threads. If `psutil` is unavailable, `bait2contig` still runs with standard-library fallback monitoring.

## Help Colors

`bait2contig --help`, `bait2contig search --help`, and `bait2contig summarize --help` use ANSI colors when stdout is a terminal. Pass `--no-color` to disable help colors. Set `CLICOLOR_FORCE=1` to force colors when output is captured; otherwise `NO_COLOR=1` disables colors.

## Command-Line Error Handling

`bait2contig` validates command-line input before starting a workflow. Missing commands, missing required arguments, unknown commands, misspelled options, invalid option values, missing input files, and output paths that point to directories are reported with concise `ERROR:` messages. When possible, the CLI suggests the nearest valid command or option and prints a minimal example command.

## Minimal Example Dataset

`contigs.fa`:

```text
>contig_001 len=12 circular=true
AAACCCGGGTTT
>contig_002 len=12
TTTGGGCCCAAA
```

`bait.fa`:

```text
>bait_001
AAACCCGGGTTT
```

`lineage.tsv`:

```text
bait_id	lineage
bait_001	d__Example;p__Example;c__;o__;f__;g__;s__
```

Expected hit TSV shape:

```text
ctg_id	bait_id	identity	aln_length	cov_bait	ctg_len	is_circular	lineage
contig_001	bait_001	1.000000	12	1.000000	12	True	d__Example;p__Example;c__;o__;f__;g__;s__
```

Expected summary TSV shape with `--best-hit --include-contigs`:

```text
bait_id	lineage	contig_count	circular_contig_count	total_ctg_len	mean_ctg_len	max_ctg_len	mean_identity	max_identity	mean_cov_bait	max_cov_bait	mean_aln_length	max_aln_length	best_ctg_id	best_identity	best_cov_bait	best_aln_length	best_ctg_len	best_is_circular	contigs
bait_001	d__Example;p__Example;c__;o__;f__;g__;s__	1	1	12	12.00	12	1.000000	1.000000	1.000000	1.000000	12.00	12	contig_001	1.000000	1.000000	12	12	True	contig_001
```

## Notes and Limitations

`bait2contig` does not perform de novo circularity detection. It does not parse SAM or BAM files, split taxonomy ranks, build databases, provide multiple aligner backends, or include a web interface. The search backend is minimap2 PAF output.

## License

`bait2contig` is released under the MIT License. See [LICENSE](LICENSE) for the full license text.
