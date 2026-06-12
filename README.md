# 🧬 bait2contig

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A lightweight, high-performance command-line tool for finding contigs that match bait or reference sequences under user-defined identity and bait-coverage thresholds. `bait2contig` maps bait sequences to contig FASTAs using `minimap2`, parses PAF output, writes hit-level TSV reports, and can summarize contigs anchored by each bait.

**Typical applications:** Full-length 16S rRNA sequences, ITS sequences, marker genes, MAG markers, viral markers, and custom reference sequences.

---

## 📖 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Overview](#command-overview)
- [Workflow](#workflow)
- [External Invocation](#external-invocation)
- [Parameters](#parameters)
- [Examples](#examples)
- [Output Formats](#output-formats)
- [Advanced Features](#advanced-features)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [License](#license)

---

## ⚙️ Installation

**Requirements:** Python 3.9 or newer

Python dependencies installed automatically:
- `rich` and `rich-argparse` – For readable help, progress, and terminal logs.
- `psutil` – For detailed CPU and memory monitoring. A standard-library fallback is still used if `psutil` is unavailable.

### From GitHub (Recommended)

```bash
pip install git+https://github.com/ypchan/bait2contig.git
```

**Force reinstall and upgrade** (useful for updating cached versions):

```bash
pip install --force-reinstall --no-cache-dir git+https://github.com/ypchan/bait2contig.git
```

### From GitHub Clone

If direct `pip install git+https://...` fails because of network, SSL, or pip Git backend issues, clone the repository first and install from the local checkout:

```bash
gh repo clone ypchan/bait2contig
cd bait2contig
python -m pip install -U .
```

Update an existing clone:

```bash
cd bait2contig
git pull
python -m pip install -U .
```

### From Local Directory

**Standard installation:**

```bash
python -m pip install .
```

**Force reinstall:**

```bash
python -m pip install --force-reinstall --no-cache-dir .
```

**Editable development mode:**

```bash
python -m pip install -e .
```

### Verify Installation

If `bait2contig` still resolves to an older installation, inspect the command path and package:

```bash
which bait2contig
python -m pip show bait2contig
```

### minimap2 Dependency

`bait2contig search` requires `minimap2` on `PATH`, or supply the path with `--minimap2`:

```bash
bait2contig search --minimap2 /path/to/minimap2 --contigs contigs.fa --bait bait.fa --out bait2contig.hits.tsv
```

**Mapping command used internally:**

```bash
minimap2 -x {preset} -t {threads} {contigs} {bait} > {tmp_paf}
```

> **Note:** Query IDs are bait IDs; target IDs are contig IDs.

### System Requirements

| Requirement | Details |
|-------------|---------|
| **Python** | 3.9+ |
| **minimap2** | 0.13+ (for `search` command) |
| **Disk space** | ~3x input FASTA size (for temporary PAF) |
| **RAM** | Typically 500MB–4GB depending on contig size and thread count |
| **CPU cores** | Scales with `--threads` (default 8) |

## 🚀 Quick Start

### Search for matching contigs:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32
```

### Summarize hits by bait:

```bash
bait2contig summarize \
  --hits bait2contig.hits.tsv \
  --out bait2contig.summary.tsv \
  --best-hit \
  --include-contigs
```

## 📋 Command Overview

`bait2contig` provides two primary subcommands:

| Command | Purpose | Inputs | Output |
|---------|---------|--------|--------|
| **`search`** | Map bait/reference sequences to contigs and write hit-level TSV | `--contigs`, `--bait`, `--out` | Hit TSV, matched FASTA*, kept PAF* |
| **`summarize`** | Summarize hit TSV output by bait | `--hits`, `--out` | Bait-level summary TSV |

*optional

### Get Help

```bash
bait2contig --help
bait2contig search --help
bait2contig summarize --help
```

---

## 🔄 Workflow

Typical use case involves two steps:

### Step 1: Search (Required)

Run `search` to map baits/references to your contigs and generate hit-level results:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait references.fa \
  --out hits.tsv \
  --identity 0.97 \
  --coverage 0.80
```

**Output:** `hits.tsv` contains one row per bait-contig match.

**Use this to:**
- Find all contigs matching your references
- Filter by strict identity/coverage thresholds
- Extract matched contig sequences
- Annotate hits with lineage information

### Step 2: Summarize (Optional)

Run `summarize` to aggregate results by bait, computing statistics and identifying best hits:

```bash
bait2contig summarize \
  --hits hits.tsv \
  --out summary.tsv \
  --best-hit \
  --include-contigs
```

**Output:** `summary.tsv` contains one row per bait with aggregate statistics.

**Use this to:**
- Get summary statistics per bait (mean/max identity, coverage, contig count)
- Identify the best contig match for each bait
- List all contigs matched by each bait in a single column

### Common Scenarios

| Scenario | Key Options |
|----------|------------|
| Find any matching contigs | Basic `search` (default) |
| Extract matched sequences | `search` + `--extract-contigs` |
| Find best hit per bait | `summarize` + `--best-hit` |
| Circular genome discovery | `search` + `--extract-mode circular` |
| Re-search with different thresholds | `search` + `--resume` (skips if already done) |
| High-confidence matching | `search` with `--identity 0.99 --coverage 0.95` |

## 🔗 External Invocation

External programs should call `bait2contig` as a command-line tool. **Python functions inside `bait2contig` are implementation details and not a stable public API.**

### Using Python subprocess

**Search command:**

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

**Summarize command:**

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

> **Tip:** Use argument lists instead of shell strings to avoid quoting errors with paths containing spaces.

### Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Runtime or validation failure |
| `2` | Command-line parsing failure |
| `130` | Interrupted process |

### Output Contract

| Item | Behavior |
|------|----------|
| **Machine-readable data** | TSV and FASTA output files |
| **Progress & errors** | stderr and plain-text log file |
| **stdout** | Reserved for help output; do not rely on stdout during normal workflows |
| **Logs** | Plain text, never gzip-compressed |
| **Help colors** | Disable with `--no-color` or `NO_COLOR=1` |

### Best Practices for External Integration

- Use **absolute paths** for inputs and outputs
- Pass `--quiet` to reduce terminal noise (warnings and errors preserved)
- Pass `--no-color` if capturing help output
- Use `--gzip` only when callers expect `.gz` output
- Check the **actual output path** after applying gzip rules
- **Check exit code** before reading output files
- **Parse TSV output**, not human-readable log messages

---

## ⚡ Parameters

### Search Command Parameters

#### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `--contigs FILE` | path | Input contig FASTA (plain or `.gz`) |
| `--bait FILE` | path | Input bait/reference FASTA (plain or `.gz`) |
| `--out FILE` | path | Output hit TSV path (with `--gzip`, `.gz` is appended if needed) |

#### Filtering Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--identity FLOAT` | 0 to 1 | `0.97` | Minimum alignment identity |
| `--coverage FLOAT` | 0 to 1 | `0.80` | Minimum bait coverage |
| `--min-aln-length INT` | ≥ 0 | `0` | Minimum alignment block length |
| `--terminal-tolerance INT` | ≥ 0 | `5` | Allowed unaligned bp at sequence ends |
| `--no-terminal-filter` | flag | off | Disable terminal-placement filtering |
| `--best-only` | flag | off | Keep only best contig per bait |

#### Annotation Arguments

| Option | Type | Description |
|--------|------|-------------|
| `--lineage FILE` | path | Optional bait lineage TSV (`bait_id`, `lineage` columns) |
| `--circular-list FILE` | path | Optional list of circular contig IDs (one per line) |

#### Mapping Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--preset STR` | string | `asm10` | minimap2 preset (`-x` flag) |
| `--threads INT` | ≥ 1 | `8` | Total minimap2 thread budget |
| `--minimap2-jobs INT` | ≥ 1 | `1` | Parallel minimap2 processes for splitting bait FASTA; each job receives part of the `--threads` budget |
| `--minimap2 PATH` | path | `minimap2` | minimap2 executable path |
| `--keep-paf` | flag | off | Keep intermediate PAF output |
| `--tmp-dir DIR` | path | output dir | Temporary PAF directory |

#### Contig Index Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--contig-index FILE` | path | `<contigs>.bait2contig.fai` | Text FASTA metadata/index file |
| `--rebuild-contig-index` | flag | off | Rebuild the index even if it matches the contig FASTA |
| `--index-threads INT` | ≥ 0 | `0` | Threads for building plain FASTA indexes; `0` follows `--threads` |
| `--no-contig-index` | flag | off | Disable indexing and read contigs directly when needed |

The index stores contig IDs, lengths, random-access offsets, circularity, and original FASTA headers in a plain-text FAI-like file. It does not copy full contig sequences, so it is much smaller than the FASTA. The first build still scans the contig FASTA once; later runs reuse the index when file size and mtime match. Plain FASTA indexing uses a mmap-based chunk scanner and can use multiple threads; by default it follows `--threads`, and `--index-threads` can override it. Gzip FASTA indexing is sequential because gzip streams are not efficiently splittable.

#### Contig Extraction Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--extract-contigs FILE` | path | not set | Write matched contigs to FASTA |
| `--extract-mode STR` | choice | `all` | Extraction mode: `all`, `best`, `circular`, `non-circular` |
| `--extract-min-identity FLOAT` | 0 to 1 | `--identity` | Min identity for extraction |
| `--extract-min-coverage FLOAT` | 0 to 1 | `--coverage` | Min bait coverage for extraction |
| `--extract-min-aln-length INT` | ≥ 0 | `--min-aln-length` | Min alignment length for extraction |
| `--extract-rename` | flag | off | Rename headers with hit metrics |
| `--extract-include-lineage` | flag | off | Include lineage in headers (requires `--extract-rename`) |
| `--extract-dedup` | flag | on | Deduplicate by contig ID |
| `--no-extract-dedup` | flag | off | Allow repeated sequences |

#### Resume & Output Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--resume` | flag | off | Skip if log has matching DONE marker and output is valid |
| `--rerun` | flag | off | Force recomputation (cannot use with `--resume`) |
| `--force` | flag | off | Allow overwriting without checks |
| `--gzip` | flag | off | Compress TSV, PAF, and FASTA outputs |
| `--log FILE` | path | `<out>.log` | Plain-text log path |

#### Runtime & Logging Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--monitor-interval INT` | ≥ 1 | `30` | Seconds between CPU and memory samples |
| `--no-color` | flag | off | Disable colored terminal output |
| `--quiet` | flag | off | Show only warnings and errors |
| `--verbose` | flag | off | Show detailed logs, including periodic resource snapshots |

---

### Summarize Command Parameters

#### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `--hits FILE` | path | Hit TSV from `bait2contig search` (plain or `.gz`) |
| `--out FILE` | path | Output summary TSV (with `--gzip`, `.gz` appended if needed) |

#### Filtering Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--min-identity FLOAT` | 0 to 1 | not set | Additional identity filter |
| `--min-coverage FLOAT` | 0 to 1 | not set | Additional bait coverage filter |
| `--min-aln-length INT` | ≥ 0 | not set | Additional alignment length filter |

#### Summary Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--best-hit` | flag | off | Add best-contig columns |
| `--include-contigs` | flag | off | Add ordered contig list |
| `--contig-sep STR` | string | `,` | Separator for contig list |

#### Resume & Output Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--resume` | flag | off | Skip if log has matching DONE marker |
| `--rerun` | flag | off | Force recomputation |
| `--force` | flag | off | Allow overwriting without checks |
| `--gzip` | flag | off | Compress output TSV |
| `--log FILE` | path | `<out>.log` | Plain-text log path |

#### Runtime & Logging Arguments

Same as `search`: `--monitor-interval`, `--no-color`, `--quiet`, `--verbose`.

---

## 📤 Output Formats & Path Rules

### Output Path Rules

`bait2contig` distinguishes between requested and actual output paths:

| Request | Actual Output |
|---------|---------------|
| `--out bait2contig.hits.tsv` | `bait2contig.hits.tsv` |
| `--out bait2contig.hits.tsv --gzip` | `bait2contig.hits.tsv.gz` |
| `--out bait2contig.hits.tsv.gz --gzip` | `bait2contig.hits.tsv.gz` |
| `--extract-contigs matched.fa --gzip` | `matched.fa.gz` |

**Default log path:** `<actual_out>.log`

Examples:
```text
bait2contig.hits.tsv.log
bait2contig.hits.tsv.gz.log
bait2contig.summary.tsv.log
```

> **Note:** Logs are always plain text, even when `--gzip` is used.
bait2contig.hits.tsv.log
bait2contig.hits.tsv.gz.log
bait2contig.summary.tsv.log
bait2contig.summary.tsv.gz.log
```

### Input FASTA Format

Plain and gzip-compressed FASTA are supported:

```fasta
>contig_001 len=15320 circular=true
ATGC...
>contig_002 len=8401
ATGC...
```

> **Note:** Sequence ID is the first token before whitespace. Lengths are calculated from sequence data.

### Lineage TSV Format

Optional lineage file (plain or gzip TSV with two columns):

```tsv
bait_id	lineage
Pace_16S_001	d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__;g__;s__
Pace_16S_002	d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__;g__;s__
```

> **Note:** Headerless input is accepted. Duplicate `bait_id` values are rejected.

### Circular List Format

Optional circular list (plain or gzip text, one contig ID per line):

```text
contig_001
contig_008
contig_109
```

> **Note:** Circularity is inferred from FASTA headers or `--circular-list` only. Keywords are case-insensitive: `circular=true`, `is_circular=true`, `circular`, `circ=true`.

### Search Output TSV

**Without lineage:**

```tsv
ctg_id	bait_id	identity	aln_length	cov_bait	ctg_len	is_circular
```

**With lineage:**

```tsv
ctg_id	bait_id	identity	aln_length	cov_bait	ctg_len	is_circular	lineage
contig_00012	Pace_16S_001	0.991234	1450	0.982394	32781	False	d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__;g__;s__
```

#### Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| `ctg_id` | string | Contig ID (first FASTA header token) |
| `bait_id` | string | Bait/reference ID (first FASTA header token) |
| `identity` | float | `residue_matches / alignment_block_length` (6 decimals) |
| `aln_length` | integer | PAF alignment block length |
| `cov_bait` | float | Aligned bait length ÷ bait length (6 decimals) |
| `ctg_len` | integer | Contig sequence length from FASTA |
| `is_circular` | boolean | `True` or `False` (inferred from `--circular-list` or FASTA headers) |
| `lineage` | string | Optional lineage (only when `--lineage` provided) |

> **Note:** If no alignments pass filters, TSV is written with only the header.

### Summary Output TSV

**Default columns:**

```tsv
bait_id	contig_count	circular_contig_count	total_ctg_len	mean_ctg_len	max_ctg_len	mean_identity	max_identity	mean_cov_bait	max_cov_bait	mean_aln_length	max_aln_length
```

If lineage is present in hits, it appears after `bait_id`. `--best-hit` adds best-hit columns; `--include-contigs` adds contig list.

#### Default Summary Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| `bait_id` | string | Bait/reference ID |
| `lineage` | string | Optional (only when input hits contain lineage) |
| `contig_count` | integer | Unique contigs anchored by bait |
| `circular_contig_count` | integer | Circular unique contigs |
| `total_ctg_len` | integer | Sum of unique contig lengths |
| `mean_ctg_len` | float | Mean unique contig length (2 decimals) |
| `max_ctg_len` | integer | Maximum unique contig length |
| `mean_identity` | float | Mean identity (6 decimals) |
| `max_identity` | float | Max identity (6 decimals) |
| `mean_cov_bait` | float | Mean bait coverage (6 decimals) |
| `max_cov_bait` | float | Max bait coverage (6 decimals) |
| `mean_aln_length` | float | Mean alignment length (2 decimals) |
| `max_aln_length` | integer | Max alignment length |

#### `--best-hit` Additional Columns

| Column | Type | Description |
|--------|------|-------------|
| `best_ctg_id` | string | Best contig (by ranking rules) |
| `best_identity` | float | Best-hit identity (6 decimals) |
| `best_cov_bait` | float | Best-hit bait coverage (6 decimals) |
| `best_aln_length` | integer | Best-hit alignment length |
| `best_ctg_len` | integer | Best-hit contig length |
| `best_is_circular` | boolean | Circularity of best hit |

#### `--include-contigs` Additional Column

| Column | Type | Description |
|--------|------|-------------|
| `contigs` | string | Unique contig IDs ordered by ranking, joined with `--contig-sep` |

---

## 🔬 Advanced Features

### Identity and Coverage Calculations

`bait2contig` parses minimap2 PAF and calculates:

```
identity = residue_matches / alignment_block_length
cov_bait = aligned_bait_length / bait_length
```

More explicitly:

```
cov_bait = (query_end - query_start) / query_len
```

**Filtering logic:**
- `identity >= --identity`
- `cov_bait >= --coverage`
- `aln_length >= --min-aln-length`

For partial alignments (`cov_bait < 1.0`), alignment must touch at least one bait end and one contig end within `--terminal-tolerance` bases (default 5 bp). Use `--no-terminal-filter` to disable.

**Best-hit ranking order:**
1. Identity (descending)
2. Bait coverage (descending)
3. Alignment length (descending)
4. Contig length (descending)

### Examples

#### Basic search:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --identity 0.97 \
  --coverage 0.80 \
  --threads 32
```

#### Search with lineage:

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

#### Search with gzip output:

```bash
bait2contig search \
  --contigs contigs.fa.gz \
  --bait bait.fa.gz \
  --out bait2contig.hits.tsv \
  --gzip \
  --threads 32
```

Actual output: `bait2contig.hits.tsv.gz`

#### Summarize hits:

```bash
bait2contig summarize \
  --hits bait2contig.hits.tsv \
  --out bait2contig.summary.tsv \
  --best-hit \
  --include-contigs
```

#### Summarize with gzip:

```bash
bait2contig summarize \
  --hits bait2contig.hits.tsv.gz \
  --out bait2contig.summary.tsv \
  --gzip \
  --best-hit \
  --include-contigs
```

Actual output: `bait2contig.summary.tsv.gz`

### Contig Extraction

Extract all matched contigs:

```bash
bait2contig search \
  --contigs contigs.fa \
  --bait bait.fa \
  --out bait2contig.hits.tsv \
  --extract-contigs matched_contigs.fa
```

Extract best contig per bait:

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

Rename extracted headers with metrics:

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

Renamed header format:

```
>{ctg_id}|bait={bait_id}|identity={identity:.6f}|cov_bait={cov_bait:.6f}|aln_length={aln_length}|ctg_len={ctg_len}|circular={is_circular}
```

### Gzip Output Examples

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
```
bait2contig.hits.tsv.gz
matched_contigs.fa.gz
bait2contig.hits.tsv.gz.log
```

### Resume and Rerun

Resume relies on plain-text logs with `[BAIT2CONTIG_DONE]` markers. It verifies output path, size, command, and parameters before skipping:

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

### Large Contig FASTA Indexing

For large contig FASTA inputs, `search` no longer loads all contig sequences before running `minimap2`. It loads bait sequences, runs `minimap2`, then looks up only hit contigs through the text FASTA index.

Default index path:

```text
<contigs>.bait2contig.fai
```

Use `--contig-index` to place the index on fast local storage, or `--rebuild-contig-index` after manually replacing a FASTA without changing its mtime. With the default `--index-threads 0`, plain FASTA indexing follows `--threads`; set `--index-threads` explicitly if you want a different index-building thread count. Plain FASTA inputs support fast offset-based sequence extraction. Gzip FASTA inputs can still be used, but sequence extraction may fall back to streaming because gzip is not efficiently seekable.

### CPU and Memory Monitoring

Resource monitoring always samples CPU and memory for final summary statistics. Periodic resource snapshot lines are written only with `--verbose`:

```
2026-06-12 09:30:00 [RESOURCE] stage=running_minimap2 elapsed=30.0s rss_mb=842.5 cpu_percent=1250.4
```

During large FASTA loading, interactive terminals show live progress bars with bytes, sequence counts, and parsed bases. CPU percent may exceed 100 because `minimap2` is multithreaded.

If `running_minimap2` still uses only one or two cores after `--threads 32`, the bait FASTA may be too small for one minimap2 process to keep all workers busy. Use `--minimap2-jobs` to split bait records across multiple minimap2 processes, for example `--threads 32 --minimap2-jobs 8`. This increases parallelism but can increase memory and I/O because each job scans the contig input.

### Help Colors

Help output uses ANSI colors when stdout is a terminal:

```bash
bait2contig --help
bait2contig search --help
bait2contig summarize --help
```

Disable colors:
```bash
bait2contig --help --no-color
NO_COLOR=1 bait2contig --help
```

Force colors:
```bash
CLICOLOR_FORCE=1 bait2contig --help
```

### Error Handling

`bait2contig` validates input and exits with concise `ERROR:` messages. It suggests the nearest valid command or option when possible, including a minimal example command.

---

## 📝 Minimal Example

### Input Files

**contigs.fa:**
```fasta
>contig_001 len=12 circular=true
AAACCCGGGTTT
>contig_002 len=12
TTTGGGCCCAAA
```

**bait.fa:**
```fasta
>bait_001
AAACCCGGGTTT
```

**lineage.tsv:**
```tsv
bait_id	lineage
bait_001	d__Example;p__Example;c__;o__;f__;g__;s__
```

### Expected Outputs

**Hit TSV shape:**
```tsv
ctg_id	bait_id	identity	aln_length	cov_bait	ctg_len	is_circular	lineage
contig_001	bait_001	1.000000	12	1.000000	12	True	d__Example;p__Example;c__;o__;f__;g__;s__
```

**Summary TSV shape** (with `--best-hit --include-contigs`):
```tsv
bait_id	lineage	contig_count	circular_contig_count	total_ctg_len	mean_ctg_len	max_ctg_len	mean_identity	max_identity	mean_cov_bait	max_cov_bait	mean_aln_length	max_aln_length	best_ctg_id	best_identity	best_cov_bait	best_aln_length	best_ctg_len	best_is_circular	contigs
bait_001	d__Example;p__Example;c__;o__;f__;g__;s__	1	1	12	12.00	12	1.000000	1.000000	1.000000	1.000000	12.00	12	contig_001	1.000000	1.000000	12	12	True	contig_001
```

---

## 🔧 Troubleshooting

### minimap2 Not Found

**Error:**
```
ERROR: minimap2 was not found. Please install minimap2 or provide its path with --minimap2.
```

**Solution:**
1. Install minimap2: `conda install -c bioconda minimap2` or `brew install minimap2`
2. Or specify the path explicitly:
   ```bash
   bait2contig search --minimap2 /usr/local/bin/minimap2 --contigs contigs.fa --bait bait.fa --out hits.tsv
   ```

### No Hits Found

**Symptom:** Output TSV contains only header, no data rows.

**Causes & Solutions:**
- Identity threshold too strict → Lower `--identity` value
- Coverage threshold too strict → Lower `--coverage` value
- Bait and contigs too divergent → Check sequence similarity manually
- Bait too short → PAF may not report short alignments; increase alignment length minimum if needed

**Debug:**
```bash
# Keep PAF to inspect raw minimap2 output
bait2contig search ... --keep-paf
```

### Slow Performance

**To improve speed:**
- Increase `--threads` (e.g., `--threads 32`)
- If minimap2 underuses CPUs with many bait records, add `--minimap2-jobs 4` or `--minimap2-jobs 8`
- For first-time plain FASTA indexing, `--index-threads 0` follows `--threads`; gzip FASTA indexing remains sequential
- Use `--best-only` to stop early per bait
- Pre-filter contigs if possible (subset input FASTA)
- Check system resources: `top` or Task Manager

**Expected times:**
- Mapping 100 baits to 1000 contigs: ~1-5 seconds (varies by contig size)
- Large metagenomics assemblies: minutes to hours depending on scale

### Out of Memory

**Solutions:**
- Reduce `--threads` (uses less memory)
- Process contigs in batches if possible
- Check available RAM: `free -h` (Linux) or Task Manager (Windows)

### Ambiguous Lineage or Missing Data

**Issue:** Lineage column empty or unexpected format.

**Check:**
- Lineage TSV must have tab-delimited columns `bait_id` and `lineage`
- Duplicate `bait_id` values are rejected
- File must be plain text or gzip (`.gz`), not other formats

### Output File Not Created

**Causes:**
- Output directory doesn't exist → `mkdir -p output_dir` first
- No write permissions → Check directory permissions
- Disk full → Free up space

---

## ❓ FAQ

### Q: What's the difference between `--identity` and `--extract-min-identity`?

**A:** `--identity` filters the initial search results. `--extract-min-identity` applies an additional filter only to contigs being extracted to FASTA. This lets you extract with looser criteria while keeping stricter hits in the TSV.

### Q: Why are some contigs marked `is_circular=False` even though they wrap around?

**A:** `bait2contig` only infers circularity from FASTA headers or `--circular-list`. It does not perform de novo detection. Add header keywords like `circular=true` or use `--circular-list`.

### Q: How do I use the lineage file?

**A:** Pass a TSV with `bait_id` and `lineage` columns:
```bash
bait2contig search --lineage lineage.tsv ...
```
Lineage is appended to hit TSV and summary TSV for easy taxonomic reference.

### Q: Can I use this with long reads (nanopore, pacbio)?

**A:** Yes, if your reference contigs are long-read assemblies. Change `--preset` to match: `--preset map-pb` (PacBio) or `--preset map-ont` (Nanopore).

### Q: What is PAF and why keep it?

**A:** PAF is minimap2's output format (Pairwise Alignment Format). Keep it with `--keep-paf` to inspect raw alignments or use for downstream analysis.

### Q: How do I resume a partially completed run?

**A:** Use `--resume`. It checks for a matching `[BAIT2CONTIG_DONE]` marker in the log and verifies the output matches the command:
```bash
bait2contig search ... --resume
```
If parameters changed, it reruns automatically.

### Q: Can I run this in parallel for multiple datasets?

**A:** Yes. Launch separate processes with different `--contigs`, `--bait`, or `--out` arguments. Use `--resume` to avoid redundant work if a run crashes.

### Q: What version of minimap2 is required?

**A:** Any recent version (0.13+). Older versions may have compatibility issues; update if you encounter problems.

### Q: How do I extract only the best contig per bait?

**A:** Use search with `--extract-mode best`:
```bash
bait2contig search ... --extract-contigs best.fa --extract-mode best
```

### Q: Can I filter the summary output further?

**A:** Yes, use `summarize` filtering options:
```bash
bait2contig summarize --hits hits.tsv --out summary.tsv \
  --min-identity 0.98 --min-coverage 0.90
```

---

## ⚠️ Notes and Limitations

- **No de novo circularity detection** – inferred from FASTA headers or `--circular-list` only
- **No SAM/BAM parsing** – uses minimap2 PAF output
- **No taxonomy rank splitting** – lineage passed through as-is
- **No database building** – operates on input files directly
- **Single aligner backend** – minimap2 only
- **No web interface** – command-line only

---

## 📄 License

Released under the **MIT License**. See [LICENSE](LICENSE) for full details.
