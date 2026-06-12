# bait2contig Development Notes

This document records implementation details that are easy to forget when maintaining `bait2contig`. It is intended for developers. User-facing usage is documented in `README.md`.

## Project Constraints

- Project name, package name, CLI command, logs, tests, README examples, and help text use `bait2contig`.
- User-facing text in source comments, docstrings, README, help messages, logs, and errors is English.
- The package intentionally has a shallow source tree:
  - `bait2contig/cli.py`
  - `bait2contig/core.py`
  - `bait2contig/fasta_index.py`
  - `bait2contig/io.py`
  - `bait2contig/log.py`
  - `bait2contig/summary.py`
- Required Python dependencies are limited to terminal/help/logging support: `rich`, `rich-argparse`, and `psutil`.
- Resource monitoring still has a standard-library fallback when `psutil` is unavailable.
- The CLI uses `argparse`, not `click` or `typer`.
- The mapping backend is `minimap2`; no alternate aligners are currently supported.

## Module Responsibilities

`bait2contig/cli.py`

- Builds the top-level parser and the `search` and `summarize` subcommands.
- Defines grouped help output.
- Adds colorized help output when stdout is a terminal.
- Implements `SmartArgumentParser` for actionable command-line errors.
- Dispatches to `run_search()` or `run_summarize()`.

`bait2contig/core.py`

- Implements the search workflow.
- Validates search arguments before running minimap2.
- Runs minimap2 in the required target/query order.
- Parses PAF through `io.parse_paf()`.
- Annotates hits with circularity and lineage.
- Applies hit filters and optional best-hit selection.
- Writes hit-level TSV output.
- Handles optional kept PAF output.
- Handles optional contig extraction.
- Writes search START, DONE, and FAILED log marker blocks.

`bait2contig/fasta_index.py`

- Builds and validates text FASTA indexes for large contig FASTA files.
- Stores contig IDs, headers, lengths, circularity, and plain-FASTA sequence offsets.
- Fetches only requested contig sequences for extraction.
- Falls back to stream-based subset reading when random access is unavailable.

`bait2contig/io.py`

- Implements `open_text()` for plain and gzip text files.
- Parses FASTA files and keeps the original header for extraction.
- Infers circularity from FASTA headers.
- Reads lineage TSV and circular contig list files.
- Parses minimap2 PAF lines into `PafHit`.
- Reads and writes TSV files.

`bait2contig/summary.py`

- Reads hit-level TSV output from `search`.
- Applies optional additional filters.
- Deduplicates repeated `bait_id` plus `ctg_id` rows by best-hit ranking.
- Computes bait-level summary statistics.
- Writes summary TSV output.
- Writes summarize START, DONE, and FAILED log marker blocks.

`bait2contig/log.py`

- Writes grep-friendly terminal and plain-text log output.
- Uses ANSI colors only for terminal output, never for log files.
- Implements standardized resume marker parsing.
- Implements lightweight runtime resource monitoring.
- Uses `psutil` if available and standard-library fallback otherwise.

## Search Data Flow

1. `cli.py` parses `bait2contig search` arguments.
2. `core.validate_search_args()` checks:
   - conflicting options such as `--resume` plus `--rerun`
   - option dependencies such as `--extract-include-lineage` requiring `--extract-rename`
   - numeric ranges
   - required input file existence
   - output path sanity
3. Actual output paths are resolved with `gzip_output_path()`.
4. Resume checks run before recomputation when `--resume` is used.
5. Bait FASTA, lineage, and circular-list inputs are loaded.
6. The contig FASTA is not loaded into memory by default.
7. minimap2 is executed as:

   ```bash
   minimap2 -x {preset} -t {threads} {contigs} {bait} > {tmp_paf}
   ```

   This order is required. Contigs are the target/reference and bait sequences are the query. With `--minimap2-jobs > 1`, bait records are split into chunks, minimap2 is run concurrently for each chunk, and PAF chunks are concatenated in chunk order. In PAF:

   - query ID is `bait_id`
   - target ID is `ctg_id`

8. PAF hits are parsed into `PafHit`.
9. If contig metadata or sequences are needed, a text contig FASTA index is reused or built.
10. PAF hits are converted into `SearchHit`.
11. Hits are filtered by:
   - identity
   - bait coverage
   - alignment length
   - partial-hit terminal placement, unless disabled
12. If `--best-only` is set, only one best contig per bait is retained.
13. The hit TSV is written even when no hit passes filters.
14. Optional contig extraction runs from the annotated hit list.
15. DONE or FAILED marker blocks are written to the plain-text log.

## PAF-Derived Metrics

PAF fields used:

| PAF index | Meaning |
| --- | --- |
| 0 | query ID, stored as `bait_id` |
| 1 | query length, stored as `bait_len` |
| 2 | query start, stored as `bait_start` |
| 3 | query end, stored as `bait_end` |
| 5 | target ID, stored as `ctg_id` |
| 6 | target length, stored as `ctg_len` unless FASTA length is available |
| 7 | target start, stored as `ctg_start` |
| 8 | target end, stored as `ctg_end` |
| 9 | residue matches |
| 10 | alignment block length, stored as `aln_length` |

Derived metrics:

```text
identity = residue_matches / alignment_block_length
cov_bait = (query_end - query_start) / query_len
```

Output formatting:

- `identity` uses six decimal places.
- `cov_bait` uses six decimal places.

## Best-Hit Ranking

The same ranking is used by search best-hit selection, extraction best mode, and summary deduplication:

1. identity descending
2. cov_bait descending
3. aln_length descending
4. ctg_len descending

The implementation uses tuple comparison:

```python
(identity, cov_bait, aln_length, ctg_len)
```

## Partial-Hit Terminal Filter

By default, partial bait hits must be terminal. This was added to reduce likely chimeric or problematic contig hits where only an internal bait fragment aligns to an internal contig region.

Implementation: `core.is_terminal_partial_hit()`.

Rules:

- If `cov_bait >= 1.0`, the hit passes this filter.
- If `bait_len <= 0`, the hit passes this filter because terminal placement cannot be evaluated.
- Otherwise, the hit must touch:
  - at least one bait end within `--terminal-tolerance`
  - and at least one contig end within `--terminal-tolerance`
- Default tolerance is 5 bp.
- `--no-terminal-filter` disables this check.

Important: The hit TSV does not currently include bait or contig alignment coordinates. If users need to audit this filter directly, use `--keep-paf` and inspect the PAF file.

## Circularity Logic

Circularity is not detected de novo.

Search uses this priority:

1. If `--circular-list` is provided, only IDs in that list are `True`.
2. Otherwise, circularity is inferred from the contig FASTA header.
3. Otherwise, `False`.

Supported header tokens are case-insensitive:

- `circular=true`
- `is_circular=true`
- `circular`
- `circ=true`

## Lineage Logic

Lineage TSV accepts header and headerless two-column input:

```text
bait_id	lineage
```

Duplicate `bait_id` values raise an error.

Lineage IDs not present in bait FASTA are warnings, not fatal errors.

When a lineage file is provided, search output adds a `lineage` column. Summary preserves lineage if that column exists in the input hits TSV.

## Contig Extraction

Extraction is only performed when `--extract-contigs` is provided.

Extraction thresholds default to the hit-table thresholds:

- `--extract-min-identity` defaults to `--identity`
- `--extract-min-coverage` defaults to `--coverage`
- `--extract-min-aln-length` defaults to `--min-aln-length`

Extraction modes:

- `all`: all annotated hits passing extraction filters
- `best`: best contig per bait
- `circular`: only circular hits
- `non-circular`: only non-circular hits

Deduplication is enabled by default. With deduplication, the same `ctg_id` is written once. If `--no-extract-dedup` is used without `--extract-rename`, warn because duplicate FASTA IDs may be written.

Renamed FASTA header format:

```text
{ctg_id}|bait={bait_id}|identity={identity:.6f}|cov_bait={cov_bait:.6f}|aln_length={aln_length}|ctg_len={ctg_len}|circular={is_circular}
```

If `--extract-include-lineage` is used, append:

```text
|lineage={lineage}
```

`--extract-include-lineage` requires `--extract-rename`.

## Contig FASTA Index

The index is enabled by default for `search`.

Default index path:

```text
<contigs>.bait2contig.fai
```

Index validity is based on:

- absolute contig FASTA path
- file size
- file mtime in nanoseconds
- index schema version

For plain FASTA, the index stores sequence offsets, so extraction can seek directly to hit contigs. Plain FASTA index building uses mmap chunk scanning and can use multiple threads; with `--index-threads 0`, the index thread budget follows `--threads`, and explicit `--index-threads` overrides it. For gzip FASTA, the index stores metadata but index building remains sequential and sequence extraction may need stream-based subset reading because gzip is not efficiently seekable.

The index intentionally does not store full sequences. It is a plain-text TSV-like file, not SQLite.

## Gzip Rules

`open_text(path, mode)` is the single helper for plain and gzip text IO.

Output path behavior:

- If `--gzip` is set and `--out` does not end with `.gz`, append `.gz`.
- If `--gzip` is set and `--out` already ends with `.gz`, do not append another suffix.
- If `--extract-contigs` is used with `--gzip`, the extraction output path follows the same rule.
- Logs are never gzip-compressed.
- Kept PAF files are gzip-compressed only when both `--keep-paf` and `--gzip` are set.

## Resume Mechanism

There are no JSON checkpoint files. The plain-text log is the checkpoint.

Standard markers:

- `[BAIT2CONTIG_START]`
- `[BAIT2CONTIG_DONE]`
- `[BAIT2CONTIG_FAILED]`

Resume checks:

1. Log file exists.
2. Output exists and is non-empty.
3. Latest relevant marker block is DONE, not START or FAILED.
4. `command=...`, `status=success`, and `exit_code=0`.
5. DONE output path matches the current actual output path.
6. START parameters match key current parameters.

Search resume parameter keys include:

- output path
- gzip flag
- input paths
- thresholds
- best-only flag
- terminal filter settings
- preset
- minimap2 version
- extraction settings

Summarize resume parameter keys include:

- output path
- gzip flag
- hits path
- additional filters
- best-hit flag
- include-contigs flag
- contig separator

If resume checks fail, the command reruns.

## Logging and Resource Monitoring

Log files remain plain text. Terminal output may be colored, but log files must never contain ANSI escape codes.

Resource monitoring behavior:

- Prefer `psutil` for detailed process and child-process CPU/RSS monitoring.
- Fall back to standard-library `resource` plus Linux `/proc` child CPU sampling when available.
- Attempt to include child process usage for minimap2.
- Sample resource usage periodically for final DONE or FAILED summary metrics.
- Write periodic `[RESOURCE]` lines only when `--verbose` is enabled.
- Resource monitoring failures must not fail the main command.
- CPU percent can exceed 100 because minimap2 is multithreaded.

## CLI Error Handling

`SmartArgumentParser` improves common command-line errors:

- missing command
- misspelled command, with suggestion
- misspelled option, with suggestion
- missing required options
- missing option value
- invalid type values
- invalid choice values

Long-option abbreviation is disabled with `allow_abbrev=False`. This is important because `argparse` would otherwise accept partial long options such as `--identit` as `--identity`, which hides user input mistakes.

Workflow-level validation then checks:

- input file existence
- output paths that point to directories
- conflicting flags
- option dependency rules
- numeric ranges

## Summary Data Flow

1. `cli.py` parses `bait2contig summarize`.
2. `summary.validate_summary_args()` checks input path and option values.
3. Resume checks run if requested.
4. `read_hit_rows()` loads the hit TSV.
5. Optional additional filters are applied.
6. Repeated `bait_id` plus `ctg_id` pairs are deduplicated by best-hit ranking.
7. Per-bait statistics are computed.
8. Summary TSV is written even when no rows pass filters.
9. DONE or FAILED marker blocks are written.

## Summary Statistics

For each `bait_id`, summary operates on unique contigs after deduplication.

Statistics:

- `contig_count`: number of unique contigs
- `circular_contig_count`: number of unique contigs marked circular
- `total_ctg_len`: sum of unique contig lengths
- `mean_ctg_len`: mean unique contig length, two decimals
- `max_ctg_len`: maximum unique contig length
- `mean_identity`: six decimals
- `max_identity`: six decimals
- `mean_cov_bait`: six decimals
- `max_cov_bait`: six decimals
- `mean_aln_length`: two decimals
- `max_aln_length`: maximum alignment length

When `--include-contigs` is used, contigs are ordered by best-hit ranking.

## External Invocation Contract

External callers should treat `bait2contig` as a command-line program, not as a stable Python API. Internal Python functions may change.

Recommended subprocess behavior:

- Pass arguments as a list, not one shell string.
- Capture stderr for logs and errors.
- Use exit code 0 as success.
- Treat exit code 1 as runtime or validation failure.
- Treat exit code 2 as command-line parsing failure.
- Use TSV outputs as machine-readable data.
- Do not parse colored help output.
- Do not parse ordinary log text except standardized marker blocks.

## Testing

Run:

```bash
python -m pytest
```

Current tests cover:

- plain and gzip IO
- FASTA parsing
- lineage parsing
- circular header parsing
- PAF parsing
- identity and coverage math
- best-hit ranking
- hit filtering
- partial-hit terminal filtering
- contig extraction modes
- renamed extraction headers
- summary statistics
- gzip hit reading
- log parser and resume checks
- CLI help and error handling

When changing CLI behavior, update tests in `tests/test_core.py`. When changing file parsing, update `tests/test_io.py`. When changing summary statistics, update `tests/test_summary.py`.

## Manual Verification

After code changes, run:

```bash
python -m pytest
bait2contig --help
bait2contig search --help
bait2contig summarize --help
```

If the package is not installed, use:

```bash
PYTHONPATH="$PWD" python -m pytest
PYTHONPATH="$PWD" python -m bait2contig.cli --help
PYTHONPATH="$PWD" python -m bait2contig.cli search --help
PYTHONPATH="$PWD" python -m bait2contig.cli summarize --help
```

To verify help colors in a non-interactive environment:

```bash
CLICOLOR_FORCE=1 bait2contig search --help
```

## Extension Guidelines

Keep changes conservative:

- Prefer extending existing modules over adding new files.
- Keep user-facing text in English.
- Preserve TSV column order unless making an intentional format change.
- Preserve log marker names and key names for resume compatibility.
- Preserve minimap2 target/query order.
- Add tests for behavior changes.
- Update README and this file when changing CLI flags, output schema, resume parameters, or filtering logic.

## Non-Goals

The current MVP does not implement:

- JSON checkpoint files
- de novo circularity detection
- multiple aligners
- SAM/BAM parsing
- taxonomy rank splitting
- complex database indexing
- a web interface
- heavy runtime dependencies
