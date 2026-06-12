"""Command-line interface for bait2contig."""

from __future__ import annotations

import argparse
import difflib
import os
import re
import sys

from rich.console import Console
from rich_argparse import RichHelpFormatter

from .core import LoggedSearchError, SearchError, run_search
from .summary import LoggedSummaryError, SummaryError, run_summarize


class HelpFormatter(RichHelpFormatter):
    """Readable argparse formatter for grouped bait2contig help."""

    use_color = False
    group_name_formatter = str

    def __init__(self, prog: str) -> None:
        console = Console(stderr=True, color_system="auto" if self.use_color else None)
        super().__init__(prog, max_help_position=32, width=100, console=console)

    def _get_help_string(self, action) -> str:
        help_text = action.help or ""
        if (
            not help_text
            or action.required
            or action.default is argparse.SUPPRESS
            or "%(default)" in help_text
            or "default:" in help_text.lower()
        ):
            return help_text
        if isinstance(action, argparse._StoreFalseAction):
            return help_text
        if isinstance(action, argparse._StoreTrueAction):
            return f"{help_text} (default: {'on' if action.default else 'off'})"
        default = "not set" if action.default is None else "%(default)s"
        return f"{help_text} (default: {default})"


class SmartArgumentParser(argparse.ArgumentParser):
    """ArgumentParser with concise, actionable error messages."""

    valid_commands = ("search", "summarize")

    def __init__(self, *args, example: str | None = None, **kwargs) -> None:
        kwargs.setdefault("allow_abbrev", False)
        super().__init__(*args, **kwargs)
        self.example = example

    def error(self, message: str) -> None:
        self.print_usage(sys.stderr)
        console = Console(stderr=True, color_system="auto" if HelpFormatter.use_color else None)
        console.print(f"[bold red]ERROR:[/bold red] {self._friendly_error(message)}")
        hint = self._error_hint(message)
        if hint:
            console.print(hint, style="yellow")
        if self.example:
            console.print(f"Example:\n  {self.example}", style="dim", soft_wrap=True)
        console.print(f"Run '{self.prog} --help' for details.")
        self.exit(2)

    def _friendly_error(self, message: str) -> str:
        if message.startswith("the following arguments are required:"):
            missing = message.split(":", 1)[1].strip()
            return f"missing required argument(s): {missing}"
        if "expected one argument" in message:
            match = re.search(r"argument ([^:]+): expected one argument", message)
            if match:
                return f"{match.group(1)} requires a value"
        invalid_choice = re.search(r"argument ([^:]+): invalid choice: '([^']+)'", message)
        if invalid_choice:
            argument, value = invalid_choice.groups()
            choices = self._choices_for(argument)
            if choices:
                return f"invalid value for {argument}: {value}. Choose one of: {', '.join(choices)}"
        invalid_typed_value = re.search(r"argument ([^:]+): invalid ([^ ]+) value: '([^']+)'", message)
        if invalid_typed_value:
            argument, expected_type, value = invalid_typed_value.groups()
            return f"invalid value for {argument}: {value}. Expected {expected_type}."
        if message.startswith("unrecognized arguments:"):
            unknown = message.split(":", 1)[1].strip()
            return f"unrecognized argument(s): {unknown}"
        return message

    def _error_hint(self, message: str) -> str:
        invalid_command = re.search(r"argument <command>: invalid choice: '([^']+)'", message)
        if invalid_command:
            suggestion = closest_match(invalid_command.group(1), self.valid_commands)
            if suggestion:
                return f"Did you mean '{suggestion}'?"

        if message.startswith("unrecognized arguments:"):
            unknown = message.split(":", 1)[1].strip().split()
            options = collect_option_strings(self)
            hints = []
            for token in unknown:
                if token.startswith("-"):
                    suggestion = closest_match(token, options)
                    if suggestion:
                        hints.append(f"Did you mean {suggestion}?")
            if hints:
                return "\n".join(hints)

        invalid_choice = re.search(r"argument ([^:]+): invalid choice: '([^']+)'", message)
        if invalid_choice:
            argument, value = invalid_choice.groups()
            suggestion = closest_match(value, self._choices_for(argument))
            if suggestion:
                return f"Did you mean '{suggestion}'?"
        return ""

    def _choices_for(self, argument: str) -> list[str]:
        if argument == "<command>":
            return list(self.valid_commands)
        for action in self._actions:
            labels = list(action.option_strings) or ([action.dest] if action.dest else [])
            if argument in labels and action.choices:
                return [str(choice) for choice in action.choices]
        return []


def collect_option_strings(parser: argparse.ArgumentParser) -> list[str]:
    """Collect option strings from one parser and its subcommands."""

    options: list[str] = []
    for action in parser._actions:
        options.extend(action.option_strings)
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                options.extend(collect_option_strings(subparser))
    return options


def closest_match(value: str, choices: list[str] | tuple[str, ...]) -> str | None:
    matches = difflib.get_close_matches(value, list(choices), n=1, cutoff=0.6)
    return matches[0] if matches else None


def should_color_help(argv: list[str] | None = None) -> bool:
    """Return whether help output should use ANSI colors."""

    argv = argv or []
    if "--no-color" in argv:
        return False
    force = os.environ.get("CLICOLOR_FORCE")
    if force and force != "0":
        return True
    if os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()


def build_parser() -> argparse.ArgumentParser:
    parser = SmartArgumentParser(
        prog="bait2contig",
        usage="bait2contig <command> [options]",
        description="Find contigs corresponding to bait/reference sequences.",
        formatter_class=HelpFormatter,
        example="bait2contig search --contigs contigs.fa --bait bait.fa --out bait2contig.hits.tsv",
    )
    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        metavar="<command>",
        parser_class=SmartArgumentParser,
    )

    search = subparsers.add_parser(
        "search",
        prog="bait2contig search",
        help="Find contigs matching bait/reference sequences.",
        description="Find contigs matching bait/reference sequences.",
        formatter_class=HelpFormatter,
        example="bait2contig search --contigs contigs.fa --bait bait.fa --out bait2contig.hits.tsv",
    )
    search.set_defaults(func=run_search)
    add_search_arguments(search)

    summarize = subparsers.add_parser(
        "summarize",
        prog="bait2contig summarize",
        help="Summarize contigs anchored by each bait.",
        description="Summarize contigs anchored by each bait.",
        formatter_class=HelpFormatter,
        example="bait2contig summarize --hits bait2contig.hits.tsv --out bait2contig.summary.tsv",
    )
    summarize.set_defaults(func=run_summarize)
    add_summarize_arguments(summarize)
    return parser


def add_search_arguments(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--contigs", metavar="FILE", required=True, help="Input contig FASTA.")
    required.add_argument("--bait", metavar="FILE", required=True, help="Input bait/reference FASTA.")
    required.add_argument("--out", metavar="FILE", required=True, help="Output hit TSV.")

    filtering = parser.add_argument_group("Filtering arguments")
    filtering.add_argument("--identity", type=float, default=0.97, metavar="FLOAT", help="Minimum alignment identity.")
    filtering.add_argument("--coverage", type=float, default=0.80, metavar="FLOAT", help="Minimum bait coverage.")
    filtering.add_argument(
        "--min-aln-length",
        type=int,
        default=0,
        metavar="INT",
        help="Minimum alignment length.",
    )
    filtering.add_argument(
        "--terminal-tolerance",
        type=int,
        default=5,
        metavar="INT",
        help="Allowed unaligned bases at sequence ends for partial bait alignments.",
    )
    filtering.add_argument(
        "--no-terminal-filter",
        dest="terminal_filter",
        action="store_false",
        help="Do not require partial bait alignments to touch bait and contig ends. Default: terminal filter enabled.",
    )
    filtering.add_argument("--best-only", action="store_true", help="Keep only the best contig for each bait.")
    parser.set_defaults(terminal_filter=True)

    annotation = parser.add_argument_group("Annotation arguments")
    annotation.add_argument(
        "--lineage",
        metavar="FILE",
        help="Optional bait lineage TSV with two columns: bait_id and lineage.",
    )
    annotation.add_argument(
        "--circular-list",
        metavar="FILE",
        help="Optional file of circular contig IDs, one ID per line.",
    )

    mapping = parser.add_argument_group("Mapping arguments")
    mapping.add_argument("--preset", default="asm10", metavar="STR", help="Minimap2 preset.")
    mapping.add_argument("--threads", type=int, default=8, metavar="INT", help="Total minimap2 thread budget.")
    mapping.add_argument(
        "--minimap2-jobs",
        type=int,
        default=1,
        metavar="INT",
        help="Parallel minimap2 processes for splitting bait FASTA.",
    )
    mapping.add_argument("--minimap2", default="minimap2", metavar="PATH", help="Path to minimap2 executable.")
    mapping.add_argument("--keep-paf", action="store_true", help="Keep intermediate PAF output.")
    mapping.add_argument("--tmp-dir", metavar="DIR", help="Temporary directory. Default: output directory.")

    index = parser.add_argument_group("Contig index arguments")
    index.add_argument(
        "--contig-index",
        metavar="FILE",
        help="Text FASTA index path. Default: <contigs>.bait2contig.fai.",
    )
    index.add_argument(
        "--rebuild-contig-index",
        action="store_true",
        help="Rebuild the contig FASTA index even when it appears current.",
    )
    index.add_argument(
        "--index-threads",
        type=int,
        default=0,
        metavar="INT",
        help="Threads for building plain FASTA indexes. Use 0 to follow --threads.",
    )
    index.add_argument(
        "--no-contig-index",
        action="store_true",
        help="Disable contig indexing and read contigs directly when needed.",
    )

    extraction = parser.add_argument_group("Contig extraction arguments")
    extraction.add_argument("--extract-contigs", metavar="FILE", help="Write matched contig sequences to this FASTA.")
    extraction.add_argument(
        "--extract-mode",
        choices=["all", "best", "circular", "non-circular"],
        default="all",
        metavar="STR",
        help="Which contigs to extract.",
    )
    extraction.add_argument(
        "--extract-min-identity",
        type=float,
        metavar="FLOAT",
        help="Minimum identity for extracted contigs. Default: --identity.",
    )
    extraction.add_argument(
        "--extract-min-coverage",
        type=float,
        metavar="FLOAT",
        help="Minimum bait coverage for extracted contigs. Default: --coverage.",
    )
    extraction.add_argument(
        "--extract-min-aln-length",
        type=int,
        metavar="INT",
        help="Minimum alignment length for extracted contigs. Default: --min-aln-length.",
    )
    extraction.add_argument(
        "--extract-rename",
        action="store_true",
        help="Rename extracted contig FASTA headers to include bait information.",
    )
    extraction.add_argument(
        "--extract-include-lineage",
        action="store_true",
        help="Include lineage in renamed FASTA headers.",
    )
    dedup = extraction.add_mutually_exclusive_group()
    dedup.add_argument(
        "--extract-dedup",
        dest="extract_dedup",
        action="store_true",
        help="Deduplicate extracted contigs. Default: enabled.",
    )
    dedup.add_argument(
        "--no-extract-dedup",
        dest="extract_dedup",
        action="store_false",
        help="Do not deduplicate extracted contigs.",
    )
    parser.set_defaults(extract_dedup=True)

    output = parser.add_argument_group("Resume and output arguments")
    output.add_argument("--resume", action="store_true", help="Resume if the log shows a previous successful run.")
    output.add_argument("--rerun", action="store_true", help="Force rerun even if output and success log exist.")
    output.add_argument("--force", action="store_true", help="Overwrite existing output.")
    output.add_argument("--gzip", action="store_true", help="Compress output TSV, kept PAF, and extracted FASTA outputs.")
    output.add_argument("--log", metavar="FILE", help="Log file. Default: <actual_out>.log.")

    runtime = parser.add_argument_group("Runtime and logging arguments")
    runtime.add_argument(
        "--monitor-interval",
        type=int,
        default=30,
        metavar="INT",
        help="Interval in seconds for sampling CPU and memory usage.",
    )
    runtime.add_argument("--no-color", action="store_true", help="Disable colored terminal output.")
    runtime.add_argument("--quiet", action="store_true", help="Only show warnings and errors on screen.")
    runtime.add_argument("--verbose", action="store_true", help="Show detailed logs, including periodic resource snapshots.")


def add_summarize_arguments(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--hits", metavar="FILE", required=True, help="Input hit TSV from bait2contig search.")
    required.add_argument("--out", metavar="FILE", required=True, help="Output bait-level summary TSV.")

    filtering = parser.add_argument_group("Filtering arguments")
    filtering.add_argument(
        "--min-identity",
        type=float,
        metavar="FLOAT",
        help="Additional identity filter. Default: no additional filter.",
    )
    filtering.add_argument(
        "--min-coverage",
        type=float,
        metavar="FLOAT",
        help="Additional bait coverage filter. Default: no additional filter.",
    )
    filtering.add_argument(
        "--min-aln-length",
        type=int,
        metavar="INT",
        help="Additional alignment length filter. Default: no additional filter.",
    )

    summary = parser.add_argument_group("Summary arguments")
    summary.add_argument("--best-hit", action="store_true", help="Report the best contig for each bait.")
    summary.add_argument("--include-contigs", action="store_true", help="Include contigs anchored by each bait.")
    summary.add_argument("--contig-sep", default=",", metavar="STR", help="Separator for contig lists.")

    output = parser.add_argument_group("Resume and output arguments")
    output.add_argument("--resume", action="store_true", help="Resume if the log shows a previous successful run.")
    output.add_argument("--rerun", action="store_true", help="Force rerun even if output and success log exist.")
    output.add_argument("--force", action="store_true", help="Overwrite existing output.")
    output.add_argument("--gzip", action="store_true", help="Compress output TSV when applicable.")
    output.add_argument("--log", metavar="FILE", help="Log file. Default: <actual_out>.log.")

    runtime = parser.add_argument_group("Runtime and logging arguments")
    runtime.add_argument(
        "--monitor-interval",
        type=int,
        default=30,
        metavar="INT",
        help="Interval in seconds for sampling CPU and memory usage.",
    )
    runtime.add_argument("--no-color", action="store_true", help="Disable colored terminal output.")
    runtime.add_argument("--quiet", action="store_true", help="Only show warnings and errors on screen.")
    runtime.add_argument("--verbose", action="store_true", help="Show detailed logs, including periodic resource snapshots.")


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    color = should_color_help(argv)
    HelpFormatter.use_color = color
    console = Console(stderr=True, color_system="auto" if color else None)
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_usage(sys.stderr)
        console.print("[bold red]ERROR:[/bold red] missing command. Choose one of: search, summarize")
        console.print("Run 'bait2contig --help' for details.")
        return 2
    try:
        return int(args.func(args))
    except (LoggedSearchError, LoggedSummaryError):
        return 1
    except (SearchError, SummaryError) as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        return 1
    except KeyboardInterrupt:
        console.print("[bold red]ERROR:[/bold red] interrupted")
        return 130
    except Exception as exc:
        console.print(f"[bold red]ERROR:[/bold red] {exc}")
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
