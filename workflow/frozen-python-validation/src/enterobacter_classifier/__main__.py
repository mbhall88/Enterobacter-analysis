"""Public entry point for the unified development assembly workflow."""
from __future__ import annotations

import argparse
import sys

from . import __version__


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="enterobacter-classifier",
        description="Development assembly evidence workflow; no positive calling policy is enabled.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("analyze", add_help=False,
                   help="Run the unified development workflow; use analyze --help for its options")
    if args and args[0] in {"assign", "build-db"}:
        parser.error(
            f"'{args[0]}' is a historical scaffold command. For historical reproduction use "
            f"python -m enterobacter_classifier.legacy_cli {args[0]}. "
            "For the current development workflow use enterobacter-classifier analyze --help."
        )
    if args and args[0] == "analyze":
        # One canonical workflow parser owns all options, defaults and validation.
        # Import lazily so top-level help/version never loads execution adapters.
        from .workflow import main as workflow_main
        workflow_main(args[1:], prog="enterobacter-classifier analyze")
        return
    parser.parse_args(args)


if __name__ == "__main__":
    main()
