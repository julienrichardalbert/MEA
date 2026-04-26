#!/usr/bin/env python3
"""Environment checks for MEApy dependencies."""

from __future__ import annotations

import argparse
import shutil
from typing import Dict, List


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check external tool availability for MEApy.")
    parser.add_argument(
        "--assay",
        choices=["chip", "rna", "wgbs", "all"],
        default="all",
        help="Limit checks to a specific assay profile.",
    )
    return parser


def check_tools(tools: List[str]) -> Dict[str, bool]:
    return {tool: (shutil.which(tool) is not None) for tool in tools}


def main() -> int:
    args = build_parser().parse_args()

    common_tools = ["samtools", "bedtools", "bedGraphToBigWig", "bcftools", "bgzip", "tabix"]
    chip_tools = ["bowtie2", "bowtie2-build", "bwa"]
    rna_tools = ["STAR", "bwa"]
    wgbs_tools = ["bismark", "bismark_methylation_extractor", "bwa"]

    if args.assay == "chip":
        requested = common_tools + chip_tools
    elif args.assay == "rna":
        requested = common_tools + rna_tools
    elif args.assay == "wgbs":
        requested = common_tools + wgbs_tools
    else:
        requested = sorted(set(common_tools + chip_tools + rna_tools + wgbs_tools))

    status = check_tools(requested)
    missing = [tool for tool, ok in status.items() if not ok]

    print("[meapy doctor] tool status")
    for tool in requested:
        marker = "OK" if status[tool] else "MISSING"
        print(f"- {tool}: {marker}")

    if missing:
        print("\n[meapy doctor] missing tools:", ", ".join(missing))
        return 1

    print("\n[meapy doctor] all required tools found.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
