#!/usr/bin/env python3
"""Validate expected outputs from a MEApy align run."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

from meapy_common import MeapyError, run_command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate MEApy output artifacts.")
    parser.add_argument("--assay", choices=["chip", "rna", "wgbs"], required=True)
    parser.add_argument("--bam-prefix", required=True, help="Same bam prefix used in align")
    parser.add_argument("--strain1", required=True)
    parser.add_argument("--strain2", required=True)
    parser.add_argument("--tracks-output-dir", required=True)
    return parser


def expect_file(path: Path, missing: List[Path]) -> None:
    if not path.is_file():
        missing.append(path)


def expected_common(run_name: str, out_dir: Path, strain1: str, strain2: str) -> List[Path]:
    return [
        out_dir / f"{run_name}_{strain1}.bedGraph",
        out_dir / f"{run_name}_{strain1}.bw",
        out_dir / f"{run_name}_{strain2}.bedGraph",
        out_dir / f"{run_name}_{strain2}.bw",
        out_dir / f"{run_name}_total.bedGraph",
        out_dir / f"{run_name}_total.bw",
        out_dir / f"{run_name}_{strain1}_projected.bedGraph",
        out_dir / f"{run_name}_{strain1}_projected.bw",
        out_dir / f"{run_name}_{strain2}_projected.bedGraph",
        out_dir / f"{run_name}_{strain2}_projected.bw",
    ]


def expected_wgbs(run_name: str, out_dir: Path, strain1: str, strain2: str, bam_prefix: str) -> List[Path]:
    prefix_dir = Path(bam_prefix).expanduser().resolve().parent
    return [
        out_dir / f"{run_name}_{strain1}_methyl.bedGraph",
        out_dir / f"{run_name}_{strain1}_methyl.bw",
        out_dir / f"{run_name}_{strain2}_methyl.bedGraph",
        out_dir / f"{run_name}_{strain2}_methyl.bw",
        out_dir / f"{run_name}_total_methyl.bedGraph",
        out_dir / f"{run_name}_total_methyl.bw",
        out_dir / f"{run_name}_{strain1}_methyl_LOC.bedGraph",
        out_dir / f"{run_name}_{strain1}_methyl_LOC.bw",
        out_dir / f"{run_name}_{strain2}_methyl_LOC.bedGraph",
        out_dir / f"{run_name}_{strain2}_methyl_LOC.bw",
        out_dir / f"{run_name}_total_methyl_LOC.bedGraph",
        out_dir / f"{run_name}_total_methyl_LOC.bw",
        prefix_dir / f"{run_name}_{strain1}_preProject.CpG_report.txt",
        prefix_dir / f"{run_name}_{strain2}_preProject.CpG_report.txt",
        prefix_dir / f"{run_name}_total.CpG_report.txt",
    ]


def bam_read_count(bam_path: Path) -> int:
    import subprocess

    completed = subprocess.run(
        ["samtools", "view", "-c", str(bam_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise MeapyError(f"Unable to count reads in BAM: {bam_path}")
    return int(completed.stdout.strip())


def main() -> int:
    args = build_parser().parse_args()
    run_name = Path(args.bam_prefix).name
    out_dir = Path(args.tracks_output_dir).expanduser().resolve()
    bam_prefix = Path(args.bam_prefix).expanduser().resolve()

    missing: List[Path] = []
    for expected in expected_common(run_name, out_dir, args.strain1, args.strain2):
        expect_file(expected, missing)
    if args.assay == "wgbs":
        for expected in expected_wgbs(run_name, out_dir, args.strain1, args.strain2, args.bam_prefix):
            expect_file(expected, missing)

    bam1 = Path(f"{bam_prefix}_{args.strain1}.bam")
    bam2 = Path(f"{bam_prefix}_{args.strain2}.bam")
    bam_total = Path(f"{bam_prefix}_total.bam")
    for bam in [bam1, bam2, bam_total]:
        expect_file(bam, missing)

    if missing:
        print("[meapy validate] missing expected files:")
        for path in missing:
            print(f"- {path}")
        return 1

    c1 = bam_read_count(bam1)
    c2 = bam_read_count(bam2)
    ctot = bam_read_count(bam_total)
    print("[meapy validate] core files present")
    print(f"- {bam1.name}: {c1} reads")
    print(f"- {bam2.name}: {c2} reads")
    print(f"- {bam_total.name}: {ctot} reads")
    if ctot <= 0:
        print("[meapy validate] warning: total BAM has zero reads")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
