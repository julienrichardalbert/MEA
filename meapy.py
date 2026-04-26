#!/usr/bin/env python3
"""Flat MEApy CLI in repository root."""

from __future__ import annotations

import argparse

from meapy_align_and_track import main as align_and_track_main
from meapy_create_genome import main as create_genome_main
from meapy_doctor import main as doctor_main
from meapy_phase_vcf import main as phase_vcf_main
from meapy_project import main as project_main
from meapy_validate import main as validate_main


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="meapy",
        description="Simplified Python interface for the MEA pipeline.",
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("phase-vcf", help="Create phased VCF")
    sub.add_parser("create-genome", help="Create in-silico genomes")
    sub.add_parser("align", help="Run align + track as one step")
    sub.add_parser("project", help="Project wig/bedGraph/bed via refmap")
    sub.add_parser("doctor", help="Check external tool dependencies")
    sub.add_parser("validate", help="Validate expected run outputs")
    return parser


def dispatch() -> int:
    parser = build_parser()
    args, remainder = parser.parse_known_args()

    if args.command == "phase-vcf":
        import sys

        sys.argv = ["meapy_phase_vcf.py", *remainder]
        return phase_vcf_main()
    if args.command == "create-genome":
        import sys

        sys.argv = ["meapy_create_genome.py", *remainder]
        return create_genome_main()
    if args.command == "align":
        import sys

        sys.argv = ["meapy_align_and_track.py", *remainder]
        return align_and_track_main()
    if args.command == "project":
        import sys

        sys.argv = ["meapy_project.py", *remainder]
        return project_main()
    if args.command == "doctor":
        import sys

        sys.argv = ["meapy_doctor.py", *remainder]
        return doctor_main()
    if args.command == "validate":
        import sys

        sys.argv = ["meapy_validate.py", *remainder]
        return validate_main()

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(dispatch())
