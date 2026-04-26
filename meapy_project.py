#!/usr/bin/env python3
"""Project signal coordinates to reference using a refmap."""

from __future__ import annotations

import argparse
import gzip
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from meapy_common import MeapyError, require_path


@dataclass
class ChromMap:
    ref1_starts: List[int]
    ref2_starts: List[int]
    lengths: List[int]

    def map_ref2_to_ref1(self, pos_ref2: int) -> Optional[int]:
        idx = bisect_right(self.ref2_starts, pos_ref2) - 1
        if idx < 0:
            return None
        start2 = self.ref2_starts[idx]
        interval_len = self.lengths[idx]
        if pos_ref2 >= start2 + interval_len:
            return None
        start1 = self.ref1_starts[idx]
        return start1 + (pos_ref2 - start2)

    def project_interval_ref2_to_ref1(
        self, start_ref2: int, end_ref2: int
    ) -> List[Tuple[int, int]]:
        if end_ref2 <= start_ref2:
            return []
        pieces: List[Tuple[int, int]] = []
        idx = bisect_right(self.ref2_starts, start_ref2) - 1
        if idx < 0:
            idx = 0
        while idx < len(self.ref2_starts):
            seg_start2 = self.ref2_starts[idx]
            seg_end2 = seg_start2 + self.lengths[idx]
            if seg_start2 >= end_ref2:
                break
            overlap_start = max(start_ref2, seg_start2)
            overlap_end = min(end_ref2, seg_end2)
            if overlap_start < overlap_end:
                seg_start1 = self.ref1_starts[idx]
                mapped_start = seg_start1 + (overlap_start - seg_start2)
                mapped_end = mapped_start + (overlap_end - overlap_start)
                pieces.append((mapped_start, mapped_end))
            idx += 1
        return pieces


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def load_refmap(refmap_path: Path) -> Dict[str, ChromMap]:
    refmap: Dict[str, ChromMap] = {}
    current_chr: Optional[str] = None
    first: List[int] = []
    second: List[int] = []
    lengths: List[int] = []

    def commit_current() -> None:
        if current_chr is None:
            return
        refmap[current_chr] = ChromMap(
            ref1_starts=list(first),
            ref2_starts=list(second),
            lengths=list(lengths),
        )

    with refmap_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit_current()
                current_chr = line[1:]
                first = []
                second = []
                lengths = []
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                raise MeapyError(f"Malformed refmap line: {line}")
            first.append(int(parts[0]))
            second.append(int(parts[1]))
            lengths.append(int(parts[2]))

    commit_current()
    return refmap


def iter_wig_records(wig_path: Path) -> Iterator[Tuple[str, int, int, float]]:
    mode: Optional[str] = None
    chrom: Optional[str] = None
    fixed_pos = 0
    fixed_step = 1
    fixed_span = 1
    variable_span = 1

    with open_text(wig_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue

            if line.startswith("fixedStep"):
                mode = "fixed"
                chrom = None
                fixed_pos = 0
                fixed_step = 1
                fixed_span = 1
                for token in line.split():
                    if "=" not in token:
                        continue
                    key, value = token.split("=", 1)
                    if key == "chrom":
                        chrom = value
                    elif key == "start":
                        fixed_pos = int(value) - 1
                    elif key == "step":
                        fixed_step = int(value)
                    elif key == "span":
                        fixed_span = int(value)
                if chrom is None:
                    raise MeapyError(f"fixedStep line missing chrom: {line}")
                continue

            if line.startswith("variableStep"):
                mode = "variable"
                chrom = None
                variable_span = 1
                for token in line.split():
                    if "=" not in token:
                        continue
                    key, value = token.split("=", 1)
                    if key == "chrom":
                        chrom = value
                    elif key == "span":
                        variable_span = int(value)
                if chrom is None:
                    raise MeapyError(f"variableStep line missing chrom: {line}")
                continue

            if mode == "fixed":
                if chrom is None:
                    raise MeapyError("fixedStep data encountered before header")
                value = float(line)
                start = fixed_pos
                end = fixed_pos + fixed_span
                fixed_pos += fixed_step
                yield chrom, start, end, value
                continue

            if mode == "variable":
                if chrom is None:
                    raise MeapyError("variableStep data encountered before header")
                parts = line.split()
                if len(parts) != 2:
                    raise MeapyError(f"Malformed variableStep data line: {line}")
                start = int(parts[0]) - 1
                end = start + variable_span
                value = float(parts[1])
                yield chrom, start, end, value
                continue

            parts = line.split()
            if len(parts) == 4:
                chr_name = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                value = float(parts[3])
                yield chr_name, start, end, value
                continue

            raise MeapyError(f"Unsupported WIG content line: {line}")


def iter_bedgraph_records(path: Path) -> Iterator[Tuple[str, int, int, float]]:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                raise MeapyError(f"Malformed bedGraph line: {line}")
            yield parts[0], int(parts[1]), int(parts[2]), float(parts[3])


def iter_bed_records(path: Path) -> Iterator[Tuple[str, int, int, float]]:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise MeapyError(f"Malformed BED line: {line}")
            value = 1.0
            if len(parts) >= 5:
                try:
                    value = float(parts[4])
                except ValueError:
                    value = 1.0
            yield parts[0], int(parts[1]), int(parts[2]), value


def detect_input_format(input_path: Path) -> str:
    name = input_path.name.lower()
    if name.endswith(".bedgraph") or name.endswith(".bedgraph.gz"):
        return "bedgraph"
    if name.endswith(".bed") or name.endswith(".bed.gz"):
        return "bed"
    if name.endswith(".wig") or name.endswith(".wig.gz"):
        return "wig"
    raise MeapyError(
        "Cannot auto-detect input format. Use --input-format wig|bedgraph|bed."
    )


def iter_input_records(input_path: Path, input_format: str):
    if input_format == "wig":
        return iter_wig_records(input_path)
    if input_format == "bedgraph":
        return iter_bedgraph_records(input_path)
    if input_format == "bed":
        return iter_bed_records(input_path)
    raise MeapyError(f"Unsupported input format: {input_format}")


def project_records(
    records: Iterable[Tuple[str, int, int, float]],
    refmap: Dict[str, ChromMap],
) -> Iterator[Tuple[str, int, int, float]]:
    merged_chr: Optional[str] = None
    merged_start = -1
    merged_end = -1
    merged_value = 0.0

    def flush_if_open():
        nonlocal merged_chr, merged_start, merged_end, merged_value
        if merged_chr is not None and merged_start >= 0:
            yield merged_chr, merged_start, merged_end, merged_value
        merged_chr = None
        merged_start = -1
        merged_end = -1
        merged_value = 0.0

    for chr_name, start, end, value in records:
        chr_map = refmap.get(chr_name)
        if chr_map is None:
            continue
        for mapped_start, mapped_end in chr_map.project_interval_ref2_to_ref1(start, end):
            if (
                merged_chr == chr_name
                and mapped_start == merged_end
                and value == merged_value
            ):
                merged_end = mapped_end
                continue
            yield from flush_if_open()
            merged_chr = chr_name
            merged_start = mapped_start
            merged_end = mapped_end
            merged_value = value

    yield from flush_if_open()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Project wig/bedGraph/bed intervals to reference coordinates."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input track file (wig/wig.gz/bedGraph/bed).",
    )
    parser.add_argument(
        "--input-format",
        choices=["auto", "wig", "bedgraph", "bed"],
        default="auto",
        help="Input format; default auto-detect by extension.",
    )
    parser.add_argument("--input-refmap", required=True, help="Input .refmap file")
    parser.add_argument("--output-bedgraph", required=True, help="Output bedGraph path")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        input_path = require_path(args.input, kind="file")
        refmap_path = require_path(args.input_refmap, kind="file")
        out_path = Path(args.output_bedgraph).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        refmap = load_refmap(refmap_path)
        input_format = (
            detect_input_format(input_path)
            if args.input_format == "auto"
            else args.input_format
        )
        records = iter_input_records(input_path, input_format)
        projected = project_records(records, refmap)

        with out_path.open("w") as output_handle:
            output_handle.write("track type=bedGraph\n")
            for chr_name, start, end, value in projected:
                output_handle.write(f"{chr_name}\t{start}\t{end}\t{value}\n")
    except MeapyError as exc:
        print(f"[meapy] error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
