#!/usr/bin/env python3
"""Phase VCF in native Python (no legacy script dependency)."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Dict, List

from meapy_common import MeapyError, require_path, run_command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create phased VCF from SHAPEIT2 haplotypes and VCF input."
    )
    parser.add_argument("--haps-dir", required=True, help="Directory containing *.haps")
    parser.add_argument("--input-vcf", required=True, help="Input unphased VCF")
    parser.add_argument(
        "--output-prefix", required=True, help="Output path prefix (without extension)"
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        help="Also create bgzip-compressed VCF and tabix index if tools are available.",
    )
    return parser


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def build_haps_dictionary(haps_dir: Path) -> Dict[str, Dict[str, List[str]]]:
    chrom_dict: Dict[str, Dict[str, List[str]]] = {}
    for hap_file in sorted(haps_dir.glob("*.haps")):
        file_name = hap_file.name
        chrom_name = file_name.split(".")[0].replace("chr", "")
        if chrom_name not in chrom_dict:
            chrom_dict[chrom_name] = {}
        with hap_file.open("r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 7:
                    continue
                pos = parts[2]
                parts[0] = chrom_name
                chrom_dict[chrom_name][pos] = parts
    return chrom_dict


def phase_vcf(haps_dir: Path, input_vcf: Path, output_prefix: Path) -> Path:
    haps_by_chrom = build_haps_dictionary(haps_dir)
    output_vcf = Path(f"{output_prefix}.vcf")

    with open_text_auto(input_vcf) as vcf_in, output_vcf.open("w") as vcf_out:
        for raw_line in vcf_in:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                vcf_out.write(f"{line}\n")
                continue
            if line.startswith("#"):
                header = line.split("\t")
                if len(header) < 10:
                    raise MeapyError("Input VCF must contain exactly one sample column.")
                header[9] = "hap1"
                header.append("hap2")
                vcf_out.write("\t".join(header) + "\n")
                continue

            variant = line.split("\t")
            if len(variant) < 10:
                continue
            chrom = variant[0]
            pos = variant[1]
            ref = variant[3]
            alt = variant[4]
            sample_fields = variant[9].split(":")
            gt = sample_fields[0]
            emitted = False

            chrom_haps = haps_by_chrom.get(chrom)
            phased_record = chrom_haps.get(pos) if chrom_haps else None
            if phased_record is not None and phased_record[3] == ref and phased_record[4] == alt:
                hap1_sample = sample_fields.copy()
                hap2_sample = sample_fields.copy()
                hap1_sample[0] = f"{phased_record[5]}/{phased_record[5]}"
                hap2_sample[0] = f"{phased_record[6]}/{phased_record[6]}"
                variant[9] = ":".join(hap1_sample)
                variant.append(":".join(hap2_sample))
                vcf_out.write("\t".join(variant) + "\n")
                emitted = True

            if not emitted and "/" in gt:
                a1, a2 = gt.split("/", 1)
                if a1 == a2:
                    variant.append(variant[9])
                    vcf_out.write("\t".join(variant) + "\n")

    return output_vcf


def compress_and_index_vcf(vcf_path: Path) -> None:
    gz_path = Path(f"{vcf_path}.gz")
    run_command(["bgzip", "-f", str(vcf_path)])
    run_command(["tabix", "-f", "-p", "vcf", str(gz_path)])


def main() -> int:
    args = build_parser().parse_args()
    try:
        haps_dir = require_path(args.haps_dir, kind="dir")
        input_vcf = require_path(args.input_vcf, kind="file")
        output_prefix = Path(args.output_prefix).expanduser().resolve()
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_vcf = phase_vcf(haps_dir, input_vcf, output_prefix)
        if args.compress:
            compress_and_index_vcf(output_vcf)
    except MeapyError as exc:
        print(f"[meapy] error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
