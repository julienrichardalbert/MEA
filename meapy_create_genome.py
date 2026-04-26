#!/usr/bin/env python3
"""Create in-silico genomes in Python without legacy config."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

from meapy_common import MeapyError, require_path, run_command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create in-silico parental genomes from phased VCF."
    )
    parser.add_argument("--strain1", required=True, help="First strain name")
    parser.add_argument("--strain2", required=True, help="Second strain name")
    parser.add_argument("--output-dir", required=True, help="Output genome directory")
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA file")
    parser.add_argument(
        "--reference-strain",
        default="C57BL6J",
        help="If strain matches this name, reference FASTA is copied directly.",
    )
    parser.add_argument("--phased-vcf", help="Single phased VCF including SNPs + indels")
    parser.add_argument("--phased-snps-vcf", help="Phased SNP VCF")
    parser.add_argument("--phased-indels-vcf", help="Phased indel VCF")
    return parser


def ensure_fasta_index(reference_fasta: Path) -> Path:
    fai = Path(f"{reference_fasta}.fai")
    if not fai.is_file():
        run_command(["samtools", "faidx", str(reference_fasta)])
    return fai


def combine_vcfs(snps_vcf: Path, indels_vcf: Path, output_dir: Path) -> Path:
    combined_path = output_dir / "combined_variants.vcf.gz"
    run_command(
        [
            "bcftools",
            "concat",
            "-a",
            "-Oz",
            "-o",
            str(combined_path),
            str(snps_vcf),
            str(indels_vcf),
        ]
    )
    run_command(["tabix", "-f", "-p", "vcf", str(combined_path)])
    return combined_path


def create_consensus_fasta(reference_fasta: Path, phased_vcf: Path, strain: str, output_fasta: Path) -> None:
    with output_fasta.open("w") as fasta_handle:
        import subprocess

        completed = subprocess.run(
            [
                "bcftools",
                "consensus",
                "-f",
                str(reference_fasta),
                "-s",
                strain,
                str(phased_vcf),
            ],
            check=False,
            stdout=fasta_handle,
        )
    if completed.returncode != 0:
        raise MeapyError(f"bcftools consensus failed for sample {strain}")


def write_concatenated_fasta(fasta1: Path, fasta2: Path, strain1: str, strain2: str, output_fasta: Path) -> None:
    with output_fasta.open("w") as out_handle:
        with fasta1.open("r") as in1:
            for line in in1:
                if line.startswith(">"):
                    out_handle.write(f">{strain1}_{line[1:]}")
                else:
                    out_handle.write(line)
        with fasta2.open("r") as in2:
            for line in in2:
                if line.startswith(">"):
                    out_handle.write(f">{strain2}_{line[1:]}")
                else:
                    out_handle.write(line)


def parse_fai_lengths(fai_path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with fai_path.open("r") as handle:
        for raw_line in handle:
            parts = raw_line.strip().split("\t")
            if len(parts) < 2:
                continue
            lengths[parts[0]] = int(parts[1])
    return lengths


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def load_variants_for_sample(phased_vcf: Path, sample_name: str) -> Dict[str, List[Tuple[int, str, str]]]:
    per_chrom: Dict[str, List[Tuple[int, str, str]]] = {}
    sample_index = -1
    with open_text_auto(phased_vcf) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                if sample_name not in header:
                    raise MeapyError(f"Sample {sample_name} not found in VCF header.")
                sample_index = header.index(sample_name)
                continue
            if sample_index < 0:
                continue
            fields = line.split("\t")
            if len(fields) <= sample_index:
                continue
            chrom = fields[0]
            pos1 = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")
            sample = fields[sample_index].split(":")[0]
            genotype = sample.replace("|", "/").split("/")
            if len(genotype) != 2 or genotype[0] != genotype[1]:
                continue
            allele_idx = int(genotype[0])
            if allele_idx == 0:
                alt = ref
            elif 1 <= allele_idx <= len(alts):
                alt = alts[allele_idx - 1]
            else:
                continue
            per_chrom.setdefault(chrom, []).append((pos1, ref, alt))
    for chrom in per_chrom:
        per_chrom[chrom].sort(key=lambda rec: rec[0])
    return per_chrom


def build_refmap_for_sample(
    chrom_lengths: Dict[str, int],
    sample_variants: Dict[str, List[Tuple[int, str, str]]],
    output_refmap: Path,
) -> None:
    with output_refmap.open("w") as handle:
        for chrom in chrom_lengths:
            chrom_length = chrom_lengths[chrom]
            variant_list = sample_variants.get(chrom, [])
            ref_cursor = 0
            strain_cursor = 0
            intervals: List[Tuple[int, int, int]] = []
            for pos1, ref, alt in variant_list:
                ref_start = pos1 - 1
                if ref_start > ref_cursor:
                    length_before = ref_start - ref_cursor
                    intervals.append((ref_cursor, strain_cursor, length_before))
                    ref_cursor += length_before
                    strain_cursor += length_before
                ref_len = len(ref)
                alt_len = len(alt)
                overlap_len = min(ref_len, alt_len)
                if overlap_len > 0:
                    intervals.append((ref_start, strain_cursor, overlap_len))
                ref_cursor = ref_start + ref_len
                strain_cursor += alt_len
            if ref_cursor < chrom_length:
                tail_len = chrom_length - ref_cursor
                intervals.append((ref_cursor, strain_cursor, tail_len))
            if not intervals:
                intervals.append((0, 0, chrom_length))
            handle.write(f">{chrom}\n")
            for ref_start, strain_start, interval_len in intervals:
                if interval_len > 0:
                    handle.write(f"{ref_start}\t{strain_start}\t{interval_len}\n")


def main() -> int:
    args = build_parser().parse_args()
    try:
        output_dir = Path(args.output_dir).expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        reference_fasta = require_path(args.reference_fasta, kind="file")
        reference_fai = ensure_fasta_index(reference_fasta)

        phased_vcf_path: Path
        if args.phased_vcf:
            phased_vcf_path = require_path(args.phased_vcf, kind="file")
        elif args.phased_snps_vcf and args.phased_indels_vcf:
            snps_vcf = require_path(args.phased_snps_vcf, kind="file")
            indels_vcf = require_path(args.phased_indels_vcf, kind="file")
            phased_vcf_path = combine_vcfs(snps_vcf, indels_vcf, output_dir)
        else:
            raise MeapyError(
                "Provide either --phased-vcf OR both --phased-snps-vcf and --phased-indels-vcf."
            )

        strain1_fasta = output_dir / f"{args.strain1}.fasta"
        strain2_fasta = output_dir / f"{args.strain2}.fasta"

        if args.strain1 == args.reference_strain:
            shutil.copyfile(reference_fasta, strain1_fasta)
            ensure_fasta_index(strain1_fasta)
            chrom_lengths = parse_fai_lengths(reference_fai)
            build_refmap_for_sample(chrom_lengths, {}, output_dir / f"{args.strain1}.fasta.refmap")
        else:
            create_consensus_fasta(reference_fasta, phased_vcf_path, args.strain1, strain1_fasta)
            ensure_fasta_index(strain1_fasta)

        if args.strain2 == args.reference_strain:
            shutil.copyfile(reference_fasta, strain2_fasta)
            ensure_fasta_index(strain2_fasta)
            chrom_lengths = parse_fai_lengths(reference_fai)
            build_refmap_for_sample(chrom_lengths, {}, output_dir / f"{args.strain2}.fasta.refmap")
        else:
            create_consensus_fasta(reference_fasta, phased_vcf_path, args.strain2, strain2_fasta)
            ensure_fasta_index(strain2_fasta)

        chrom_lengths = parse_fai_lengths(reference_fai)
        variants_strain1 = (
            {} if args.strain1 == args.reference_strain else load_variants_for_sample(phased_vcf_path, args.strain1)
        )
        variants_strain2 = (
            {} if args.strain2 == args.reference_strain else load_variants_for_sample(phased_vcf_path, args.strain2)
        )
        build_refmap_for_sample(chrom_lengths, variants_strain1, output_dir / f"{args.strain1}.fasta.refmap")
        build_refmap_for_sample(chrom_lengths, variants_strain2, output_dir / f"{args.strain2}.fasta.refmap")

        concat_fasta = output_dir / f"{args.strain1}_{args.strain2}.fasta"
        write_concatenated_fasta(strain1_fasta, strain2_fasta, args.strain1, args.strain2, concat_fasta)
        ensure_fasta_index(concat_fasta)
        reference_link = output_dir / "reference.fasta"
        if not reference_link.exists():
            reference_link.symlink_to(reference_fasta)
    except MeapyError as exc:
        print(f"[meapy] error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
