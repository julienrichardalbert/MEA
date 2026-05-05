#!/usr/bin/env python3
"""Combined align+track wrapper (replaces separate user-facing steps)."""

from __future__ import annotations

import argparse
import shutil
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from meapy_common import MeapyError, require_path, run_command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run alignment and track creation as one command."
    )
    parser.add_argument(
        "--read-layout",
        choices=["single", "paired"],
        required=True,
        help="Read layout for alignment input",
    )
    parser.add_argument("--reads1", required=True, help="FASTQ/BAM input read 1")
    parser.add_argument("--reads2", help="FASTQ/BAM input read 2 (paired only)")
    parser.add_argument(
        "--genome-input",
        required=True,
        help="Concatenated genome input expected by chosen aligner",
    )
    parser.add_argument("--strain1", required=True, help="First strain name")
    parser.add_argument("--strain2", required=True, help="Second strain name")
    parser.add_argument(
        "--bam-prefix",
        required=True,
        help="Output BAM prefix including desired output directory",
    )
    parser.add_argument(
        "--refmap1",
        help="Refmap for strain1. Defaults to <genome-dir>/<strain1>.fasta.refmap",
    )
    parser.add_argument(
        "--refmap2",
        help="Refmap for strain2. Defaults to <genome-dir>/<strain2>.fasta.refmap",
    )
    parser.add_argument(
        "--tracks-output-dir",
        help="Track output dir. Defaults to directory of --bam-prefix",
    )
    parser.add_argument(
        "--aligner",
        choices=["auto", "bwa", "bowtie2", "star", "tophat2", "bismark", "minimap2"],
        default="auto",
        help="Aligner for native alignment path (default auto by assay).",
    )
    parser.add_argument(
        "--long",
        action="store_true",
        help=(
            "Enable long-read mode for alignment. For --assay rna this selects "
            "a minimap2 + samtools primary/high-confidence filtering path."
        ),
    )
    parser.add_argument(
        "--reference-genome",
        help="Reference FASTA for total alignment.",
    )
    parser.add_argument(
        "--quick-start",
        action="store_true",
        help=(
            "Infer strain FASTAs/refmaps from --genome-input directory and infer "
            "--reference-genome from --reference-fasta or --reference-genome symlink there."
        ),
    )
    parser.add_argument(
        "--reference-fasta",
        help="Reference FASTA path used for quick-start inference.",
    )
    parser.add_argument(
        "--assay",
        choices=["chip", "rna", "wgbs"],
        default="chip",
        help="Assay type for alignment and track generation.",
    )
    parser.add_argument(
        "--chrom-sizes",
        help="Chrom sizes file for track generation.",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=None,
        help=(
            "Minimum MAPQ for filtering BAM before track generation. "
            "If unset, defaults by aligner: STAR=255, Bowtie2=30, Bismark=1."
        ),
    )
    parser.add_argument(
        "--filter-flag",
        type=int,
        default=1540,
        help="SAM flag mask to exclude reads before track generation.",
    )
    parser.add_argument(
        "--min-depth",
        type=int,
        default=1,
        help="Minimum depth used for WGBS methylation bedGraph output.",
    )
    parser.add_argument(
        "--se-extension",
        type=int,
        default=300,
        help="Single-end read extension length for ChIP-style genome coverage (legacy MEA behavior).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Thread count for aligners/indexers and samtools sort where supported.",
    )
    return parser


def default_refmaps(genome_input: str, strain1: str, strain2: str) -> tuple[str, str]:
    genome_path = Path(genome_input).expanduser().resolve()
    genome_dir = genome_path.parent
    refmap1 = str(genome_dir / f"{strain1}.fasta.refmap")
    refmap2 = str(genome_dir / f"{strain2}.fasta.refmap")
    return refmap1, refmap2


def infer_reference_genome(genome_input: str, explicit_reference_fasta: Optional[str]) -> Optional[str]:
    if explicit_reference_fasta:
        return str(Path(explicit_reference_fasta).expanduser().resolve())
    genome_path = Path(genome_input).expanduser().resolve()
    genome_dir = genome_path.parent
    candidates = [
        genome_dir / "reference.fasta",
        genome_dir / "reference.fa",
        genome_dir / "reference_genome.fasta",
        genome_dir / "reference_genome.fa",
        genome_dir / "reference.fasta.symlink",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    return None


def ensure_chrom_sizes(
    chrom_sizes_arg: Optional[str],
    quick_start: bool,
    reference_fasta: Path,
    genome_input: str,
) -> str:
    if chrom_sizes_arg:
        return str(Path(chrom_sizes_arg).expanduser().resolve())
    if not quick_start:
        raise MeapyError("--chrom-sizes is required (or use --quick-start).")

    genome_dir = Path(genome_input).expanduser().resolve().parent
    for candidate in [genome_dir / "reference.chrom.sizes", genome_dir / "chrom.sizes"]:
        if candidate.is_file():
            return str(candidate.resolve())

    fai = Path(f"{reference_fasta}.fai")
    if not fai.is_file():
        run_command(["samtools", "faidx", str(reference_fasta)])
    chrom_sizes_out = genome_dir / "reference.chrom.sizes"
    with fai.open("r") as in_handle, chrom_sizes_out.open("w") as out_handle:
        for raw_line in in_handle:
            parts = raw_line.strip().split("\t")
            if len(parts) >= 2:
                out_handle.write(f"{parts[0]}\t{parts[1]}\n")
    return str(chrom_sizes_out.resolve())


def run_shell_pipeline(command: str) -> None:
    completed = subprocess.run(
        ["bash", "-c", command],
        cwd=str(Path(__file__).resolve().parent),
        check=False,
    )
    if completed.returncode != 0:
        raise MeapyError(f"Pipeline command failed: {command}")


def run_command_in_dir(cmd: List[str], working_dir: Path, env_overrides: Optional[dict[str, str]] = None) -> None:
    if not cmd:
        raise MeapyError("Empty command provided.")
    env = None
    if env_overrides:
        env = dict(**subprocess.os.environ)
        env.update(env_overrides)
    completed = subprocess.run(
        cmd,
        cwd=str(working_dir),
        env=env,
        check=False,
    )
    if completed.returncode != 0:
        raise MeapyError(f"Command failed with code {completed.returncode}: {cmd[0]}")


def ensure_bwa_index(fasta_path: Path) -> None:
    bwt_path = Path(f"{fasta_path}.bwt")
    if bwt_path.is_file():
        print(f"[meapy] using existing BWA index for {fasta_path}")
        return
    print(f"[meapy] BWA index missing; building beside FASTA: {fasta_path}")
    run_command(["bwa", "index", str(fasta_path)])


def ensure_bowtie2_index(fasta_path: Path, index_prefix: Path, threads: int) -> None:
    if (Path(f"{index_prefix}.1.bt2").is_file()) or (Path(f"{index_prefix}.1.bt2l").is_file()):
        print(f"[meapy] using existing Bowtie2 index prefix {index_prefix}")
        return
    print(f"[meapy] Bowtie2 index missing; building beside FASTA with prefix {index_prefix}")
    run_command(
        ["bowtie2-build", "--threads", str(max(1, threads)), str(fasta_path), str(index_prefix)]
    )


def ensure_star_index(fasta_path: Path, index_dir: Path, threads: int) -> None:
    required = [
        index_dir / "Genome",
        index_dir / "SA",
        index_dir / "SAindex",
        index_dir / "chrLength.txt",
        index_dir / "genomeParameters.txt",
    ]
    if all(p.is_file() for p in required):
        print(f"[meapy] using existing STAR genomeDir {index_dir}")
        return
    print(f"[meapy] STAR genomeDir missing/incomplete; building beside FASTA at {index_dir}")
    index_dir.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            str(index_dir),
            "--genomeFastaFiles",
            str(fasta_path),
            "--runThreadN",
            str(max(1, threads)),
        ]
    )


def align_bwa_sorted_bam(
    fasta_path: Path,
    reads1: Path,
    reads2: Optional[Path],
    output_bam: Path,
    threads: int,
) -> None:
    ensure_bwa_index(fasta_path)
    if reads2 is None:
        command = (
            f"bwa mem -t {int(max(1, threads))} {shlex.quote(str(fasta_path))} "
            f"{shlex.quote(str(reads1))} | "
            f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
        )
    else:
        command = (
            f"bwa mem -t {int(max(1, threads))} {shlex.quote(str(fasta_path))} {shlex.quote(str(reads1))} "
            f"{shlex.quote(str(reads2))} | "
            f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
        )
    run_shell_pipeline(command)
    run_command(["samtools", "index", str(output_bam)])


def align_bowtie2_sorted_bam(
    fasta_path: Path,
    reads1: Path,
    reads2: Optional[Path],
    output_bam: Path,
    threads: int,
) -> None:
    index_prefix = fasta_path.parent / f"{fasta_path.stem}.bowtie2_index"
    ensure_bowtie2_index(fasta_path, index_prefix, threads=threads)
    if reads2 is None:
        command = (
            f"bowtie2 -p {int(max(1, threads))} "
            f"-x {shlex.quote(str(index_prefix))} "
            f"-U {shlex.quote(str(reads1))} | "
            f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
        )
    else:
        command = (
            f"bowtie2 -p {int(max(1, threads))} "
            f"-x {shlex.quote(str(index_prefix))} "
            f"-1 {shlex.quote(str(reads1))} -2 {shlex.quote(str(reads2))} | "
            f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
        )
    run_shell_pipeline(command)
    run_command(["samtools", "index", str(output_bam)])


def align_star_sorted_bam(
    fasta_path: Path,
    reads1: Path,
    reads2: Optional[Path],
    output_bam: Path,
    threads: int,
) -> None:
    index_dir = fasta_path.parent / f"{fasta_path.stem}.star_index"
    ensure_star_index(fasta_path, index_dir, threads=threads)
    star_prefix = str(output_bam.parent / f"{output_bam.stem}_star_")
    cmd = [
        "STAR",
        "--runMode",
        "alignReads",
        "--genomeDir",
        str(index_dir),
        "--outFileNamePrefix",
        star_prefix,
        "--readFilesIn",
        str(reads1),
    ]
    if reads2 is not None:
        cmd.append(str(reads2))
    cmd.extend(
        [
            "--runThreadN",
            str(max(1, threads)),
            "--outSAMtype",
            "SAM",
        ]
    )
    reads_are_gz = str(reads1).endswith(".gz") or (
        reads2 is not None and str(reads2).endswith(".gz")
    )
    if reads_are_gz:
        cmd.extend(["--readFilesCommand", "gunzip", "-c"])
    run_command(cmd)
    star_final = Path(f"{star_prefix}Log.final.out")
    if star_final.is_file():
        input_reads = None
        with star_final.open("r") as handle:
            for raw_line in handle:
                if "Number of input reads" in raw_line:
                    pieces = raw_line.strip().split("|")
                    if len(pieces) >= 2:
                        try:
                            input_reads = int(pieces[-1].strip())
                        except ValueError:
                            input_reads = None
                    break
        if input_reads == 0:
            raise MeapyError(
                "STAR reported zero input reads. This STAR runtime appears broken for read parsing in "
                f"this environment (log: {star_final})."
            )
    star_sam = Path(f"{star_prefix}Aligned.out.sam")
    if not star_sam.is_file():
        raise MeapyError(f"STAR SAM output missing: {star_sam}")
    run_shell_pipeline(
        f"samtools view -bShu {shlex.quote(str(star_sam))} | "
        f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
    )
    run_command(["samtools", "index", str(output_bam)])


def align_tophat2_sorted_bam(
    fasta_path: Path,
    reads1: Path,
    reads2: Optional[Path],
    output_bam: Path,
    threads: int,
) -> None:
    python2_bin = shutil.which("python2")
    tophat_script = shutil.which("tophat")
    if not python2_bin or not tophat_script:
        raise MeapyError(
            "TopHat2 requires a Python 2 runtime. This environment does not provide `python2`/`tophat` "
            "in a compatible form."
        )
    index_prefix = fasta_path.parent / f"{fasta_path.stem}.bowtie2_index"
    ensure_bowtie2_index(fasta_path, index_prefix, threads=threads)
    tophat_out = output_bam.parent / f"{output_bam.stem}_tophat2"
    tophat_out.mkdir(parents=True, exist_ok=True)
    cmd = [
        python2_bin,
        tophat_script,
        "--read-mismatches",
        "0",
        "--read-gap-length",
        "0",
        "--read-edit-dist",
        "0",
        "--no-sort-bam",
        "--no-convert-bam",
        "-o",
        str(tophat_out),
        "-p",
        str(max(1, threads)),
        str(index_prefix),
        str(reads1),
    ]
    if reads2 is not None:
        cmd.append(str(reads2))
    run_command(cmd)
    accepted_hits = tophat_out / "accepted_hits.bam"
    if not accepted_hits.is_file():
        raise MeapyError(f"TopHat2 BAM output missing: {accepted_hits}")
    run_command(
        ["samtools", "sort", "-@", str(max(1, threads)), "-o", str(output_bam), str(accepted_hits)]
    )
    run_command(["samtools", "index", str(output_bam)])


def align_minimap2_rna_sorted_bam(
    fasta_path: Path,
    reads1: Path,
    reads2: Optional[Path],
    output_bam: Path,
    min_mapq: int,
    threads: int,
) -> None:
    if reads2 is not None:
        raise MeapyError("Long-read RNA minimap2 mode only supports --read-layout single.")
    cmd = (
        f"minimap2 -t {int(max(1, threads))} -ax splice {shlex.quote(str(fasta_path))} {shlex.quote(str(reads1))} | "
        f"samtools view -h -F 0x900 -q {int(min_mapq)} - | "
        f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
    )
    run_shell_pipeline(cmd)
    run_command(["samtools", "index", str(output_bam)])


def run_python_alignment(
    reads1: Path,
    reads2: Optional[Path],
    strain1_name: str,
    strain2_name: str,
    pseudogenome_fasta: Path,
    reference_fasta: Path,
    bam_prefix: str,
    aligner: str,
    long_mode: bool = False,
    threads: int = 4,
    split_min_mapq: int = 0,
) -> None:
    concat_bam = Path(f"{bam_prefix}_{strain1_name}_{strain2_name}.bam").expanduser().resolve()
    bam1 = Path(f"{bam_prefix}_{strain1_name}.bam").expanduser().resolve()
    bam2 = Path(f"{bam_prefix}_{strain2_name}.bam").expanduser().resolve()
    bam_total = Path(f"{bam_prefix}_total.bam").expanduser().resolve()
    bam1.parent.mkdir(parents=True, exist_ok=True)

    if long_mode:
        if aligner != "minimap2":
            raise MeapyError("Long-read mode currently supports aligner minimap2 only.")
        align_minimap2_rna_sorted_bam(
            pseudogenome_fasta, reads1, reads2, concat_bam, min_mapq=20, threads=threads
        )
        align_minimap2_rna_sorted_bam(
            reference_fasta, reads1, reads2, bam_total, min_mapq=20, threads=threads
        )
    elif aligner == "bwa":
        align_bwa_sorted_bam(pseudogenome_fasta, reads1, reads2, concat_bam, threads=threads)
        align_bwa_sorted_bam(reference_fasta, reads1, reads2, bam_total, threads=threads)
    elif aligner == "bowtie2":
        align_bowtie2_sorted_bam(
            pseudogenome_fasta,
            reads1,
            reads2,
            concat_bam,
            threads=threads,
        )
        align_bowtie2_sorted_bam(
            reference_fasta,
            reads1,
            reads2,
            bam_total,
            threads=threads,
        )
    elif aligner == "star":
        align_star_sorted_bam(pseudogenome_fasta, reads1, reads2, concat_bam, threads=threads)
        align_star_sorted_bam(reference_fasta, reads1, reads2, bam_total, threads=threads)
    elif aligner == "tophat2":
        align_tophat2_sorted_bam(
            pseudogenome_fasta, reads1, reads2, concat_bam, threads=threads
        )
        align_tophat2_sorted_bam(reference_fasta, reads1, reads2, bam_total, threads=threads)
    else:
        raise MeapyError(f"Unsupported Python aligner: {aligner}")

    split_allelic_bams_from_concat(
        concat_bam=concat_bam,
        strain1=strain1_name,
        strain2=strain2_name,
        reference_fasta=reference_fasta,
        output_bam1=bam1,
        output_bam2=bam2,
        min_mapq=split_min_mapq,
        exact_mapq=None,
    )


def split_allelic_bams_from_concat(
    concat_bam: Path,
    strain1: str,
    strain2: str,
    reference_fasta: Path,
    output_bam1: Path,
    output_bam2: Path,
    min_mapq: int = 0,
    exact_mapq: Optional[int] = None,
) -> None:
    try:
        import pysam
    except ModuleNotFoundError as exc:
        raise MeapyError(
            "Missing required Python package `pysam` for align/split workflow. "
            "Install it in your active environment (for example: conda install -n meapy -c conda-forge pysam)."
        ) from exc

    reference_fai = Path(f"{reference_fasta}.fai")
    if not reference_fai.is_file():
        run_command(["samtools", "faidx", str(reference_fasta)])

    ref_lengths: List[Tuple[str, int]] = []
    with reference_fai.open("r") as fai_handle:
        for raw_line in fai_handle:
            parts = raw_line.strip().split("\t")
            if len(parts) >= 2:
                ref_lengths.append((parts[0], int(parts[1])))
    if not ref_lengths:
        raise MeapyError(f"Reference FASTA index has no contigs: {reference_fai}")

    out_header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": l} for c, l in ref_lengths]}

    name_sorted_bam = concat_bam.with_suffix(".qname.tmp.bam")
    run_command(["samtools", "sort", "-n", "-o", str(name_sorted_bam), str(concat_bam)])

    unsorted_bam1 = output_bam1.with_suffix(".unsorted.bam")
    unsorted_bam2 = output_bam2.with_suffix(".unsorted.bam")
    prefixes = {strain1: f"{strain1}_", strain2: f"{strain2}_"}

    def read_passes(read: "pysam.AlignedSegment") -> bool:
        if read.is_unmapped:
            return False
        if read.is_secondary or read.is_supplementary:
            return False
        if exact_mapq is not None:
            return read.mapping_quality == exact_mapq
        return read.mapping_quality >= min_mapq

    def remap_read_to_reference(
        read: "pysam.AlignedSegment",
        source_bam: "pysam.AlignmentFile",
        out_bam: "pysam.AlignmentFile",
        strain: str,
    ) -> Optional["pysam.AlignedSegment"]:
        prefix = prefixes[strain]
        ref_name = source_bam.get_reference_name(read.reference_id)
        if ref_name is None or not ref_name.startswith(prefix):
            return None
        target_ref = ref_name[len(prefix) :]
        target_tid = out_bam.get_tid(target_ref)
        if target_tid < 0:
            return None
        read.reference_id = target_tid
        if read.next_reference_id >= 0:
            next_ref_name = source_bam.get_reference_name(read.next_reference_id)
            if next_ref_name and next_ref_name.startswith(prefix):
                target_next = next_ref_name[len(prefix) :]
                next_tid = out_bam.get_tid(target_next)
                read.next_reference_id = next_tid if next_tid >= 0 else -1
            else:
                read.next_reference_id = -1
        return read

    with pysam.AlignmentFile(str(name_sorted_bam), "rb") as in_bam, pysam.AlignmentFile(
        str(unsorted_bam1), "wb", header=out_header
    ) as out1, pysam.AlignmentFile(str(unsorted_bam2), "wb", header=out_header) as out2:
        current_qname: Optional[str] = None
        current_group: List["pysam.AlignedSegment"] = []

        def flush_group(group: List["pysam.AlignedSegment"]) -> None:
            if not group:
                return
            kept: List[tuple[str, "pysam.AlignedSegment"]] = []
            strains_seen = set()
            for rec in group:
                if not read_passes(rec):
                    continue
                ref_name = in_bam.get_reference_name(rec.reference_id)
                if ref_name is None:
                    continue
                if ref_name.startswith(prefixes[strain1]):
                    strains_seen.add(strain1)
                    kept.append((strain1, rec))
                elif ref_name.startswith(prefixes[strain2]):
                    strains_seen.add(strain2)
                    kept.append((strain2, rec))
            # Critical: drop ambiguous reads that map to both haplotypes.
            if len(strains_seen) != 1:
                return
            chosen = next(iter(strains_seen))
            for strain_name, rec in kept:
                if strain_name != chosen:
                    continue
                if chosen == strain1:
                    remapped = remap_read_to_reference(rec, in_bam, out1, strain1)
                    if remapped is not None:
                        out1.write(remapped)
                else:
                    remapped = remap_read_to_reference(rec, in_bam, out2, strain2)
                    if remapped is not None:
                        out2.write(remapped)

        for read in in_bam.fetch(until_eof=True):
            qname = read.query_name
            if current_qname is None:
                current_qname = qname
            if qname != current_qname:
                flush_group(current_group)
                current_group = []
                current_qname = qname
            current_group.append(read)
        flush_group(current_group)

    run_command(["samtools", "sort", "-o", str(output_bam1), str(unsorted_bam1)])
    run_command(["samtools", "index", str(output_bam1)])
    run_command(["samtools", "sort", "-o", str(output_bam2), str(unsorted_bam2)])
    run_command(["samtools", "index", str(output_bam2)])
    unsorted_bam1.unlink(missing_ok=True)
    unsorted_bam2.unlink(missing_ok=True)
    name_sorted_bam.unlink(missing_ok=True)


def run_bismark_alignment(
    genome_folder: Path,
    reads1: Path,
    reads2: Optional[Path],
    output_bam: Path,
    output_name: str,
    threads: int,
) -> None:
    output_dir = output_bam.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    bismark_index_dir = genome_folder / "Bisulfite_Genome"
    if bismark_index_dir.is_dir():
        print(f"[meapy] using existing Bismark index in {genome_folder}")
    else:
        raise MeapyError(
            "Missing Bismark genome preparation output in "
            f"{genome_folder} (expected directory {bismark_index_dir}). "
            f"Build it with: bismark_genome_preparation --bowtie2 {genome_folder}"
        )
    bowtie2_bin = shutil.which("bowtie2")
    samtools_bin = shutil.which("samtools")
    if not bowtie2_bin or not samtools_bin:
        raise MeapyError("Bismark mode requires bowtie2 and samtools on PATH.")
    bowtie2_dir = str(Path(bowtie2_bin).parent)
    samtools_dir = str(Path(samtools_bin).parent)

    cmd = [
        "bismark",
        "--bowtie2",
        "--sam",
        "-p",
        str(max(1, threads)),
        "--path_to_bowtie",
        bowtie2_dir,
        "--samtools_path",
        samtools_dir,
        "--temp_dir",
        str(output_dir),
        "--basename",
        output_name,
        "-o",
        str(output_dir),
        str(genome_folder),
    ]
    if reads2 is None:
        cmd.append(str(reads1))
    else:
        cmd.extend(["-1", str(reads1), "-2", str(reads2)])
    # Run from output directory so Bismark temp/intermediate files are created there.
    bismark_env = {"TMPDIR": str(output_dir), "TMP": str(output_dir), "TEMP": str(output_dir)}
    run_command_in_dir(cmd, output_dir, env_overrides=bismark_env)
    pe_bam_path = output_dir / f"{output_name}_pe.bam"
    bam_path = output_dir / f"{output_name}.bam"
    pe_sam_path = output_dir / f"{output_name}_pe.sam"
    sam_path = output_dir / f"{output_name}.sam"
    if pe_bam_path.is_file():
        run_command(
            ["samtools", "sort", "-@", str(max(1, threads)), "-o", str(output_bam), str(pe_bam_path)]
        )
    elif bam_path.is_file():
        run_command(
            ["samtools", "sort", "-@", str(max(1, threads)), "-o", str(output_bam), str(bam_path)]
        )
    elif pe_sam_path.is_file():
        run_shell_pipeline(
            f"samtools view -bShu {shlex.quote(str(pe_sam_path))} | "
            f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
        )
    elif sam_path.is_file():
        run_shell_pipeline(
            f"samtools view -bShu {shlex.quote(str(sam_path))} | "
            f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
        )
    else:
        # Parallel runs can force default Bismark file naming.
        fallback_pe_bam = sorted(output_dir.glob("*_pe.bam"), key=lambda p: p.stat().st_mtime, reverse=True)
        fallback_bam = sorted(output_dir.glob("*.bam"), key=lambda p: p.stat().st_mtime, reverse=True)
        fallback_pe_sam = sorted(output_dir.glob("*_pe.sam"), key=lambda p: p.stat().st_mtime, reverse=True)
        fallback_sam = sorted(output_dir.glob("*.sam"), key=lambda p: p.stat().st_mtime, reverse=True)
        if fallback_pe_bam:
            run_command(
                [
                    "samtools",
                    "sort",
                    "-@",
                    str(max(1, threads)),
                    "-o",
                    str(output_bam),
                    str(fallback_pe_bam[0]),
                ]
            )
        elif fallback_bam:
            run_command(
                [
                    "samtools",
                    "sort",
                    "-@",
                    str(max(1, threads)),
                    "-o",
                    str(output_bam),
                    str(fallback_bam[0]),
                ]
            )
        elif fallback_pe_sam:
            run_shell_pipeline(
                f"samtools view -bShu {shlex.quote(str(fallback_pe_sam[0]))} | "
                f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
            )
        elif fallback_sam:
            run_shell_pipeline(
                f"samtools view -bShu {shlex.quote(str(fallback_sam[0]))} | "
                f"samtools sort -@ {int(max(1, threads))} -o {shlex.quote(str(output_bam))}"
            )
        else:
            raise MeapyError(
                f"Bismark alignment output not found for {output_name} "
                f"(checked: {pe_bam_path.name}, {bam_path.name}, {pe_sam_path.name}, {sam_path.name})"
            )
    run_command(["samtools", "index", str(output_bam)])


def run_bismark_methyl_extractor(
    input_bam: Path,
    genome_folder: Path,
    output_dir: Path,
    output_prefix: str,
    is_paired: bool,
    threads: int,
) -> Path:
    extractor_input = input_bam
    if is_paired:
        # Bismark methylation extractor expects paired reads to be adjacent.
        # Ensure query-name sort even if upstream BAM is coordinate-sorted.
        extractor_input = input_bam.with_suffix(".qname.bam")
        run_command(
            [
                "samtools",
                "sort",
                "-@",
                str(max(1, threads)),
                "-n",
                "-o",
                str(extractor_input),
                str(input_bam),
            ]
        )

    samtools_bin = shutil.which("samtools")
    if not samtools_bin:
        raise MeapyError("Bismark methylation extractor requires samtools on PATH.")
    samtools_dir = str(Path(samtools_bin).parent)

    cmd = [
        "bismark_methylation_extractor",
        "-p" if is_paired else "-s",
        "--comprehensive",
        "--cytosine_report",
        "--samtools_path",
        samtools_dir,
        "-o",
        str(output_dir),
        "--genome_folder",
        str(genome_folder),
    ]
    cmd.append(str(extractor_input))
    # Run from output directory so Bismark temp/intermediate files are created there.
    bismark_env = {"TMPDIR": str(output_dir), "TMP": str(output_dir), "TEMP": str(output_dir)}
    run_command_in_dir(cmd, output_dir, env_overrides=bismark_env)
    expected_candidates = [
        output_dir / f"{extractor_input.stem}.CpG_report.txt",
        output_dir / f"{input_bam.stem}.CpG_report.txt",
    ]
    cpg_report = None
    for candidate in expected_candidates:
        if candidate.is_file():
            cpg_report = candidate
            break
    if cpg_report is None:
        # Fallback: accept any CpG_report produced for this run prefix.
        glob_candidates = sorted(output_dir.glob(f"{output_prefix}*.CpG_report.txt"))
        if glob_candidates:
            cpg_report = glob_candidates[0]
    if cpg_report is None:
        checked = ", ".join(str(p) for p in expected_candidates)
        raise MeapyError(f"CpG report not created. Checked: {checked}")
    renamed = output_dir / f"{output_prefix}.CpG_report.txt"
    if renamed != cpg_report:
        cpg_report.rename(renamed)
    return renamed


def split_cpg_by_strain(
    combined_cpg_report: Path,
    strain1: str,
    strain2: str,
    output1: Path,
    output2: Path,
) -> None:
    with combined_cpg_report.open("r") as in_handle, output1.open("w") as out1, output2.open(
        "w"
    ) as out2:
        for raw_line in in_handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 1:
                continue
            chrom = parts[0]
            if chrom.startswith(f"{strain1}_"):
                parts[0] = chrom[len(strain1) + 1 :]
                out1.write("\t".join(parts) + "\n")
            elif chrom.startswith(f"{strain2}_"):
                parts[0] = chrom[len(strain2) + 1 :]
                out2.write("\t".join(parts) + "\n")


def sort_bedgraph_file(path: Path) -> None:
    rows: List[tuple[str, int, int, str]] = []
    with path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("track"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    rows.sort(key=lambda row: (row[0], row[1], row[2]))
    with path.open("w") as handle:
        for chrom, start, end, value in rows:
            handle.write(f"{chrom}\t{start}\t{end}\t{value}\n")


def ensure_nonempty_bedgraph(path: Path, chrom_sizes: str) -> None:
    if path.stat().st_size > 0:
        return
    chrom_file = Path(chrom_sizes).expanduser().resolve()
    with chrom_file.open("r") as handle:
        for raw_line in handle:
            parts = raw_line.strip().split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            length = int(parts[1])
            if length > 0:
                with path.open("w") as out:
                    out.write(f"{chrom}\t0\t1\t0\n")
                return
    raise MeapyError(f"Cannot build placeholder bedGraph from chrom sizes: {chrom_file}")


def project_bedgraph(input_bedgraph: Path, input_refmap: Path, output_bedgraph: Path) -> None:
    run_command(
        [
            sys.executable,
            str((Path(__file__).resolve().parent / "meapy_project.py")),
            "--input",
            str(input_bedgraph),
            "--input-format",
            "bedgraph",
            "--input-refmap",
            str(input_refmap),
            "--output-bedgraph",
            str(output_bedgraph),
        ]
    )


def count_bam_alignments(input_bam: Path) -> int:
    completed = subprocess.run(
        ["samtools", "view", "-c", str(input_bam)],
        cwd=str(Path(__file__).resolve().parent),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise MeapyError(f"Failed to count alignments in BAM: {input_bam}")
    try:
        return int(completed.stdout.strip() or "0")
    except ValueError as exc:
        raise MeapyError(f"Unexpected samtools count output for {input_bam}: {completed.stdout!r}") from exc


def write_rpm_scaled_bedgraph(input_bedgraph: Path, output_bedgraph: Path, rpm_scale: float) -> None:
    with input_bedgraph.open("r") as in_handle, output_bedgraph.open("w") as out_handle:
        for raw_line in in_handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("track") or line.startswith("#"):
                out_handle.write(raw_line if raw_line.endswith("\n") else raw_line + "\n")
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                value = float(parts[3]) * rpm_scale
            except ValueError:
                continue
            out_handle.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{value}\n")


def make_tracks_from_bam(
    bam_prefix: str,
    strain_name: str,
    chrom_sizes: str,
    output_dir: str,
    min_mapq: int,
    filter_flag: int,
    split_reads: bool,
    single_end_extension: Optional[int],
) -> Path:
    output_prefix = Path(output_dir).expanduser().resolve() / (
        f"{Path(bam_prefix).name}_{strain_name}"
    )
    input_bam = Path(f"{bam_prefix}_{strain_name}.bam").expanduser().resolve()
    if not input_bam.is_file():
        raise MeapyError(f"Expected BAM not found: {input_bam}")
    filtered_bam = Path(f"{output_prefix}_F{filter_flag}_q{min_mapq}.bam")
    bedgraph_path = Path(f"{output_prefix}.bedGraph")
    bw_path = Path(f"{output_prefix}.bw")

    run_command(
        [
            "samtools",
            "view",
            "-bh",
            "-F",
            str(filter_flag),
            "-q",
            str(min_mapq),
            str(input_bam),
            "-o",
            str(filtered_bam),
        ]
    )
    genomecov_cmd = [
        "bedtools",
        "genomecov",
        "-bg",
        "-ibam",
        str(filtered_bam),
        "-g",
        str(Path(chrom_sizes).expanduser().resolve()),
    ]
    if single_end_extension is not None and single_end_extension > 0:
        genomecov_cmd.extend(["-fs", str(single_end_extension)])
    if split_reads:
        genomecov_cmd.insert(3, "-split")
    with bedgraph_path.open("w") as bedgraph_handle:
        completed = subprocess.run(
            genomecov_cmd,
            cwd=str(Path(__file__).resolve().parent),
            check=False,
            stdout=bedgraph_handle,
        )
    if completed.returncode != 0:
        raise MeapyError("bedtools genomecov failed.")
    sort_bedgraph_file(bedgraph_path)
    ensure_nonempty_bedgraph(bedgraph_path, chrom_sizes)
    run_command(
        [
            "bedGraphToBigWig",
            str(bedgraph_path),
            str(Path(chrom_sizes).expanduser().resolve()),
            str(bw_path),
        ]
    )
    return bedgraph_path


def generate_tracks_python(
    bam_prefix: str,
    strain1: str,
    strain2: str,
    refmap1: str,
    refmap2: str,
    tracks_output_dir: str,
    chrom_sizes: str,
    min_mapq: int,
    filter_flag: int,
    assay: str,
    read_layout: str,
    aligner: str,
    se_extension: int,
) -> None:
    split_reads = assay == "rna"
    single_end_extension: Optional[int] = None
    if (
        read_layout == "single"
        and assay != "rna"
        and aligner in {"bwa", "bowtie2"}
    ):
        single_end_extension = se_extension
    Path(tracks_output_dir).expanduser().resolve().mkdir(parents=True, exist_ok=True)

    bed1 = make_tracks_from_bam(
        bam_prefix,
        strain1,
        chrom_sizes,
        tracks_output_dir,
        min_mapq,
        filter_flag,
        split_reads,
        single_end_extension,
    )
    bed2 = make_tracks_from_bam(
        bam_prefix,
        strain2,
        chrom_sizes,
        tracks_output_dir,
        min_mapq,
        filter_flag,
        split_reads,
        single_end_extension,
    )
    total_bed = make_tracks_from_bam(
        bam_prefix,
        "total",
        chrom_sizes,
        tracks_output_dir,
        min_mapq,
        filter_flag,
        split_reads,
        single_end_extension,
    )

    projected1 = Path(tracks_output_dir).expanduser().resolve() / f"{Path(bam_prefix).name}_{strain1}_projected.bedGraph"
    projected2 = Path(tracks_output_dir).expanduser().resolve() / f"{Path(bam_prefix).name}_{strain2}_projected.bedGraph"
    project_bedgraph(bed1, Path(refmap1).expanduser().resolve(), projected1)
    project_bedgraph(bed2, Path(refmap2).expanduser().resolve(), projected2)

    run_command(
        [
            "bedGraphToBigWig",
            str(projected1),
            str(Path(chrom_sizes).expanduser().resolve()),
            str(projected1).replace(".bedGraph", ".bw"),
        ]
    )
    run_command(
        [
            "bedGraphToBigWig",
            str(projected2),
            str(Path(chrom_sizes).expanduser().resolve()),
            str(projected2).replace(".bedGraph", ".bw"),
        ]
    )
    output_root = Path(tracks_output_dir).expanduser().resolve()
    total_filtered_bam = output_root / (
        f"{Path(bam_prefix).name}_total_F{filter_flag}_q{min_mapq}.bam"
    )
    total_alignments = count_bam_alignments(total_filtered_bam)
    rpm_scale = 0.0 if total_alignments <= 0 else 1_000_000.0 / float(total_alignments)
    print(f"[meapy] RPM scaling factor from total filtered alignments: {rpm_scale}")

    strain1_rpm_bed = output_root / f"{Path(bam_prefix).name}_{strain1}_RPM.bedGraph"
    strain2_rpm_bed = output_root / f"{Path(bam_prefix).name}_{strain2}_RPM.bedGraph"
    total_rpm_bed = output_root / (
        f"{Path(bam_prefix).name}_total_F{filter_flag}_q{min_mapq}_RPM.bedGraph"
    )
    write_rpm_scaled_bedgraph(projected1, strain1_rpm_bed, rpm_scale)
    write_rpm_scaled_bedgraph(projected2, strain2_rpm_bed, rpm_scale)
    write_rpm_scaled_bedgraph(total_bed, total_rpm_bed, rpm_scale)
    sort_bedgraph_file(strain1_rpm_bed)
    sort_bedgraph_file(strain2_rpm_bed)
    sort_bedgraph_file(total_rpm_bed)
    ensure_nonempty_bedgraph(strain1_rpm_bed, chrom_sizes)
    ensure_nonempty_bedgraph(strain2_rpm_bed, chrom_sizes)
    ensure_nonempty_bedgraph(total_rpm_bed, chrom_sizes)
    run_command(
        [
            "bedGraphToBigWig",
            str(strain1_rpm_bed),
            str(Path(chrom_sizes).expanduser().resolve()),
            str(strain1_rpm_bed).replace(".bedGraph", ".bw"),
        ]
    )
    run_command(
        [
            "bedGraphToBigWig",
            str(strain2_rpm_bed),
            str(Path(chrom_sizes).expanduser().resolve()),
            str(strain2_rpm_bed).replace(".bedGraph", ".bw"),
        ]
    )
    run_command(
        [
            "bedGraphToBigWig",
            str(total_rpm_bed),
            str(Path(chrom_sizes).expanduser().resolve()),
            str(total_rpm_bed).replace(".bedGraph", ".bw"),
        ]
    )
    print(f"[meapy] created total track: {total_bed}")


def merge_two_strand_methylation(cpg_report: Path, output_site_report: Path) -> None:
    pending_plus: Optional[Tuple[str, int, int, int, str, str]] = None
    with cpg_report.open("r") as in_handle, output_site_report.open("w") as out_handle:
        for raw_line in in_handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            strand = parts[2]
            methyl = int(parts[3])
            unmethyl = int(parts[4])
            context = parts[5]
            tri = parts[6]

            if strand == "+":
                pending_plus = (chrom, pos, methyl, unmethyl, context, tri)
                continue

            if strand == "-" and pending_plus is not None:
                plus_chrom, plus_pos, plus_methyl, plus_unmethyl, plus_context, plus_tri = pending_plus
                if plus_chrom == chrom and pos == plus_pos + 1:
                    total_methyl = plus_methyl + methyl
                    total_unmethyl = plus_unmethyl + unmethyl
                    total_depth = total_methyl + total_unmethyl
                    if total_depth > 0:
                        methyl_pct = (total_methyl / total_depth) * 100.0
                        out_handle.write(
                            f"{chrom}\t{plus_pos}\t{pos}\t{methyl_pct:.6f}\t"
                            f"{total_methyl}\t{total_unmethyl}\t{context}\t{plus_tri}\n"
                        )
                    else:
                        out_handle.write(
                            f"{chrom}\t{plus_pos}\t{pos}\tNaN\t"
                            f"{total_methyl}\t{total_unmethyl}\t{context}\t{plus_tri}\n"
                        )
                pending_plus = None


def methylation_site_report_to_bedgraph(
    site_report: Path, output_bedgraph: Path, min_depth: int
) -> None:
    with site_report.open("r") as in_handle, output_bedgraph.open("w") as out_handle:
        for raw_line in in_handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1]) - 1
            end = int(parts[2])
            value = parts[3]
            methyl = int(parts[4])
            unmethyl = int(parts[5])
            if methyl + unmethyl >= max(1, min_depth):
                out_handle.write(f"{chrom}\t{start}\t{end}\t{value}\n")


def methyl_to_loc_bedgraph(input_bedgraph: Path, output_bedgraph: Path) -> None:
    with input_bedgraph.open("r") as in_handle, output_bedgraph.open("w") as out_handle:
        for raw_line in in_handle:
            line = raw_line.strip()
            if not line or line.startswith("track"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            out_handle.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t1\n")


def generate_wgbs_tracks_python(
    bam_prefix: str,
    strain1: str,
    strain2: str,
    refmap1: str,
    refmap2: str,
    tracks_output_dir: str,
    chrom_sizes: str,
    min_mapq: int,
    filter_flag: int,
    min_depth: int,
) -> None:
    output_dir = Path(tracks_output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_name = Path(bam_prefix).name

    make_tracks_from_bam(
        bam_prefix, strain1, chrom_sizes, tracks_output_dir, min_mapq, filter_flag, False, None
    )
    make_tracks_from_bam(
        bam_prefix, strain2, chrom_sizes, tracks_output_dir, min_mapq, filter_flag, False, None
    )
    make_tracks_from_bam(
        bam_prefix, "total", chrom_sizes, tracks_output_dir, min_mapq, filter_flag, False, None
    )

    pre1 = Path(f"{bam_prefix}_{strain1}_preProject.CpG_report.txt").expanduser().resolve()
    pre2 = Path(f"{bam_prefix}_{strain2}_preProject.CpG_report.txt").expanduser().resolve()
    pre_total = Path(f"{bam_prefix}_total.CpG_report.txt").expanduser().resolve()
    if not pre1.is_file() or not pre2.is_file() or not pre_total.is_file():
        raise MeapyError(
            "Expected WGBS CpG report files from alignReads were not found. "
            "Run with Bismark mode and verify alignReads outputs."
        )

    site1 = output_dir / f"{run_name}_{strain1}_preProject.CpG_site_report.txt"
    site2 = output_dir / f"{run_name}_{strain2}_preProject.CpG_site_report.txt"
    site_total = output_dir / f"{run_name}_total.CpG_site_report.txt"
    merge_two_strand_methylation(pre1, site1)
    merge_two_strand_methylation(pre2, site2)
    merge_two_strand_methylation(pre_total, site_total)

    pre_bed1 = output_dir / f"{run_name}_{strain1}_preProject_methyl.bedGraph"
    pre_bed2 = output_dir / f"{run_name}_{strain2}_preProject_methyl.bedGraph"
    methylation_site_report_to_bedgraph(site1, pre_bed1, min_depth)
    methylation_site_report_to_bedgraph(site2, pre_bed2, min_depth)
    sort_bedgraph_file(pre_bed1)
    sort_bedgraph_file(pre_bed2)

    methyl1 = output_dir / f"{run_name}_{strain1}_methyl.bedGraph"
    methyl2 = output_dir / f"{run_name}_{strain2}_methyl.bedGraph"
    project_bedgraph(pre_bed1, Path(refmap1).expanduser().resolve(), methyl1)
    project_bedgraph(pre_bed2, Path(refmap2).expanduser().resolve(), methyl2)
    sort_bedgraph_file(methyl1)
    sort_bedgraph_file(methyl2)
    ensure_nonempty_bedgraph(methyl1, chrom_sizes)
    ensure_nonempty_bedgraph(methyl2, chrom_sizes)

    methyl_total = output_dir / f"{run_name}_total_methyl.bedGraph"
    methylation_site_report_to_bedgraph(site_total, methyl_total, min_depth)
    sort_bedgraph_file(methyl_total)
    ensure_nonempty_bedgraph(methyl_total, chrom_sizes)

    methyl1_loc = output_dir / f"{run_name}_{strain1}_methyl_LOC.bedGraph"
    methyl2_loc = output_dir / f"{run_name}_{strain2}_methyl_LOC.bedGraph"
    methyl_total_loc = output_dir / f"{run_name}_total_methyl_LOC.bedGraph"
    methyl_to_loc_bedgraph(methyl1, methyl1_loc)
    methyl_to_loc_bedgraph(methyl2, methyl2_loc)
    methyl_to_loc_bedgraph(methyl_total, methyl_total_loc)
    ensure_nonempty_bedgraph(methyl1_loc, chrom_sizes)
    ensure_nonempty_bedgraph(methyl2_loc, chrom_sizes)
    ensure_nonempty_bedgraph(methyl_total_loc, chrom_sizes)

    for bed in [methyl1, methyl2, methyl_total, methyl1_loc, methyl2_loc, methyl_total_loc]:
        run_command(
            [
                "bedGraphToBigWig",
                str(bed),
                str(Path(chrom_sizes).expanduser().resolve()),
                str(bed).replace(".bedGraph", ".bw"),
            ]
        )


def main() -> int:
    args = build_parser().parse_args()
    try:
        require_path(args.reads1, kind="file")
        if args.threads < 1:
            raise MeapyError("--threads must be >= 1.")
        require_path(args.genome_input, kind="file")
        bam_prefix_path = Path(args.bam_prefix).expanduser().resolve()
        tracks_output_dir = (
            args.tracks_output_dir
            if args.tracks_output_dir
            else str(bam_prefix_path.parent)
        )

        reads1 = Path(args.reads1).expanduser().resolve()
        reads2: Optional[Path] = None
        if args.read_layout == "paired":
            if not args.reads2:
                raise MeapyError("--reads2 is required when --read-layout paired.")
            reads2 = Path(require_path(args.reads2, kind="file"))

        refmap1, refmap2 = default_refmaps(args.genome_input, args.strain1, args.strain2)
        refmap1 = args.refmap1 or refmap1
        refmap2 = args.refmap2 or refmap2
        require_path(refmap1, kind="file")
        require_path(refmap2, kind="file")

        selected_aligner = args.aligner
        if selected_aligner == "auto":
            if args.assay == "wgbs":
                selected_aligner = "bismark"
            elif args.assay == "rna":
                selected_aligner = "minimap2" if args.long else "star"
            else:
                selected_aligner = "bowtie2"

        if args.min_mapq is None:
            default_min_mapq_by_aligner = {
                "star": 255,
                "bowtie2": 10,
                "bismark": 1,
                "bwa": 1,
                "tophat2": 1,
                "minimap2": 20 if args.long else 1,
            }
            effective_min_mapq = default_min_mapq_by_aligner.get(selected_aligner, 1)
        else:
            effective_min_mapq = args.min_mapq
        if effective_min_mapq < 0:
            raise MeapyError("--min-mapq must be >= 0.")

        if args.long:
            if args.assay != "rna":
                raise MeapyError("Long-read mode (--long) is currently supported only with --assay rna.")
            if args.read_layout != "single":
                raise MeapyError("Long-read RNA mode currently requires --read-layout single.")
            if selected_aligner != "minimap2":
                raise MeapyError("Long-read RNA mode requires aligner minimap2.")

        reference_genome_arg = args.reference_genome
        if args.quick_start and not reference_genome_arg:
            reference_genome_arg = infer_reference_genome(
                args.genome_input, args.reference_fasta
            )
        if not reference_genome_arg:
            raise MeapyError(
                "--reference-genome is required (or use --quick-start with --reference-fasta)."
            )
        reference_fasta = Path(require_path(reference_genome_arg, kind="file"))
        pseudogenome_fasta = Path(require_path(args.genome_input, kind="file"))

        if args.assay == "wgbs":
            if selected_aligner != "bismark":
                raise MeapyError("WGBS alignment requires aligner bismark.")
            bam_prefix_path = Path(args.bam_prefix).expanduser().resolve()
            out_dir = bam_prefix_path.parent
            run_name = bam_prefix_path.name
            concat_folder = Path(args.genome_input).expanduser().resolve().parent
            ref_folder = reference_fasta.parent

            run_bismark_alignment(
                genome_folder=concat_folder,
                reads1=reads1,
                reads2=reads2,
                output_bam=Path(f"{args.bam_prefix}_{args.strain1}_{args.strain2}.bam").expanduser().resolve(),
                output_name=f"{run_name}_{args.strain1}_{args.strain2}",
                threads=args.threads,
            )
            run_bismark_alignment(
                genome_folder=ref_folder,
                reads1=reads1,
                reads2=reads2,
                output_bam=Path(f"{args.bam_prefix}_total.bam").expanduser().resolve(),
                output_name=f"{run_name}_total",
                threads=args.threads,
            )
            combined_bam = Path(f"{args.bam_prefix}_{args.strain1}_{args.strain2}.bam").expanduser().resolve()
            split_allelic_bams_from_concat(
                concat_bam=combined_bam,
                strain1=args.strain1,
                strain2=args.strain2,
                reference_fasta=reference_fasta,
                output_bam1=Path(f"{args.bam_prefix}_{args.strain1}.bam").expanduser().resolve(),
                output_bam2=Path(f"{args.bam_prefix}_{args.strain2}.bam").expanduser().resolve(),
                min_mapq=effective_min_mapq,
                exact_mapq=None,
            )
            combined_cpg = run_bismark_methyl_extractor(
                input_bam=combined_bam,
                genome_folder=concat_folder,
                output_dir=out_dir,
                output_prefix=f"{run_name}_{args.strain1}_{args.strain2}",
                is_paired=reads2 is not None,
                threads=args.threads,
            )
            split_cpg_by_strain(
                combined_cpg_report=combined_cpg,
                strain1=args.strain1,
                strain2=args.strain2,
                output1=out_dir / f"{run_name}_{args.strain1}_preProject.CpG_report.txt",
                output2=out_dir / f"{run_name}_{args.strain2}_preProject.CpG_report.txt",
            )
            total_bam = Path(f"{args.bam_prefix}_total.bam").expanduser().resolve()
            run_bismark_methyl_extractor(
                input_bam=total_bam,
                genome_folder=ref_folder,
                output_dir=out_dir,
                output_prefix=f"{run_name}_total",
                is_paired=reads2 is not None,
                threads=args.threads,
            )
        else:
            run_python_alignment(
                reads1=reads1,
                reads2=reads2,
                strain1_name=args.strain1,
                strain2_name=args.strain2,
                pseudogenome_fasta=pseudogenome_fasta,
                reference_fasta=reference_fasta,
                bam_prefix=args.bam_prefix,
                aligner=selected_aligner,
                long_mode=args.long,
                threads=args.threads,
                split_min_mapq=effective_min_mapq,
            )

        chrom_sizes = ensure_chrom_sizes(
            args.chrom_sizes, args.quick_start, reference_fasta, args.genome_input
        )
        require_path(chrom_sizes, kind="file")
        if args.assay == "wgbs":
            generate_wgbs_tracks_python(
                bam_prefix=args.bam_prefix,
                strain1=args.strain1,
                strain2=args.strain2,
                refmap1=refmap1,
                refmap2=refmap2,
                tracks_output_dir=tracks_output_dir,
                chrom_sizes=chrom_sizes,
                min_mapq=effective_min_mapq,
                filter_flag=args.filter_flag,
                min_depth=args.min_depth,
            )
        else:
            generate_tracks_python(
                bam_prefix=args.bam_prefix,
                strain1=args.strain1,
                strain2=args.strain2,
                refmap1=refmap1,
                refmap2=refmap2,
                tracks_output_dir=tracks_output_dir,
                chrom_sizes=chrom_sizes,
                min_mapq=effective_min_mapq,
                filter_flag=args.filter_flag,
                assay=args.assay,
                read_layout=args.read_layout,
                aligner=selected_aligner,
                se_extension=args.se_extension,
            )
    except MeapyError as exc:
        print(f"[meapy] error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
