# MEA: a Methylomic and Epigenomic Allele-specific analysis pipeline

![MEA Pipeline](doc/MEA-diagram.png)


Building on our previous work on the ALEA software package, a computational toolbox for allele-specific epigenomics analysis incorporating allelic variation data, we detail here ALEA’s successor - MEA. This package dramatically increases the functionality and usability of ALEA, incorporating allelic variation with existing resources, allowing for the identification of significant associations of epigenetic modifications and specific allelic variants in human and mouse cells. Similar to ALEA, MEA provides a customizable pipeline for allele-specific analysis of next-generation sequencing data which takes raw sequencing data for ChIP-seq, RNA-seq and DNA Methylation analysis, producing a UCSC track hub. MEA takes advantage of the available genomic resources for human (The 1000 Genomes Project Consortium) and mouse (The Mouse Genome Project) to reconstruct diploid in silico genomes for human samples or hybrid mouse samples. Then, for each accompanying ChIP-seq, RNA-seq or DNA Methylation dataset, MEA generates two wig files from short reads aligned differentially to each haplotype. This pipeline has been validated using human and hybrid mouse ChIPseq, RNAseq and DNA Methylation data (See [Test Data](#test-data)). 

## Please refer to the MEA user guide (PDF) for installation and run instructions
  
* [Docker installation](https://github.com/julienrichardalbert/MEA/raw/master/docker)  
* [distribution package](https://github.com/julienrichardalbert/MEA/raw/master/dist/mea.1.0.tar.gz)
* [documentation (PDF)](https://github.com/julienrichardalbert/MEA/raw/master/doc/MEA_UserGuide_v1.0.pdf)
* test data: ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/

## MEApy (in progress on `MEApy` branch)

The branch `MEApy` includes a flat Python interface at repository root:

- `meapy.py`
- `meapy_phase_vcf.py`
- `meapy_create_genome.py`
- `meapy_align_and_track.py`
- `meapy_project.py`

### Install dependencies (one go)

Use Conda/Mamba for MEApy toolchain installation.

```bash
mamba env create -f environment.yml
conda activate meapy
```

If `mamba` is not installed:

```bash
conda env create -f environment.yml
conda activate meapy
```

Current intent:

- simplify pipeline usage with a native Python interface
- combine alignment + track creation into one command (`align`)
- support a Python projection command (`project`) for refmap + wig/bedGraph/bed inputs
- support Python track generation for ChIP/RNA/WGBS
- defer `create_report` migration until core pipeline workflow is stable

Current dependency status:

- `phase-vcf` is native Python and does not use `mea.config`
- `project` is native Python and does not use `mea.config`
- `create-genome` is native Python and does not use `mea.config` (requires `bcftools` and `samtools`)
- `align` supports Python alignment for ChIP/RNA (`bwa`) and WGBS (`bismark`)
- aligner defaults by assay: ChIP `bowtie2`, RNA `STAR`, WGBS `bismark`; `bwa` remains an optional override and RNA can also use `tophat2`

### TopHat2 legacy environment (optional)

TopHat2 is Python-2-era software. On Python 3-only environments it can fail at runtime.
If you want `--aligner tophat2`, create a separate legacy env and prepend its `bin`
directory to `PATH` while running `meapy`:

```bash
mamba create -n tophat2-legacy -c conda-forge -c bioconda python=2.7 tophat bowtie2 samtools -y
LEGACY_BIN="$(conda run -n tophat2-legacy python -c 'import sys, os; print(os.path.dirname(sys.executable))')"
PATH="$LEGACY_BIN:$PATH" python meapy.py align ... --assay rna --aligner tophat2
```

This keeps the main `meapy` environment modern while enabling optional TopHat2 runs.

### Prebuilt index requirement

`meapy.py align` does not build aligner indexes automatically. If an index is missing, it fails fast and prints the exact command to build it.

One-time index build examples:

```bash
# BWA
bwa index /path/to/genome.fa

# Bowtie2
bowtie2-build /path/to/genome.fa /path/to/genome.bowtie2_index

# STAR
STAR --runMode genomeGenerate --genomeDir /path/to/star_index --genomeFastaFiles /path/to/genome.fa --runThreadN 4

# Bismark
bismark_genome_preparation --bowtie2 /path/to/genome_folder
```

### Daily paper/news/weather digest script

`daily_research_digest.py` is a standalone Python script that compiles:

- recent papers by hard-coded keywords (OpenAlex)
- recent papers by followed scientists (OpenAlex author IDs, optional name-to-ID map)
- optional recent papers by keywords and followed authors (PubMed)
- Google Scholar alerts ingestion from your inbox (IMAP; preferred for quota)
- top stories for Paris, Canada, France, plus one Onion story (RSS)
- daily weather (Open-Meteo)

Setup:

```bash
cp daily_digest_config.example.json daily_digest_config.json
```

Edit `daily_digest_config.json` with your keywords, followed authors, location, and email settings.
For security, use environment-variable fields (`*_env`) for credentials instead of raw passwords in JSON.
For OpenAlex keyword search, tune `research.openalex_per_keyword`; set `research.openalex_biorxiv_only=true` to restrict keyword results to bioRxiv.
For OpenAlex followed authors, `research.followed_openalex_authors` can be either a list of IDs or a dict of `Author Name: OpenAlex ID`.
For PubMed, set `research.pubmed_enabled=true`; PubMed author tracking defaults to names from the OpenAlex author dict keys.
Optionally set `research.followed_pubmed_authors` to override with a custom PubMed author list.
For Google Scholar alerts mode (recommended), set:
- `google_scholar.imap_alerts.username_env` and `google_scholar.imap_alerts.password_env`
- optional `label` to a Gmail label containing your scholar alert emails

For Brevo, use: host `smtp-relay.brevo.com`, port `587`, `smtp.username_env` and `smtp.password_env`.

Set environment variables before running:

```bash
export SCHOLAR_IMAP_USERNAME="your_alert_mailbox@example.com"
export SCHOLAR_IMAP_PASSWORD="your_imap_app_password"
export DIGEST_SMTP_USERNAME="your_brevo_smtp_login"
export DIGEST_SMTP_PASSWORD="your_brevo_smtp_key"
```

Dry-run (prints report only):

```bash
python3 daily_research_digest.py --config daily_digest_config.json --dry-run
```

Dry-run and save digest to file:

```bash
python3 daily_research_digest.py --config daily_digest_config.json --dry-run --save-digest ./daily_digest_latest.txt
```

Send email:

```bash
python3 daily_research_digest.py --config daily_digest_config.json
```

Cron at 07:00 every day:

```bash
0 7 * * * cd /Users/jra/github/MEA && /usr/bin/python3 daily_research_digest.py --config /Users/jra/github/MEA/daily_digest_config.json
```

Launchd (recommended on macOS) at 07:00 every day:

1) Ensure runner is executable:

```bash
chmod +x /Users/jra/github/MEA/run_daily_digest.sh
```

2) Copy plist into your user LaunchAgents folder:

```bash
mkdir -p ~/Library/LaunchAgents
cp /Users/jra/github/MEA/com.jra.daily-research-digest.plist ~/Library/LaunchAgents/
```

3) Load and start schedule:

```bash
launchctl unload ~/Library/LaunchAgents/com.jra.daily-research-digest.plist 2>/dev/null || true
launchctl load ~/Library/LaunchAgents/com.jra.daily-research-digest.plist
```

4) Trigger once immediately (optional):

```bash
launchctl start com.jra.daily-research-digest
```

5) Check logs:

```bash
tail -f /Users/jra/github/MEA/daily_digest.log /Users/jra/github/MEA/daily_digest.err.log
```

Examples:

- `python3 meapy.py phase-vcf --haps-dir <dir> --input-vcf <vcf> --output-prefix <prefix> --compress`
- `python3 meapy.py create-genome --reference-fasta <ref.fa> --phased-vcf <phased.vcf.gz> --strain1 hap1 --strain2 hap2 --output-dir <dir>`
- `python3 meapy.py align --read-layout single --reads1 <reads.fastq.gz> --genome-input <genome.fa> --reference-genome <ref.fa> --strain1 hap1 --strain2 hap2 --bam-prefix <out/prefix> --assay chip --aligner bowtie2 --chrom-sizes <build.chrom.sizes>`
- `python3 meapy.py align --read-layout paired --reads1 <R1.fastq.gz> --reads2 <R2.fastq.gz> --genome-input <genome.fa> --reference-genome <ref.fa> --strain1 hap1 --strain2 hap2 --bam-prefix <out/prefix> --assay wgbs --aligner bismark --chrom-sizes <build.chrom.sizes>`
- `python3 meapy.py align --quick-start --reference-fasta <ref.fa> --read-layout single --reads1 <reads.fastq.gz> --genome-input <create_genome_dir/hap1_hap2.fasta> --strain1 hap1 --strain2 hap2 --bam-prefix <out/prefix> --assay chip`
- `python3 meapy.py project --input <track.bedGraph> --input-format bedgraph --input-refmap <strain.fasta.refmap> --output-bedgraph <projected.bedGraph>`
- `python3 meapy.py doctor --assay all`
- `python3 meapy.py validate --assay chip --bam-prefix <out/prefix> --strain1 hap1 --strain2 hap2 --tracks-output-dir <out>`



## Credits
Hamid Younesy, Torsten Möller, Alireza Heravi-Moussavi, Jeffrey B. Cheng, Joseph F. Costello, Matthew C. Lorincz, Mohammad M. Karimi and Steven J. M. Jones, "ALEA: a toolbox for allele-specific epigenomics analysis." Bioinformatics 30.8 (2014): 1172-1174. [[link to paper](http://bioinformatics.oxfordjournals.org/content/30/8/1172.long)]

Julien Richard Albert, Tasuku Koike, Hamid Younesy, Richard Thompson, Aaron B. Bogutz, Mohammad M. Karimi and Matthew C. Lorincz, "Development and application of an integrated allele-specific pipeline for methylomic and epigenomic analysis (MEA)." BMC Genomics201819:463 (https://doi.org/10.1186/s12864-018-4835-2)




