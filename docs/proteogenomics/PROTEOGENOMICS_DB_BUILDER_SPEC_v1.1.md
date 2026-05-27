# Feature Spec: Proteogenomics Database Builder

> **Version**: 1.1 — May 2026 (post-validation; codebase-reconciled May 21, 2026)
> **Author**: Brett Phinney / UC Davis Proteomics Core
> **Priority**: P2 (after Proteogenomics Glossary addendum)
> **Prereqs**: `DIANN_HPC_INTEGRATION_SPEC.md`, `CLAUDE_EXPORT_PROTEOGENOMICS_ADDENDUM.md`
> **Pairs with**: SLIMS data delivery (DNA Technologies Core, UC Davis Genome Center)
> **Status**: Pipeline architecture validated end-to-end on Hive (May 20, 2026). Final
> FASTA produced: 67,386 entries (66,046 REF + 1,340 NOVEL_GENE) on mouse test dataset.
> All 8 pipeline stages verified. 11 spec lessons captured. File names and packaging
> conventions reconciled against actual DE-LIMP code on May 21, 2026 — see §12
> (Files to Create / Modify) for the corrected paths.

---

## Changes from v1.0

This revision incorporates findings from the May 20, 2026 end-to-end validation run on
Hive HPC. The pipeline architecture is unchanged in shape but several specific choices
have been corrected based on what actually worked:

1. **STAR replaces HISAT2** as the default aligner — STAR indices are pre-staged by
   bioinfocore-grp, STAR is the proteogenomics-community standard (Galaxy-P, Jagtap),
   and Hive's 2 TB RAM nodes make HISAT2's lower-memory argument irrelevant.
2. **Reference data paths corrected** to actual Hive locations (`/quobyte/bioinfocore-grp/genomes/`
   and `/quobyte/proteomics-grp/de-limp/references/`). The original spec's `/share/proteomics/`
   paths don't exist on Hive.
3. **Apptainer container deferred / unnecessary** — Hive's central modules cover the
   entire pipeline (fastp, bowtie2, STAR, samtools, stringtie via conda, gffread, transdecoder,
   diamond, seqkit). One small conda env (`proteog_helpers`) covers the Python helper
   needs (gffutils, biopython). No custom container is required.
4. **rRNA pre-filter step added as mandatory** — discovered during validation when
   the test dataset showed 73% multi-mapping. rRNA filtering is essential for any library
   that isn't strictly polyA-enriched.
5. **Adaptive STAR thresholds based on read length** — STAR's defaults are calibrated for
   150bp+ reads (which is what DNA Tech Core actually delivers). Validation tested with
   92bp reads which need relaxed thresholds. The pipeline must detect and adapt.
6. **`gffcompare` step added** for NOVEL_ISOFORM detection — without it, `stringtie --merge`
   collapses novel-isoform discoveries into REF entries. This is the difference between
   discovering "1,340 novel genes" vs "1,340 novel genes + N hundred novel isoforms of
   known genes."
7. **Header format updated** to the validated `sp|<protein_id>|<symbol>_<TAG> source=... ORF_type=... ...`
   structure with 7 key=value metadata fields. Original Jagtap-style `_u_/_c_/_i_/_o_`
   suffix codes are replaced by explicit `source=REF/NOVEL_GENE/NOVEL_ISOFORM/UNPARSED` tags.
8. **Quality gates added** — pipeline halts on (a) species mismatch detected pre-alignment,
   (b) uniquely-mapped rate below read-length-tiered threshold, (c) parse-malformed
   headers in rewriter output.
9. **SRA / ENA accession verification step added** — for the case where users provide
   accessions instead of SLIMS URLs. Validation surfaced two real bugs caused by skipping
   this check (wrong species, wrong library type).
10. **StringTie module choice** — Hive's native `stringtie/2.2.1` module has a fatal bgzf
    assertion crash; use `conda/stringtie/3.0.3` instead. Bug filed to HPC@UCD.
11. **Conda channel order** must be `-c conda-forge -c bioconda --strict-channel-priority`,
    not the reverse. Reverse order causes the solver to pick python-2-only versions of
    biopython.

---

## Table of Contents

1. [Overview & Value Proposition](#1-overview--value-proposition)
2. [Background: DNA Tech Core Data Delivery](#2-background-dna-tech-core)
3. [Architecture Overview](#3-architecture)
4. [Phase 1: SLIMS Ingestion Tab](#4-phase-1-slims-ingestion)
5. [Phase 2: RNA-seq Pipeline](#5-phase-2-rnaseq-pipeline)
6. [Phase 3: Variant Caller (Optional)](#6-phase-3-variant-caller)
7. [Phase 4: FASTA Assembly & Registry](#7-phase-4-fasta-assembly)
8. [Phase 5: Integration with DIA-NN New Search Tab](#8-phase-5-diann-integration)
9. [Quality Gates](#9-quality-gates)
10. [UI Layout](#10-ui-layout)
11. [Reactive Values & Module Placement](#11-reactive-values)
12. [Files to Create / Modify](#12-files)
13. [Dependencies](#13-dependencies)
14. [Cost & Time Estimates](#14-cost-and-time)
15. [Testing Checklist](#15-testing-checklist)
16. [Validation Lessons Captured](#16-validation-lessons)

---

## 1. Overview & Value Proposition

DE-LIMP can already search MS data against a proteogenomics-expanded FASTA when one is
provided. This spec adds the upstream half of that workflow: **building the FASTA from
sample-matched RNA-seq data delivered by the UC Davis DNA Technologies Core.**

The user journey:

```
1. Submit RNA-seq samples to DNA Tech Core (same biological samples used for proteomics)
2. Receive SLIMS notification email with download URL
3. In DE-LIMP: paste SLIMS URL into the new "Build Database" tab
4. Click "Build Proteogenomics FASTA"
5. DE-LIMP downloads, runs QC, aligns, assembles transcripts, predicts ORFs, and merges
   with the reference proteome
6. The resulting FASTA appears in the New Search tab's pre-staged-databases dropdown,
   tagged with the project name
7. Run DIA-NN as usual — results flow through the existing proteogenomics-aware pipeline
   (proteogenomics_glossary.txt ships in the Claude export, etc.)
```

A bench scientist who has never run a single bioinformatics command line tool can now
produce a proteogenomics database from their own RNA-seq data, indistinguishable in
quality from what the Jagtap lab or Galaxy-P produce manually.

---

## 2. Background: DNA Tech Core

### Data delivery format

The DNA Tech Core delivers Illumina sequencing data via SLIMS at
`slimsdata.genomecenter.ucdavis.edu`. Each submission produces a URL with a random
ID that acts as both the location and the access token:

```
http://slimsdata.genomecenter.ucdavis.edu/Data/{RANDOM_ID}/Unaligned/
├── Project_{Name}/
│   ├── Sample_001_R1.fastq.gz
│   ├── Sample_001_R2.fastq.gz
│   ├── ...
├── Reports/
│   └── demultiplex_stats.html
└── checksums.md5
```

Key facts:
- Files are gzip-compressed FASTQ, demultiplexed by sample, dual-indexed
- Each sample has R1 + R2 (paired-end is the default for RNA-seq)
- md5 checksums are provided; we verify automatically
- Data is free for 1 month then may be deleted → DE-LIMP downloads on the user's behalf
- The URL is the auth token; user just pastes it, no separate login

Recommended download: `wget -r -nH -nc -R 'index.html*'`. SLIMS is reachable from Hive
login nodes (verified during validation). Compute nodes may not have outbound HTTP,
so the wget step runs on a login node, not as a SLURM job.

### Library types

DNA Tech Core offers:
- **Standard mRNA-Seq** (poly-A enriched) — most common for differential expression
- **Total RNA-seq with rRNA depletion** — captures non-coding RNAs, often preferred
  for proteogenomics because it captures more lincRNA and small ORF transcripts
- **Stranded RNA-seq** — directional; tells us which strand was transcribed
- **3'-Tag-Seq** — short 3' tags; **not suitable for proteogenomics** (doesn't cover
  full transcripts). DE-LIMP must detect and refuse this.
- **miRNA-seq** — **not suitable**

### Read length expectations

DNA Tech Core's default RNA-seq output is **PE150** (150bp paired-end). This is the
target for STAR's default thresholds. The pipeline must adapt to shorter reads when
encountered (e.g., re-analyzing older PE100 or PE92 data from ENCODE/SRA).

---

## 3. Architecture

### Overall data flow

```
DNA Tech Core SLIMS URL  (or SRA/ENA accession)
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  Phase 1: Ingestion (DE-LIMP on Hive login)                      │
│  - Paste URL → wget mirror to /quobyte/proteomics-grp/.../rnaseq/│
│  - Verify md5 checksums                                          │
│  - Detect sample structure (count R1/R2 pairs)                   │
│  - For SRA accessions: query ENA metadata, verify species/lib    │
│  - Library-type questionnaire (refuse 3'-Tag-Seq, miRNA-Seq)     │
└─────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  Phase 2: RNA-seq Pipeline (SLURM on Hive)                       │
│                                                                   │
│  For each sample (parallel sbatch array):                         │
│    fastp        → trim + filter, detect read length              │
│    bowtie2      → filter vs organism-specific rRNA FASTA         │
│                   (keep unmapped reads only)                     │
│    STAR         → align with read-length-adaptive thresholds     │
│    [QC gate]    → halt if uniquely-mapped rate < tier threshold  │
│                                                                   │
│  Then once:                                                       │
│    stringtie --merge → unified transcript GTF                    │
│    gffcompare        → annotate against reference for NOVEL_ISOFORM │
│    gffread           → extract transcript FASTA from genome      │
│    TransDecoder      → predict ORFs (--single_best_only)         │
│                                                                   │
│  Finally:                                                         │
│    rewrite_headers.py → produce sp|ID|SYM_TAG ... format         │
│                                                                   │
│  Output: predicted_orfs.fasta + merged.gtf (kept for coord lookup)│
└─────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  Phase 3: Optional Variant Encoding                              │
│  If user provides matched DNA-seq VCF → variant proteoforms      │
│  (Off by default; skipped if no VCF)                             │
└─────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  Phase 4: FASTA Assembly                                         │
│  Concatenate: predicted_orfs.fasta + UniProt reference +         │
│               variants (if any) + cRAP contaminants              │
│  Deduplicate by sequence (seqkit rmdup -s)                       │
│  Register in registry.json with project metadata                 │
│  Path: /quobyte/proteomics-grp/de-limp/databases/proteogenomics/ │
└─────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  Phase 5: DIA-NN New Search integration                          │
│  - FASTA appears in pre-staged-databases dropdown with 🧬 tag    │
│  - Selecting it triggers auto-warning + parameter recommendations│
│  - Proteogenomics-aware result flow activates downstream         │
└─────────────────────────────────────────────────────────────────┘
```

### SLURM dependency chain

```
job_fastp_array     (N parallel, one per sample, ~10 min)
  ↓
job_rrna_array      (N parallel, one per sample, ~5 min)
  ↓
job_star_array      (N parallel, one per sample, ~15-30 min depending on depth)
  ↓
job_merge           (1 job, depends on all star, ~10 min)
  ↓
job_gffcompare      (1 job, ~2 min)
  ↓
job_orf_predict     (1 job, ~30 min)
  ↓
job_assemble_fasta  (1 job, ~5 min)
```

End-to-end for 12 samples on `high` partition with `--account=genome-center-grp`:
typically 3–6 hours wall time. The user can close the browser; status polls via
`squeue` and persists in SQLite (Core Facility mode) or in-memory.

### Account flags

All sbatch scripts include:
```bash
#SBATCH --account=genome-center-grp   # Brett's account for high partition
#SBATCH --partition=high              # 30-day walltime, fast start
```

For users without `genome-center-grp` access, fall back to:
```bash
#SBATCH --account=publicgrp
#SBATCH --partition=low               # 7-day walltime, preemptible
#SBATCH --qos=publicgrp-low-qos
```

DE-LIMP detects available accounts at startup via `sacctmgr show user $USER` and
picks the appropriate one. Surface to user as "Submitting with `genome-center-grp`
access (priority)" or "Submitting with `publicgrp` access (low priority — jobs may
be preempted)."

---

## 4. Phase 1: Ingestion

### Two input modes

The Build Database tab accepts two ingestion sources:

**Mode A — SLIMS URL** (primary use case, UC Davis users)
```
SLIMS URL:  [http://slimsdata.genomecenter.ucdavis.edu/Data/i5om268pkp/Unaligned/]
```

**Mode B — SRA/ENA accession** (re-analysis, external collaborators)
```
Accession:  [SRR1303776, SRR1303777]   (comma-separated, max 24)
```

### SLIMS URL workflow

```r
scan_slims_url <- function(slims_url) {
  if (!grepl("^https?://slimsdata\\.genomecenter\\.ucdavis\\.edu/Data/[a-z0-9]+/",
             slims_url)) {
    return(list(success = FALSE,
                error = "URL does not match SLIMS format."))
  }

  index_html <- tryCatch(
    httr::content(httr::GET(slims_url), as = "text"),
    error = function(e) NULL
  )
  if (is.null(index_html)) {
    return(list(success = FALSE,
                error = paste0("Could not reach SLIMS from Hive login node. ",
                               "If this persists, check that ",
                               "slimsdata.genomecenter.ucdavis.edu is reachable.")))
  }

  file_pattern <- "[A-Za-z0-9_.-]+\\.fastq\\.gz"
  files <- unique(regmatches(index_html,
    gregexpr(file_pattern, index_html))[[1]])

  if (length(files) == 0) {
    return(list(success = FALSE,
                error = "No .fastq.gz files found. Data may have been deleted (1-month retention)."))
  }

  r1_files <- grep("_R1[._]", files, value = TRUE)
  r2_files <- grep("_R2[._]", files, value = TRUE)
  is_paired <- length(r1_files) > 0 && length(r2_files) > 0

  list(
    success      = TRUE,
    url          = slims_url,
    n_samples    = if (is_paired) length(r1_files) else length(files),
    sample_names = if (is_paired) gsub("_R1[._].*", "", r1_files) else gsub("\\.fastq\\.gz$", "", files),
    files        = files,
    is_paired    = is_paired,
    has_md5      = grepl("checksums\\.md5", index_html)
  )
}
```

Download runs on the Hive **login node**, not as a SLURM job, since compute nodes
may not have outbound HTTP. From DE-LIMP's R session (which is itself running on
Hive via the de-limp.sif container):

```r
launch_slims_download <- function(slims_url, project_name) {
  dest_dir <- file.path("/quobyte/proteomics-grp/de-limp/rnaseq",
                        project_name)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  # Background process so it doesn't block the Shiny session
  cmd <- sprintf(
    "cd %s && wget -r -nH --cut-dirs=3 -nc -R 'index.html*' '%s' > wget.log 2>&1 && md5sum -c checksums.md5 > md5_verify.log 2>&1",
    dest_dir, slims_url
  )

  system(sprintf("nohup bash -c %s > %s/download.log 2>&1 &",
                 shQuote(cmd), dest_dir))
}
```

Poll `download.log` for completion. Surface progress to UI.

### SRA/ENA accession workflow

**Critical step learned in validation:** verify ENA metadata BEFORE downloading.
The validation run was given SRR accessions claimed to be K562 human RNA-seq;
they turned out to be Mus musculus RNA-seq with heavy rRNA contamination. This
wasted ~30 minutes of compute time that a 2-second metadata check would have prevented.

```r
verify_sra_accession <- function(accession) {
  url <- sprintf("https://www.ebi.ac.uk/ena/browser/api/xml/%s", accession)
  resp <- tryCatch(httr::GET(url), error = function(e) NULL)
  if (is.null(resp) || httr::status_code(resp) != 200) {
    return(list(success = FALSE, error = "Could not reach ENA metadata API."))
  }

  xml <- xml2::read_xml(httr::content(resp, as = "text"))

  list(
    success         = TRUE,
    accession       = accession,
    scientific_name = xml2::xml_text(xml2::xml_find_first(xml, ".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME")),
    library_strategy= xml2::xml_text(xml2::xml_find_first(xml, ".//LIBRARY_STRATEGY")),
    library_source  = xml2::xml_text(xml2::xml_find_first(xml, ".//LIBRARY_SOURCE")),
    library_selection = xml2::xml_text(xml2::xml_find_first(xml, ".//LIBRARY_SELECTION")),
    instrument      = xml2::xml_text(xml2::xml_find_first(xml, ".//INSTRUMENT_MODEL")),
    layout          = if (length(xml2::xml_find_all(xml, ".//PAIRED")) > 0) "paired" else "single"
  )
}
```

Display to user before download proceeds:
```
Verifying accessions...
   SRR1303776
     Species:     Mus musculus
     Library:     RNA-Seq (TRANSCRIPTOMIC, cDNA)
     Selection:   PolyA
     Instrument:  Illumina HiSeq 2500
     Layout:      paired

Selected reference: Homo sapiens (GRCh38) — DOES NOT MATCH

[ Change reference to Mus musculus ]  [ Cancel ]
```

**Refuse to proceed if** library_strategy is `Tag-Seq`, `miRNA-Seq`, `OTHER`, or
`small-RNA` — these are unsuitable for proteogenomics novel-ORF discovery.

Download uses ENA URL streaming (not sra-toolkit, which requires `vdb-config --interactive`
that doesn't work in sbatch). ENA mirrors essentially all of SRA.

```bash
# ENA URL pattern (one per read):
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR130/006/SRR1303776/SRR1303776_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR130/006/SRR1303776/SRR1303776_2.fastq.gz
```

Stream-subsample if needed (use `curl URL | zcat | head -n N | gzip > out.fq.gz`,
with `set +o pipefail` to tolerate the SIGPIPE that `head` produces when it closes
the upstream — this was a real bug discovered in validation).

---

## 5. Phase 2: RNA-seq Pipeline

All sbatch scripts use Hive's central modules. No custom container required.

### Stage 2.1 — fastp (adapter trim, quality filter, read-length detection)

```bash
#!/bin/bash -l
#SBATCH --job-name=fastp_${PROJECT}
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --array=1-${N_SAMPLES}
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

module load fastp/0.23.4

SAMPLES=(${SAMPLE_LIST})
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

fastp \
  --in1  rnaseq/${SAMPLE}_R1.fastq.gz \
  --in2  rnaseq/${SAMPLE}_R2.fastq.gz \
  --out1 fastp_out/${SAMPLE}_R1.fastq.gz \
  --out2 fastp_out/${SAMPLE}_R2.fastq.gz \
  --json fastp_out/${SAMPLE}.json \
  --html fastp_out/${SAMPLE}.html \
  --thread 8 \
  --detect_adapter_for_pe \
  --length_required 36
```

After fastp completes, DE-LIMP samples 100 reads from `${SAMPLE}_R1.fastq.gz` and
computes the median length. This drives the read-length tier selection for STAR
(Stage 2.3).

### Stage 2.2 — rRNA pre-filter (MANDATORY)

This step was not in v1.0 of the spec. Validation showed it is essential — even
nominally polyA-enriched libraries can have 5–15% rRNA, and total-RNA libraries
can have 60–90%. Without filtering, StringTie's transcript assembly is overwhelmed
by ribosomal multi-mappers.

```bash
#!/bin/bash -l
#SBATCH --job-name=rrna_filter_${PROJECT}
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --array=1-${N_SAMPLES}
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

module load bowtie2/2.5.2

SAMPLES=(${SAMPLE_LIST})
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

# rRNA index path is per-organism, built once and registered
RRNA_INDEX=${RRNA_INDEX_PREFIX}  # e.g. /quobyte/proteomics-grp/de-limp/references/rrna_index/mm39/rrna

bowtie2 -x ${RRNA_INDEX} \
  -1 fastp_out/${SAMPLE}_R1.fastq.gz \
  -2 fastp_out/${SAMPLE}_R2.fastq.gz \
  --very-sensitive-local \
  --un-conc-gz rrna_filt/${SAMPLE}_norrna_R%.fq.gz \
  -S /dev/null \
  -p 8 \
  --no-unal 2> rrna_filt/${SAMPLE}_rrna_filter.log
```

The rRNA FASTAs are pre-staged by bioinfocore-grp at:
```
/quobyte/bioinfocore-grp/genomes/human/GRCh38.p14/rRNA_human_03-12-2026.fasta
/quobyte/bioinfocore-grp/genomes/mouse/GRCm39/rRNA_mouse_03-12-2026.fasta
```

DE-LIMP builds bowtie2 indices from these once and caches them at:
```
/quobyte/proteomics-grp/de-limp/references/rrna_index/{organism}/
```

The build is tiny (rRNA FASTAs are 1.7 MB; index is ~10 MB) and takes seconds.
Index build is idempotent — if it already exists, skip.

Capture rRNA% from the bowtie2 log and surface in the pipeline summary. This is a
useful quality indicator for users (high rRNA% = library prep issue worth knowing about).

### Stage 2.3 — STAR alignment with adaptive thresholds

**Read-length tier selection** (computed from median post-fastp read length):

```r
select_star_params <- function(read_length) {
  if (read_length >= 130) {
    list(
      tier = "default",
      flags = c(
        "--outFilterMismatchNoverLmax 0.04",  # STAR default
        "--outFilterScoreMinOverLread 0.66",  # STAR default
        "--outFilterMatchNminOverLread 0.66", # STAR default
        "--outFilterMultimapNmax 20"
      ),
      qc_gate_unique_pct = 60
    )
  } else if (read_length >= 100) {
    list(
      tier = "mildly_relaxed",
      flags = c(
        "--outFilterMismatchNoverLmax 0.06",
        "--outFilterScoreMinOverLread 0.50",
        "--outFilterMatchNminOverLread 0.50",
        "--outFilterMultimapNmax 20"
      ),
      qc_gate_unique_pct = 45
    )
  } else if (read_length >= 60) {
    list(
      tier = "significantly_relaxed",
      flags = c(
        "--outFilterMismatchNoverLmax 0.10",
        "--outFilterScoreMinOverLread 0.30",
        "--outFilterMatchNminOverLread 0.30",
        "--outFilterMultimapNmax 20"
      ),
      qc_gate_unique_pct = 25,
      warning = "Read length below 100bp produces noisier transcript assembly. Consider commissioning a higher-quality RNA-seq run."
    )
  } else {
    list(
      tier = "refuse",
      error = "Read length below 60bp is unsuitable for proteogenomics novel-ORF discovery. Please consult the Proteomics Core about appropriate library design."
    )
  }
}
```

Surface tier choice in UI before submission:
```
Median read length: 92bp
Tier: significantly_relaxed
   ⚠ Reads shorter than 100bp produce noisier transcript assembly.
   ⚠ Pipeline will proceed with relaxed thresholds. QC gate: ≥25% uniquely mapped.
[ Proceed ]  [ Cancel ]
```

Sbatch:
```bash
#!/bin/bash -l
#SBATCH --job-name=star_${PROJECT}
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --array=1-${N_SAMPLES}
#SBATCH --time=4:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16

module load star/2.7.11a

SAMPLES=(${SAMPLE_LIST})
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

STAR \
  --runMode alignReads \
  --genomeDir ${STAR_INDEX} \
  --readFilesIn rrna_filt/${SAMPLE}_norrna_R1.fq.gz rrna_filt/${SAMPLE}_norrna_R2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMattributes Standard XS NH \
  --runThreadN 16 \
  --outFileNamePrefix star_out/${SAMPLE}_ \
  ${ADAPTIVE_FLAGS}  # from tier selection
```

**The `--outSAMstrandField intronMotif` flag is critical** — without it StringTie
cannot determine strand for spliced reads and assembly quality drops sharply.

After STAR completes, parse `${SAMPLE}_Log.final.out` for "Uniquely mapped reads %"
and apply the QC gate.

### Stage 2.3a — QC Gate (NEW)

```r
check_alignment_quality <- function(star_log_path, tier_params) {
  log_lines <- readLines(star_log_path)
  unique_pct <- as.numeric(gsub(".*Uniquely mapped reads % \\|\\t(.+)%.*", "\\1",
    grep("Uniquely mapped reads %", log_lines, value = TRUE)))

  if (unique_pct < tier_params$qc_gate_unique_pct) {
    list(
      pass = FALSE,
      unique_pct = unique_pct,
      gate = tier_params$qc_gate_unique_pct,
      tier = tier_params$tier,
      message = sprintf(
        "Alignment quality below threshold. Uniquely mapped: %.1f%%, required: %d%%. Possible causes (in likelihood order): (1) wrong reference genome, (2) heavy contamination beyond rRNA (mitochondrial, host cell line, bacterial), (3) library type unsuited to this pipeline (Ribo-Seq, CLIP-Seq, 3'-Tag-Seq), (4) severe sample degradation. Pipeline halted.",
        unique_pct, tier_params$qc_gate_unique_pct
      )
    )
  } else {
    list(pass = TRUE, unique_pct = unique_pct)
  }
}
```

If gate fails, halt the SLURM dependency chain (cancel pending stringtie/merge/orf jobs)
and surface the diagnostic to the user. Do not produce a partial FASTA.

### Stage 2.4 — Per-sample StringTie

```bash
#!/bin/bash -l
#SBATCH --job-name=stringtie_${PROJECT}
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --array=1-${N_SAMPLES}
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

# IMPORTANT: native stringtie/2.2.1 module has a fatal bgzf assertion crash.
# Use the conda module instead.
module load conda/stringtie/3.0.3

SAMPLES=(${SAMPLE_LIST})
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

stringtie star_out/${SAMPLE}_Aligned.sortedByCoord.out.bam \
  -G ${REFERENCE_GTF} \
  -o stringtie_out/${SAMPLE}.gtf \
  -p 8 \
  ${STRAND_FLAG}  # --rf for reverse-stranded (TruSeq default), --fr for forward, blank for unstranded
```

Strand flag is determined by either:
- User-confirmed choice in the UI (recommended for SLIMS data — DNA Tech Core knows
  the prep type and the user should know which they submitted), or
- `rseqc/5.0.4`'s `infer_experiment.py` against a BED12 derived from the reference GTF
  (auto-detect mode)

Validation default if uncertain: `--rf` (TruSeq stranded, most common).

### Stage 2.5 — Merge per-sample GTFs

```bash
#!/bin/bash -l
#SBATCH --job-name=stringtie_merge_${PROJECT}
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

module load conda/stringtie/3.0.3

ls stringtie_out/*.gtf > stringtie_out/gtf_list.txt
stringtie --merge \
  -G ${REFERENCE_GTF} \
  -o stringtie_out/merged.gtf \
  stringtie_out/gtf_list.txt
```

**Important behavior of `stringtie --merge -G`**: transcripts overlapping reference
exons inherit the reference `gene_id` (so they retain ENSMUST/ENSMUSG IDs and look
like REF). Only fully intergenic loci get fresh `MSTRG` IDs. This means NOVEL_ISOFORM
detection requires the gffcompare step below.

### Stage 2.5a — gffcompare (NEW — for NOVEL_ISOFORM detection)

```bash
module load gffcompare/0.12.6

gffcompare \
  -r ${REFERENCE_GTF} \
  -o stringtie_out/gffcmp \
  stringtie_out/merged.gtf
```

This produces `gffcmp.merged.gtf.tmap` with class codes for every transcript:
- `=` exact match to reference
- `c` contained in reference
- `j` multi-exon novel isoform of a reference gene ← **NOVEL_ISOFORM candidate**
- `e` single-exon novel isoform ← **NOVEL_ISOFORM candidate**
- `i` intronic
- `o` opposite-strand overlap
- `u` unknown (intergenic) ← **NOVEL_GENE candidate** (matches MSTRG.* IDs)

The header rewriter (Stage 2.8) uses this table to assign `source=NOVEL_ISOFORM`
to transcripts with class code `j` or `e`, and `source=NOVEL_GENE` to those with
class code `u`. Everything else becomes `source=REF`.

### Stage 2.6 — Extract transcript FASTA

```bash
module load gffread/0.12.7

gffread -w stringtie_out/merged_transcripts.fa \
  -g ${GENOME_FASTA} \
  stringtie_out/merged.gtf
```

`${GENOME_FASTA}` for mouse is `/quobyte/bioinfocore-grp/genomes/mouse/GRCm39/GRCm39.primary_assembly.genome.fa`.

For human (hg38), bioinfocore-grp does **not** stage a FASTA alongside the GTF —
only the STAR index. DE-LIMP must download it once and cache at:
```
/quobyte/proteomics-grp/de-limp/references/genomes/hg38_GRCh38.p14/genome.fa
```
URL: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz`
(~900 MB compressed, ~3 GB uncompressed, 5–10 min download on Hive)

The download runs from the login node, not as a SLURM job (compute nodes may lack
outbound HTTP).

### Stage 2.7 — TransDecoder ORF prediction

```bash
#!/bin/bash -l
#SBATCH --job-name=transdecoder_${PROJECT}
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

module load transdecoder/5.7.1 diamond/2.1.7

cd transdecoder_out

TransDecoder.LongOrfs \
  -t ../stringtie_out/merged_transcripts.fa \
  -m 100  # minimum ORF length in aa (default 100, good for novel-ORF discovery)

# OPTIONAL homology filter — skip if no DIAMOND DB available
if [ -f "${DIAMOND_UNIPROT_DB}" ]; then
  diamond blastp \
    --query merged_transcripts.fa.transdecoder_dir/longest_orfs.pep \
    --db ${DIAMOND_UNIPROT_DB} \
    --max-target-seqs 1 \
    --outfmt 6 \
    --evalue 1e-5 \
    --threads 16 \
    --out blastp.outfmt6

  TransDecoder.Predict \
    -t ../stringtie_out/merged_transcripts.fa \
    --retain_blastp_hits blastp.outfmt6 \
    --single_best_only
else
  TransDecoder.Predict \
    -t ../stringtie_out/merged_transcripts.fa \
    --single_best_only
fi
```

`--single_best_only` ensures one ORF per transcript (with `.p1` suffix). Without
it, the same transcript can produce multiple ORF candidates and the FASTA bloats.

### Stage 2.8 — Header rewrite (Python helper)

The TransDecoder output FASTA has headers like:
```
>Gene.1::STRG.1017.1::g.1::m.1 type:complete len:124 gc:universal STRG.1017.1:192-563(+)
```

These are not parser-friendly for DE-LIMP. The rewrite step produces:
```
>sp|ENSMUST00000000001.5.p1|Gnai3_MM39TEST source=REF ORF_type=complete strand=+ len=354 coords=ENSMUST00000000001.5:142-1206(+) parent_gene=ENSMUSG00000000001.5 transcript=ENSMUST00000000001.5
>sp|MSTRG.10029.5.p2|MSTRG.10029_MM39TEST source=NOVEL_GENE ORF_type=5prime_partial strand=- len=112 coords=MSTRG.10029.5:344-682(-) parent_gene=MSTRG.10029 transcript=MSTRG.10029.5
```

This format was validated on May 20, 2026 with 67,386 headers produced, **zero malformed**.

The rewriter (`rewrite_transdecoder_headers.py`) lives at:
```
/quobyte/proteomics-grp/de-limp/scripts/rewrite_transdecoder_headers.py
```

It runs in the `proteog_helpers` conda env (biopython + gffutils). Key fields per header:

| Field | Source | Example |
|-------|--------|---------|
| `source` | gffcompare class code | `REF` / `NOVEL_GENE` / `NOVEL_ISOFORM` / `UNPARSED` |
| `ORF_type` | TransDecoder | `complete` / `5prime_partial` / `3prime_partial` / `internal` |
| `strand` | merged.gtf | `+` / `-` |
| `len` | TransDecoder | integer (aa) |
| `coords` | TransDecoder | `transcript_id:start-end(strand)` |
| `parent_gene` | merged.gtf | ENSMUSG... or MSTRG... |
| `transcript` | merged.gtf | ENSMUST... or MSTRG.... |

Project tag (`_MM39TEST` in the example) is glued to the symbol with `_`, set per project.

**Genomic coordinates are NOT in the header** — they're looked up on demand from
the preserved `merged.gtf` using the `parent_gene` field. This is a deliberate design
choice (headers should be immutable; coordinate lookups need to be fresh).

The rewriter must produce **100% parse-clean output**. If any header fails to
parse, it's tagged `source=UNPARSED` and a warning is surfaced. Validation requires
zero UNPARSED entries.

---

## 6. Phase 3: Variant Encoding (Optional, Off by Default)

Unchanged from v1.0. See original spec §6.

---

## 7. Phase 4: FASTA Assembly & Registry

```r
assemble_proteogenomics_fasta <- function(
  project_name,
  uniprot_fasta,
  predicted_orfs_fasta,
  merged_gtf,                # preserved for coord lookup
  variant_fasta = NULL,
  contaminants_fasta = NULL,
  output_dir,
  dedupe = TRUE
) {
  out_path <- file.path(output_dir,
    sprintf("%s_proteogenomics_%s.fasta",
            project_name, format(Sys.Date(), "%Y_%m")))

  components <- c(uniprot_fasta, predicted_orfs_fasta)
  if (!is.null(variant_fasta) && file.exists(variant_fasta)) {
    components <- c(components, variant_fasta)
  }
  if (!is.null(contaminants_fasta) && file.exists(contaminants_fasta)) {
    components <- c(components, contaminants_fasta)
  }

  if (dedupe) {
    tmp_concat <- tempfile(fileext = ".fasta")
    system2("cat", args = c(components, ">", tmp_concat))
    system2("seqkit", args = c("rmdup", "-s", "-o", out_path, tmp_concat))
    file.remove(tmp_concat)
  } else {
    system2("cat", args = c(components, ">", out_path))
  }

  # Count entries by source class
  composition <- count_proteog_classes(out_path)

  # Preserve the merged GTF alongside for coordinate lookups
  gtf_dest <- sub("\\.fasta$", "_merged.gtf", out_path)
  file.copy(merged_gtf, gtf_dest)

  register_proteogenomics_fasta(
    path = out_path,
    merged_gtf_path = gtf_dest,
    project_name = project_name,
    composition = composition
  )

  out_path
}

count_proteog_classes <- function(fasta_path) {
  # Parse source= tags from headers
  headers <- system2("grep", c("^>", fasta_path), stdout = TRUE)
  sources <- gsub(".*source=([A-Z_]+).*", "\\1", headers)
  list(
    total          = length(headers),
    REF            = sum(sources == "REF"),
    NOVEL_GENE     = sum(sources == "NOVEL_GENE"),
    NOVEL_ISOFORM  = sum(sources == "NOVEL_ISOFORM"),
    UNPARSED       = sum(sources == "UNPARSED"),
    UNIPROT        = sum(!grepl("source=", headers))  # canonical UniProt entries have no source= tag
  )
}
```

The registry now stores both the FASTA and the merged GTF path:

```json
{
  "my_experiment_2026_05": {
    "path": "/quobyte/proteomics-grp/de-limp/databases/proteogenomics/my_experiment_2026_05_proteogenomics_2026_05.fasta",
    "merged_gtf_path": "/quobyte/proteomics-grp/de-limp/databases/proteogenomics/my_experiment_2026_05_proteogenomics_2026_05_merged.gtf",
    "project_name": "my_experiment_2026_05",
    "organism": "Mus musculus",
    "reference_build": "GRCm39",
    "rnaseq_n_samples": 12,
    "rnaseq_total_reads": 348000000,
    "rrna_pct_mean": 4.2,
    "uniquely_mapped_pct_mean": 78.3,
    "read_length_tier": "default",
    "composition": {
      "total": 75432,
      "UNIPROT": 18203,
      "REF": 51289,
      "NOVEL_GENE": 5821,
      "NOVEL_ISOFORM": 119,
      "UNPARSED": 0
    },
    "created": "2026-05-20T14:52:00",
    "created_by": "brettsp",
    "pipeline_version": "1.1"
  }
}
```

---

## 8. Phase 5: DIA-NN Integration

Mostly unchanged from v1.0 — the FASTA appears in the pre-staged dropdown with a
🧬 tag, selecting it triggers an auto-warning notification. Two additions:

1. **Auto-recommend `--relaxed-prot-inf`** for the DIA-NN search since protein
   inference is harder with many similar novel-ORF sequences.
2. **Auto-bump min peptide length** from 7 to 8 since longer peptides have lower
   false-positive rates against the expanded search space.

```r
observeEvent(input$prestaged_fasta, {
  req(input$prestaged_fasta)
  registry <- load_proteog_registry()
  entry <- registry[[basename(input$prestaged_fasta)]]

  if (!is.null(entry)) {
    showNotification(
      tags$div(
        tags$p(strong("Proteogenomics database selected"), " for project ",
               tags$em(entry$project_name)),
        tags$p(sprintf(
          "Composition: %s UniProt + %s reference + %s novel genes + %s novel isoforms",
          format(entry$composition$UNIPROT, big.mark = ","),
          format(entry$composition$REF, big.mark = ","),
          format(entry$composition$NOVEL_GENE, big.mark = ","),
          format(entry$composition$NOVEL_ISOFORM, big.mark = ",")
        )),
        tags$p("Expected search time is roughly 3-8× longer than a standard UniProt search."),
        tags$p("After search completes, the Proteogenomics tab will activate and the ",
               tags$strong("Proteogenomics Glossary"), " will ship with the Claude export.")
      ),
      type = "message", duration = 12
    )

    updateTextInput(session, "diann_extra_flags",
      value = paste(input$diann_extra_flags, "--relaxed-prot-inf"))
    updateNumericInput(session, "min_pep_len", value = 8)
  }
})
```

---

## 9. Quality Gates

DE-LIMP halts the pipeline and surfaces a clear diagnostic if any of these fire:

| Gate | Stage | Threshold | Likely causes |
|------|-------|-----------|----------------|
| Species mismatch | Phase 1 | ENA metadata organism ≠ selected reference | wrong reference, wrong accession |
| Library type unsuitable | Phase 1 | strategy in {Tag-Seq, miRNA-Seq, OTHER} | use a different ingestion workflow |
| Read length too short | Phase 2.3 | median < 60 bp | re-sequence at PE100 or PE150 |
| rRNA contamination extreme | Phase 2.2 | >50% rRNA filtered | library prep failure; flag but allow proceed |
| Uniquely-mapped low | Phase 2.3a | <tier threshold | wrong reference, non-rRNA contamination, unsuitable library |
| Header parse failure | Phase 2.8 | any UNPARSED entries | pipeline bug; halt and surface |

Each gate produces a diagnostic block explaining (a) what was measured, (b) what
was expected, (c) likely causes ranked by frequency, (d) suggested next steps.
**Never silently produce a partial or noisy FASTA.**

---

## 10. UI Layout

(Unchanged from v1.0 in structure. See original spec §9. Add a "Read length tier"
indicator after fastp completes, and an "rRNA contamination %" indicator after
the rRNA filter completes. Both surface as informational chips; only the QC
gates halt the pipeline.)

---

## 11. Reactive Values

Unchanged from v1.0.

---

## 12. Files to Create / Modify

| File | Change |
|------|--------|
| `R/server_proteog_builder.R` | **New** — orchestration logic, status polling |
| `R/helpers_slims.R` | **New** — SLIMS URL scanning, ENA accession verification |
| `R/helpers_rnaseq.R` | **New** — sbatch generators with adaptive STAR thresholds |
| `R/helpers_proteog_assembly.R` | **New** — `assemble_proteogenomics_fasta()`, `count_proteog_classes()`, registry I/O |
| `R/helpers_proteog_qc.R` | **New** — `check_alignment_quality()`, gate enforcement |
| `R/server_search.R` | Modified — auto-warning observer when proteogenomics FASTA selected; uses extended `scan_prestaged_databases()` |
| `R/helpers_search.R` | Modified — `scan_prestaged_databases()` extended to read proteogenomics `registry.json` and emit 🧬-labeled choices |
| `R/ui.R` | Modified — convert top-level "New Search" `nav_panel` into `nav_menu("New Search", …)` dropdown containing **Run Search** + **Build Database 🧬** sub-panels (mirrors the Comparator pattern under Analysis). HPC-gated visibility for Build Database. |
| `app.R` | Modified — new `reactiveValues` for `is_proteogenomics`, `protein_classification`, `proteog_build_jobs`, `proteog_active_fasta` |
| `scripts/rewrite_transdecoder_headers.py` | **New** — header rewriter (production-grade). Mirror of the Hive-deployed copy at `/quobyte/proteomics-grp/de-limp/references/scripts/`. Not in `inst/` because DE-LIMP isn't an installed R package. |
| `scripts/proteogenomics_glossary.txt` | **New** — Glossary text shipped in Claude export ZIP; resolved relative to `app_dir` |
| `scripts/setup_references.sh` | **New** — admin script to pre-stage genome FASTAs, build rRNA indices |
| `CLAUDE.md` | Modified — add `"search_tab"` and `"build_database_tab"` to the "Tab values that MUST NOT change" protected list |

---

## 13. Dependencies

### R packages

| Package | Purpose | Likely present? |
|---------|---------|-----------------|
| `httr` | SLIMS/ENA HTTP | ✓ |
| `xml2` | ENA XML parsing | ✓ (used elsewhere) |
| `jsonlite` | Registry I/O | ✓ |
| `glue` | sbatch templating | ✓ |

No new R packages needed.

### Hive central modules (all verified available May 20, 2026)

```
fastp/0.23.4
bowtie2/2.5.2
star/2.7.11a
samtools/1.19.2
conda/stringtie/3.0.3   # NOT native stringtie/2.2.1 — that module has a bgzf bug
gffcompare/0.12.6
gffread/0.12.7
transdecoder/5.7.1
diamond/2.1.7           # already used by Cascadia
seqkit                  # via conda
conda/rseqc/5.0.4       # optional, for strand auto-detection
multiqc/1.33            # optional, for QC summary
```

### Custom conda env

One small env for the header rewriter:
```
/quobyte/proteomics-grp/de-limp/envs/proteog_helpers/
  python=3.11 biopython gffutils
```

Build with:
```bash
module load conda
conda create -p /quobyte/proteomics-grp/de-limp/envs/proteog_helpers \
  -c conda-forge -c bioconda --strict-channel-priority -y \
  python=3.11 biopython gffutils
```

The channel order matters. Reverse order causes the solver to pick python-2-only
biopython and fail.

### Reference data (pre-staged, admin-managed)

Available from bioinfocore-grp (read-only):
```
/quobyte/bioinfocore-grp/genomes/human/GRCh38.p14/
  ├── GCF_000001405.40_GRCh38.p14_genomic.gtf
  ├── STAR_2.7.11b_index/
  └── rRNA_human_03-12-2026.fasta
  (no genome FASTA — DE-LIMP downloads to its own tree)

/quobyte/bioinfocore-grp/genomes/mouse/GRCm39/
  ├── gencode.vM38.basic.annotation.gtf
  ├── STAR_GRCm39_vM38/
  ├── GRCm39.primary_assembly.genome.fa
  └── rRNA_mouse_03-12-2026.fasta
```

DE-LIMP's own tree (writable, populated by `setup_references.sh`):
```
/quobyte/proteomics-grp/de-limp/references/
  ├── registry.json
  ├── genomes/
  │   └── hg38_GRCh38.p14/genome.fa            # downloaded from NCBI
  ├── rrna_index/
  │   ├── hg38/                                 # bowtie2 indices we built
  │   └── mm39/
  └── (symlinks to bioinfocore-grp resources where appropriate)
```

### No Apptainer container needed

The original v1.0 spec recommended an Apptainer container bundling all tools.
Validation showed this is unnecessary — Hive's central modules cover everything,
and the modules are maintained by HPC@UCD rather than us. The `de-limp.sif`
container at `/quobyte/proteomics-grp/de-limp/containers/de-limp.sif` is
unchanged (still ships the DE-LIMP app and DIA-NN). The proteogenomics pipeline
itself is module-driven sbatch.

---

## 14. Cost and Time Estimates

For a typical experiment of 12 PE150 RNA-seq samples (~30M read pairs each), based
on validation data scaled up:

| Stage | Wall time | CPU-hours | Storage |
|-------|-----------|-----------|---------|
| SLIMS download (login node) | 30–60 min | 1 | 60 GB raw FASTQ |
| fastp (parallel array) | 20 min | 32 | +40 GB trimmed |
| bowtie2 rRNA filter (parallel) | 15 min | 24 | +35 GB filtered |
| STAR (parallel) | 60 min | 192 | +80 GB BAMs |
| stringtie per-sample (parallel) | 20 min | 32 | +2 GB GTFs |
| stringtie --merge | 10 min | 1 | +200 MB |
| gffcompare | 2 min | 0.2 | +50 MB |
| gffread | 5 min | 0.2 | +500 MB transcript FASTA |
| TransDecoder | 30 min | 8 | +800 MB ORFs |
| Header rewrite | 5 min | 0.2 | +200 MB final FASTA |
| Assembly + register | 5 min | 0.5 | (already counted) |
| **Total** | **~3-4 hours wall** | **~290 CPU-hours** | **~220 GB peak** |

After completion, intermediates can be auto-purged, leaving:
- Final FASTA (~200 MB depending on novel-ORF count)
- Preserved merged.gtf (~200 MB)
- QC report (~1 MB)
- Reproducibility log (~50 KB)

Total persistent storage per project: ~400 MB.

---

## 15. Testing Checklist

### Phase 1
- [ ] `scan_slims_url()` correctly parses a real SLIMS URL
- [ ] Non-SLIMS URL produces clear error
- [ ] R1/R2 pairing detection correct
- [ ] md5 verification runs and reports mismatches
- [ ] `verify_sra_accession()` correctly retrieves ENA metadata
- [ ] Species mismatch is caught BEFORE download starts
- [ ] Library type Tag-Seq/miRNA-Seq is refused with clear message
- [ ] ENA stream-subsample handles SIGPIPE correctly (set +o pipefail in sbatch)

### Phase 2
- [ ] fastp produces JSON with read length info
- [ ] Median read length detection from FASTQ works
- [ ] STAR tier selection picks correct flags for 92bp / 100bp / 150bp / 250bp reads
- [ ] rRNA bowtie2 index builds idempotently (skips if present)
- [ ] rRNA filtering keeps unmapped reads, drops mapped
- [ ] STAR completes with chosen flags
- [ ] **QC gate fires correctly** when uniquely-mapped is below threshold
- [ ] QC gate diagnostic mentions all 4 likely-cause categories
- [ ] gffcompare runs and produces class codes
- [ ] stringtie --merge inherits ref gene_ids for overlapping transcripts (expected)
- [ ] NOVEL_ISOFORM detection via gffcompare class codes `j` and `e` works
- [ ] gffread extracts transcripts cleanly
- [ ] TransDecoder produces .pep file
- [ ] Header rewriter produces **zero UNPARSED** entries
- [ ] Header field count is exactly 8 per entry (sanity check)

### Phase 3 (variant encoding, optional)
- [ ] Unchanged from v1.0

### Phase 4
- [ ] FASTA concatenation produces valid multi-FASTA
- [ ] seqkit rmdup removes sequence-identical entries
- [ ] Composition counts match source= tag counts
- [ ] Registry JSON includes merged_gtf_path field
- [ ] Older registry entries without merged_gtf_path don't crash on load

### Phase 5
- [ ] FASTA appears in pre-staged dropdown with 🧬 icon
- [ ] Selecting triggers notification with composition breakdown
- [ ] DIA-NN search runs successfully
- [ ] Results load with `is_proteogenomics = TRUE`
- [ ] Proteogenomics Glossary appears in Claude export ZIP

### Quality gates
- [ ] Pipeline halts on species mismatch
- [ ] Pipeline halts on too-short reads (<60bp)
- [ ] Pipeline halts on low uniquely-mapped rate
- [ ] Pipeline halts on header parse failures
- [ ] Halting produces clear, actionable diagnostic message
- [ ] No partial FASTAs produced after halt

### End-to-end
- [ ] Submit a real SLIMS URL or SRA accession from start
- [ ] Pipeline completes with all QC gates passed
- [ ] DIA-NN search against produced FASTA finds at least some novel-ORF peptides
- [ ] Methods text in Claude export correctly cites HISAT2 → STAR (Kim 2019),
      StringTie (Pertea 2015), gffcompare (Pertea 2020), TransDecoder, DIAMOND
- [ ] DNA Tech Core acknowledgment included (RRID:SCR_017740, NIH grant 1S10OD010786-01)

### Failure modes
- [ ] SLIMS URL expired (>1 month) → clear error
- [ ] Wrong reference selected → species check catches it OR low mapping triggers halt
- [ ] Network unreachable from Hive → fallback message
- [ ] Pipeline interrupted mid-array → resumable from last completed stage
- [ ] Out-of-disk during intermediates → graceful failure with cleanup suggestion

---

## 16. Validation Lessons Captured

All eleven lessons from the May 20, 2026 validation are now encoded as requirements:

1. ✅ **Conda channel ordering** — §13, env build command uses `-c conda-forge -c bioconda --strict-channel-priority`
2. ✅ **SRA accession verification** — §4 (Mode B), `verify_sra_accession()` displays metadata before download
3. ✅ **ENA streaming preferred** — §4, sra-toolkit is fallback only
4. ✅ **Reference completeness tiering** — §13, registry distinguishes complete vs index-only references
5. ✅ **rRNA pre-filtering mandatory** — §5.2, bowtie2 step inserted between fastp and STAR
6. ✅ **Adaptive STAR thresholds** — §5.3, `select_star_params()` with four tiers
7. ✅ **Unique-mapping QC gate** — §5.3a, halts pipeline below tier threshold
8. ✅ **StringTie module choice** — §13, uses `conda/stringtie/3.0.3`, not native module
9. ✅ **gffcompare for novel isoforms** — §5.5a, enables NOVEL_ISOFORM detection
10. ✅ **Header rewriter is integration linchpin** — §5.8, format validated, 100% parse-clean required
11. ✅ **No Apptainer container needed** — §13, all tools available as Hive modules

---

## Methods boilerplate (auto-generated)

For the Claude export Manuscript template, append to Methods generator:

```r
generate_proteogenomics_methods <- function(build_metadata) {
  glue::glue("
RNA-sequencing data for proteogenomics database construction was generated at the
UC Davis DNA Technologies Core (RRID:SCR_017740) on an Illumina {build_metadata$instrument}
platform using {build_metadata$library_type} library preparation.
Sequencing reads were quality- and adapter-trimmed using fastp v0.23.4
(Chen et al., 2018), then filtered against an organism-specific ribosomal RNA
reference using bowtie2 v2.5.2 (Langmead and Salzberg, 2012) with the
--very-sensitive-local preset. Remaining reads were aligned to the
{build_metadata$genome_version} reference genome with STAR v2.7.11a (Dobin et al., 2013)
using read-length-appropriate thresholds (tier: {build_metadata$tier}). Aligned reads
were assembled into per-sample transcripts with StringTie v3.0.3 (Pertea et al., 2015),
then merged across samples to produce a unified transcript catalog. Merged transcripts
were classified against the reference annotation using gffcompare v0.12.6 (Pertea and
Pertea, 2020), enabling identification of both novel intergenic loci and novel splice
isoforms of known genes. Open reading frames were predicted from the merged transcripts
using TransDecoder v5.7.1, with --single_best_only retaining one ORF per transcript.
Predicted protein sequences were concatenated with the {build_metadata$uniprot_release}
UniProt reference proteome and cRAP contaminant sequences to form the proteogenomics
search database ({build_metadata$total_entries} total entries: {build_metadata$canonical}
canonical, {build_metadata$ref_novel} reference-derived predicted, {build_metadata$novel_gene}
novel genes, {build_metadata$novel_isoform} novel isoforms). The expanded database was
used as input to DIA-NN as described above.
")
}
```

References:
- Chen S et al. (2018). fastp. *Bioinformatics* 34:i884–i890.
- Langmead B, Salzberg SL (2012). bowtie2. *Nat Methods* 9:357–359.
- Dobin A et al. (2013). STAR. *Bioinformatics* 29:15–21.
- Pertea M et al. (2015). StringTie. *Nat Biotechnol* 33:290–295.
- Pertea G, Pertea M (2020). GFF Utilities: GffRead and GffCompare. *F1000Research* 9:304.
- Haas BJ. TransDecoder. https://github.com/TransDecoder/TransDecoder
- Buchfink B et al. (2021). DIAMOND. *Nat Methods* 18:366–368.
- Acknowledgment: "RNA sequencing was performed at the DNA Technologies and
  Expression Analysis Core at the UC Davis Genome Center (RRID:SCR_017740),
  supported by NIH Shared Instrumentation Grant 1S10OD010786-01."

---

*Spec version 1.1 — Brett Phinney / UC Davis Proteomics Core — May 2026 (post-validation)*
