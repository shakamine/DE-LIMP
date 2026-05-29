# ASMS Poster Content Brief — De Novo Sequencing & Proteogenomics in DE-LIMP

> **Purpose:** Source material for designing an ASMS conference poster.
> **For:** Claude (poster design / layout).
> **From:** Brett Phinney, UC Davis Proteomics Core.
> **Audience:** ASMS attendees — mass spectrometrists, proteomics core staff, method developers. Technical, but the poster should still read in 60 seconds from 6 feet away.
> **Level:** Highlights — headline story + key results + figure ideas, not exhaustive method detail.
> **Date:** 2026-05-22

---

## 0. The one-sentence story

**DE-LIMP** — an open-source Shiny platform for DIA proteomics — now closes the gap on peptides that conventional reference-database search misses, via two complementary additions: **sample-matched proteogenomics database building** (expand the search space from the sample's own RNA-seq) and **de novo sequencing** (identify peptides with no reference at all), both orchestrated on HPC from a point-and-click interface.

---

## 1. Suggested poster title options

- *"Beyond the Reference: Integrated Proteogenomics and De Novo Sequencing in DE-LIMP, an Open-Source DIA Analysis Platform"*
- *"Finding the Peptides You're Missing: Sample-Matched Proteogenomic Databases and De Novo Sequencing in a Point-and-Click DIA Pipeline"*
- *"DE-LIMP: From DIA-NN to Novel Peptides — Proteogenomics and De Novo Sequencing Without Leaving the GUI"*

(Authors / affiliation: Brett Phinney et al., UC Davis Proteomics Core / Genome Center. Fill in co-authors.)

---

## 2. Motivation / Introduction (poster intro block)

Standard DIA workflows (DIA-NN, Spectronaut) identify peptides **only if they already exist in a reference proteome.** Two large classes of real peptides are invisible to this approach:

1. **Sample-specific sequences** — novel genes, novel isoforms, and variants that aren't in the public reference for that organism/sample.
2. **Anything from an organism with no good reference at all** — non-model species, paleoproteomics (feathers, fossils), environmental samples.

DE-LIMP adds two orthogonal solutions, both runnable by a non-bioinformatician through the existing GUI, with all heavy compute pushed to an HPC SLURM cluster:

- **Proteogenomics DB Builder** — *expand* the database from the sample's own transcriptome.
- **De Novo Sequencing** — *eliminate* the need for a database.

**Hook for the poster:** these are usually bioinformatics-heavy, command-line workflows. DE-LIMP makes them a few clicks in a web app, with self-describing methods output for reproducibility.

---

## 3. Feature 1 — Proteogenomics Database Builder

### Concept
Submit the *same biological samples* to both the Proteomics Core (MS) and the DNA Technologies Core (RNA-seq). Build a custom search FASTA = **reference proteome + novel ORFs translated from that sample's own RNA-seq.** Search the DIA data against the expanded FASTA → recover peptides absent from any public reference.

### Pipeline (HPC, 11-stage SLURM dependency chain)
```
fastp → rRNA filter → STAR align → QC gate → StringTie assemble
   → merge → gffcompare → gffread → TransDecoder (ORF calling)
   → header rewrite → assemble FASTA
```
Every ORF is tagged in its FASTA header by source class:
**REF · NOVEL_GENE · NOVEL_ISOFORM · VARIANT · UNIPROT**

### Headline result (validated end-to-end on UC Davis "Hive" HPC, May 2026)
- **67,386-entry custom FASTA** produced from a mouse test dataset:
  - **66,046 reference** entries
  - **1,340 novel-gene** entries discovered from sample RNA-seq
- All 11 pipeline stages verified; species-aware quality gates enforced.

### Method highlights worth calling out (these are the "we learned this the hard way" credibility points)
- **STAR** aligner (proteogenomics-community standard; pre-staged indices).
- **Mandatory rRNA pre-filter** — without it, validation data showed 73% multi-mapping.
- **Adaptive STAR thresholds by read length** — defaults assume ≥150 bp; pipeline detects read length and relaxes thresholds for shorter reads (tiers at 130 / 100 / 60 bp).
- **gffcompare step** distinguishes novel *isoforms* of known genes from collapsing into reference (the difference between "1,340 novel genes" and "1,340 novel genes + N novel isoforms").
- **Quality gates that HALT, not guess:** species mismatch (pre-alignment), uniquely-mapped rate below read-length-tiered threshold, malformed headers. On failure the pipeline surfaces candidate causes (wrong reference / contamination / wrong library type / degraded sample) rather than silently proceeding.

### Suggested figure
- **Pipeline flow diagram** (the 11 stages, with the FASTA-class legend) — the centerpiece schematic.
- **Stacked bar / donut** of FASTA composition: 66,046 REF vs 1,340 NOVEL_GENE.
- (If available) a volcano or table showing example novel-gene peptides recovered in the MS search that were absent from reference-only search.

---

## 4. Feature 2 — De Novo Sequencing (reference-free)

### Concept
De novo sequencing reads peptide sequence **directly from the spectrum** — no database. Two uses: (1) orthogonal validation of database-search IDs, and (2) the killer app — **identifying samples with no reference proteome.**

### Engines integrated (multi-engine, multi-acquisition)
- **Cascadia** — de novo for **DIA** data, including timsTOF ion-mobility.
- **Casanovo** — de novo for **DDA** data (Orbitrap + timsTOF).
- **Sage** — fast DDA database search, run alongside for comparison.
- **DIAMOND BLAST** — maps de novo peptides to species/proteins (vs SwissProt).

All run as **GPU SLURM jobs**; results auto-downloaded and visualized in 12 dedicated panels (tables, BLAST/species breakdowns, FDR, cross-species Venn, coverage maps, modification tracking).

### Flagship application — paleoproteomics (great poster visual)
**Bird-feather and fossil species identification with no reference proteome:**
- Workflow: .raw → msconvert → mzML → **Sage (91k PSMs)** + MGF → **Casanovo (28k PSMs)** → classify → **DIAMOND BLAST (56k hits)**.
- Top species hits on feather data: **pigeon (COLLI), chicken (CHICK), mallard duck (ANAPL)** — biologically sensible.
- Real samples: Bonaparte's Gull & Whooping Crane feathers, *Aratasaurus* fossil collagen.
- **Deamidation (N/Q) tracking** included as a paleoproteomics authenticity signal.

### Benchmarking results worth a panel
- **Sage vs DIA-NN, feather Orbitrap:** 91,238 vs ~14,000 PSMs; comparable protein-group counts (1,421 vs 1,409) despite 6× PSM difference.
- **Sage vs DIA-NN, timsTOF HeLa DDA:** 80,040 vs 42,000 PSMs (Sage wins on timsTOF DDA).

### Model-development angle (optional "methods development" panel — strong for ASMS)
We fine-tuned de novo models for **timsTOF ion mobility (IM)**:
- Demonstrated **IM is highly discriminative** (8.79× discriminative ratio between sequences) — strong rationale for IM-aware de novo.
- Found that **uniform-LR fine-tuning fails** to learn IM features (89.7% of IM weights stayed ~0) and causes **catastrophic forgetting** of Orbitrap capability — motivating a **differential-LR / zero-init** training recipe. (Honest negative + positive result; ASMS audiences appreciate this.)

### Suggested figures
- **Feather-to-species workflow schematic** (.raw → de novo → BLAST → species call) — the narrative figure.
- **Species breakdown bar / cross-species Venn** for a feather sample.
- **Sage vs DIA-NN PSM comparison bar chart** (the 91k vs 14k contrast is visually punchy).
- **Per-residue confidence heatmap** on an example peptide (shows the de novo confidence concept).

---

## 5. How the two features fit together (one unifying schematic)

```
   RNA-seq of same sample ──► Proteogenomics DB Builder ──► expanded FASTA ──┐
                                                                              ▼
                                                                       DIA-NN search ──► DE-LIMP DE analysis
                                                                              │
   MS spectra ───────────────────────────────────────────────► De Novo Sequencing (no DB needed)
                                                                              │
                                                                              ▼
                                                                  novel / no-reference peptides
```
**Unifying message:** both recover peptides the standard reference search misses — one by *expanding the database from the sample's own RNA*, the other by *needing no database at all*.

---

## 6. Platform / "by the numbers" sidebar (optional poster corner)

DE-LIMP is open-source, single-developer, built rapidly with AI-assisted coding:
- **700 commits over ~16 weeks** (Jan–May 2026), active 81% of days.
- **~47K lines of R**, 39 GitHub releases.
- Deployed on GitHub + Hugging Face; runs locally, in Docker, or on HPC (Apptainer + SLURM proxy).
- (Two pre-made poster figures already exist in UC Davis Aggie Blue/Gold: a commit-activity timeline and a "by the numbers" stat panel — `~/Downloads/DE-LIMP_commits_timeline.png` and `DE-LIMP_stat_tiles.png`.)

> Use this only if there's room — it's a "meta" angle (rapid open-source tool development) that some ASMS software sessions like, but the science (§3–§4) is the main event.

---

## 7. Conclusions (poster wrap-up bullets)

- DE-LIMP brings **proteogenomic database construction** and **de novo sequencing** into a single point-and-click DIA platform, with HPC compute hidden behind the GUI.
- Sample-matched proteogenomics recovered **1,340 novel-gene entries** beyond the reference on a validated mouse run.
- De novo sequencing enables **species identification with no reference proteome** (paleoproteomics: feathers, fossils), with multi-engine (Cascadia/Casanovo/Sage) cross-validation.
- All outputs are **self-describing for reproducibility** — methods text is generated from the pipeline object, not hardcoded.
- Free, open-source, and deployable from laptop to HPC.

---

## 8. Practical notes for the poster designer (Claude)

- **Two columns ≈ two features**, joined by the unifying schematic (§5) as the visual bridge.
- The **strongest single visuals**: the proteogenomics 11-stage pipeline diagram, and the feather→species de novo workflow. Lead with these.
- **Punchy numbers to make large:** 67,386 FASTA entries (1,340 novel) · 91k vs 14k PSMs (Sage vs DIA-NN) · 56k BLAST hits · 12 de novo result panels.
- UC Davis palette: **Aggie Blue `#022851`**, **Aggie Gold `#FFBF00`**.
- Keep jargon captioned — even an ASMS audience spans method-developers to core-facility generalists.
- Ask me for any specific result table, screenshot, or figure export and I can generate it.
```
```
