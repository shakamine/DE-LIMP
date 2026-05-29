# ASMS Poster Content Brief — Run Comparator & Multi-Omics Integration in DE-LIMP

> **Purpose:** Source material for designing an ASMS conference poster.
> **For:** Claude (poster design / layout).
> **From:** Brett Phinney, UC Davis Proteomics Core.
> **Audience:** ASMS attendees — mass spectrometrists, proteomics core staff, method developers.
> **Level:** Highlights — headline story + key results + figure ideas.
> **Date:** 2026-05-22
> **Companion brief:** `DESIGN_BRIEF_DENOVO_PROTEOGENOMICS.md` (de novo + proteogenomics). These two could be one poster or two.

---

## 0. The one-sentence story

**DE-LIMP** adds two analyst-facing capabilities on top of its DIA differential-expression core: a **Run Comparator** that diagnoses *why two software pipelines disagree* on the same dataset (DE-LIMP vs Spectronaut / FragPipe / a second DE-LIMP run), and **MOFA2 multi-omics integration** that finds shared and view-specific variation across up to six data layers — both with AI-assisted interpretation built in.

---

## 1. Suggested poster title options

- *"Why Don't My Tools Agree? An Automated Cross-Software Diagnostic Comparator for DIA Proteomics"*
- *"From Disagreement to Diagnosis: Comparing DIA-NN, Spectronaut and FragPipe in DE-LIMP"*
- *"Beyond a Single Pipeline: Cross-Tool Comparison and Multi-Omics Factor Integration in DE-LIMP"*

(Authors / affiliation: Brett Phinney et al., UC Davis Proteomics Core / Genome Center — fill in co-authors.)

---

## 2. Feature 1 — Run Comparator

### The problem (poster intro hook)
Run the *same* DIA dataset through DIA-NN and through Spectronaut and you get **different protein counts, different precursors, and different DE calls.** Core facilities field this question constantly: *"which tool is right, and why do they disagree?"* Answering it normally means hours of manual cross-referencing. The Run Comparator automates the diagnosis.

### What it does
Loads two analyses of the **same samples** and aligns them protein-by-protein and precursor-by-precursor, then reports *where* and *why* they differ. Three modes:
- **DE-LIMP vs Spectronaut** (parses the full Spectronaut export ZIP — Pivot + Normal report + Candidates + RunSummaries).
- **DE-LIMP vs FragPipe** (FragPipe-Analyst or raw `combined_protein.tsv`).
- **DE-LIMP vs DE-LIMP** (same pipeline, different parameter settings).

### Diagnostic layers
1. **Settings Diff** — every pipeline parameter side-by-side (FDR cutoffs, normalization, MBR, library type, TopN). Color-coded by severity; flags structurally different and *problematic* settings (e.g. Spectronaut's "Use All MS-Level Quantities" / Quant3, which silently doubles t-test sample size).
2. **Protein Universe** — Venn + per-tool ID stats (proteins / precursors / peptides, totals and per-sample), plus **per-sample Venn diagrams** of protein overlap.
3. **Quantification** — intensity correlation on shared proteins, systematic-bias detection, TopN-limitation effect.
4. **DE Concordance** — a 3×3 (Up / Down / NS) concordance matrix, a hypothesis engine (9 rules) that assigns a likely *cause* to each discordant protein, and a "rescue" analysis for proteins untestable in one tool but significant in the other.

### Worked example (real result this session — dog plasma DIA, 60 SPD timsTOF)
Comparing DE-LIMP (DIA-NN/DPC-Quant) vs Spectronaut on the same 8-sample dataset:
- **Proteins:** 2,205 (DE-LIMP) vs 2,456 (Spectronaut) — **+11.4%**
- **Precursors:** 24,528 vs 32,111 — **+30.9%**
- **Peptides:** 20,488 vs 23,248 — **+13.5%**
- **Median peptides/protein:** 4.0 vs 2.0
- **1,896 shared proteins.**

**Interpretation the tool surfaces:** Spectronaut's larger precursor count (+31%) doesn't translate into proportionally more proteins (+11%) — its median peptides-per-protein is *half* DE-LIMP's. So the extra precursors spread thin across more single/double-peptide protein hits rather than deepening coverage. The comparator gives the analyst the primary per-precursor data (charge, q-value, proteotypic status) to test whether those extras are near-threshold IDs, extra charge states, or genuine new coverage.

### AI-assisted interpretation (a distinctive angle for ASMS)
The comparator exports a **primary-data bundle + an instruction-rich prompt** for an LLM (Claude/Gemini): per-precursor long tables for both tools, protein-intensity matrices, settings diff, and a prompt that lists specific analyses to run and "smoking-gun patterns" to check (settings mismatch, near-threshold IDs, charge-coverage gap, single-peptide hits, sample dropouts). Design philosophy: **ship primary data + how-to-analyze, not pre-baked summaries** — so the model can answer questions we didn't anticipate.

### Suggested figures
- **The Venn + ID-stats panel** (2,205 vs 2,456, with the per-sample bar chart of proteins/precursors/peptides) — the headline visual.
- **Settings Diff table excerpt** showing a flagged "severe" row (e.g. Quant3 enabled) — makes the diagnostic concrete.
- **3×3 DE concordance matrix** with the hypothesis-engine legend.
- **Median peptides-per-protein contrast** (4.0 vs 2.0) as a simple, punchy two-bar figure — it's the crux of the dog-data story.

---

## 3. Feature 2 — Multi-Omics Integration (MOFA2)

### What it does
Wraps **MOFA2** (Multi-Omics Factor Analysis) so an analyst can integrate up to **6 data views** (proteomics + transcriptomics + metabolomics + phospho + …) and discover **latent factors** — unsupervised axes of variation that are either *shared* across views or *view-specific*. It's PCA generalized to multiple omics layers, telling you which biological signal lives in which data type.

### How it works in DE-LIMP
- The current proteomics result auto-registers as the first view ("Global Proteomics").
- Add up to 6 views via a card UI; each view gets a data-type tag (continuous / count / etc.) and its own upload or in-app source.
- "Scale views" option (recommended ON) equalizes contributions across omics layers.
- Trains MOFA2 in an **isolated subprocess** (basilisk/Python) so a Python crash can't take down the Shiny app — a robustness detail worth a footnote.
- Two built-in example datasets for demo: **Mouse Brain (2-view)** and **TCGA Breast (3-view)**.

### Result panels
- **Variance-explained heatmap** — factor × view, the core MOFA readout (which factor is shared vs view-specific).
- **Factor weights plot** — top features driving each factor.
- **Factor scores plot** — samples in latent factor space (clustering / group separation).
- **Top features table** — the loadings, exportable.
- **Factor ↔ DE correlation** — links unsupervised MOFA factors back to the supervised differential-expression contrasts, so a factor can be tied to the experimental design.

### Cross-feature reuse (nice integration story)
The Run Comparator can *also* invoke MOFA2: it treats Run A and Run B as two views and decomposes their joint vs tool-specific variance — using the multi-omics engine to quantify how much of the variation is "real biology" vs "tool artifact." One engine, two uses.

### Suggested figures
- **Variance-explained heatmap** (factor × view) — the signature MOFA figure; lead with it.
- **Factor scores scatter** colored by experimental group — shows the integration separating conditions.
- **The view-card UI** screenshot — communicates "up to 6 omics, point-and-click."

---

## 4. Unifying message (if these share a poster)

Both features are about **interpretation, not just identification**:
- The **Comparator** answers *"do my tools agree, and if not, why?"* — turning a black-box discrepancy into a ranked set of causes.
- **MOFA2** answers *"where does my signal live across omics layers?"* — turning multiple data tables into a small set of interpretable factors.

Both also lean on **AI-assisted analysis** (Comparator's LLM export prompt; the broader DE-LIMP AI Summary feature), positioning DE-LIMP as a platform that helps the analyst *reason about* results, not just produce them.

---

## 5. Platform / "by the numbers" sidebar (optional, shared with companion poster)

- Open-source, single-developer Shiny platform: **700 commits / ~16 weeks**, **~47K lines of R**, 39 releases.
- Deploys local / Docker / HPC (Apptainer + SLURM proxy); GitHub + Hugging Face.
- Pre-made UC Davis Aggie Blue/Gold poster figures exist (`~/Downloads/DE-LIMP_commits_timeline.png`, `DE-LIMP_stat_tiles.png`).

---

## 6. Conclusions (poster wrap-up bullets)

- The **Run Comparator** automates cross-software diagnosis (DIA-NN vs Spectronaut vs FragPipe), aligning data from protein down to precursor level and assigning likely *causes* to every disagreement.
- On real dog-plasma timsTOF data it showed Spectronaut's +31% precursors yielding only +11% proteins — extra IDs spread across shallow-coverage proteins, not deeper coverage.
- **MOFA2 integration** brings unsupervised multi-omics factor analysis (up to 6 views) into the same GUI, with factors linked back to differential-expression contrasts.
- Both features emphasize **interpretation**, including AI-assisted analysis via structured primary-data exports.
- Free, open-source, laptop-to-HPC.

---

## 7. Practical notes for the poster designer (Claude)

- **Two columns ≈ two features.** The Comparator is the stronger, more novel story — give it the larger share; MOFA2 is the complementary "and we also integrate omics" panel.
- **Strongest single visuals:** the Comparator's Venn + ID-stats + per-sample bars, and MOFA2's variance-explained heatmap.
- **Punchy numbers to enlarge:** 2,205 vs 2,456 proteins · +30.9% precursors · median 4.0 vs 2.0 peptides/protein · 9-rule hypothesis engine · up to 6 omics views.
- UC Davis palette: **Aggie Blue `#022851`**, **Aggie Gold `#FFBF00`**.
- The dog-plasma comparison is real and reproducible in-app — I can export any specific figure (Venn, per-sample bars, settings-diff excerpt, concordance matrix) on request.
- Caption jargon (DPC-Quant, MBR, Quant3, MOFA factor) even for an ASMS crowd.
