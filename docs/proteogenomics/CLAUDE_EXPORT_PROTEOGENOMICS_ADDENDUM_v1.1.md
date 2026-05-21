# Claude Export Addendum: Proteogenomics Context Injection

> **Version**: 1.1 — May 2026 (post-validation; codebase-reconciled May 21, 2026)
> **Addendum to**: `CLAUDE_EXPORT_PROMPTS_SPEC.md`
> **Author**: Brett Phinney / UC Davis Proteomics Core
> **Prereqs**: `CLAUDE_EXPORT_PROMPTS_SPEC.md` implemented (all three templates)
> **Pairs with**: `PROTEOGENOMICS_DB_BUILDER_SPEC.md` v1.1
> **Status**: Header format validated end-to-end May 20, 2026 against real
> pipeline output (67,386 entries, 100% parse-clean). Glossary file location
> and `write_proteogenomics_glossary()` signature reconciled against actual
> DE-LIMP layout (not an installed R package) on May 21, 2026 — glossary lives
> at `scripts/proteogenomics_glossary.txt`, resolved via `app_dir`.

---

## Changes from v1.0

This revision aligns the Glossary content and prompt blocks with the actual header
format produced by the validated pipeline.

1. **Header format updated** — replaces hypothetical Jagtap-style `_u_/_c_/_i_/_o_`
   suffix codes with the validated `sp|ID|SYMBOL_TAG source=REF/NOVEL_GENE/NOVEL_ISOFORM/UNPARSED ...` 
   format. All metadata fields (source, ORF_type, strand, len, coords, parent_gene,
   transcript) are explicit key=value pairs, no decoding needed.
2. **ORF_type confidence tiering added** — Glossary now distinguishes `complete`
   (high confidence) from `5prime_partial` / `3prime_partial` (lower confidence)
   from `internal` (lowest confidence). This was implicit in v1.0 but now explicit.
3. **NOVEL_ISOFORM class added** — v1.0 didn't include this because the pipeline
   plan didn't yet have the gffcompare step. v1.1 of the builder spec adds it,
   so the Glossary must describe it.
4. **Variant proteoform discussion deferred** — INDEL_ENSP*/SNV_ENSP* entries
   come from the optional variant-encoding Phase 3. Most builds won't have them.
   Glossary text now treats them as conditional rather than co-equal with novel ORFs.
5. **Classifier function updated** — `classify_protein_id()` now reads the explicit
   `source=` tag from the header description rather than inferring class from the
   accession prefix. Simpler and more robust.

---

## Overview

When the loaded session was searched against a proteogenomics-expanded FASTA, the
protein table contains identifiers that Claude — and most readers — won't recognize:

```
sp|ENSMUST00000000001.5.p1|Gnai3_MM39TEST         ← REF (reference protein)
sp|MSTRG.10029.5.p2|MSTRG.10029_MM39TEST          ← NOVEL_GENE (StringTie novel locus)
sp|MSTRG.10075.2.p1|Trim25_MM39TEST               ← NOVEL_ISOFORM (novel isoform of known gene)
INDEL_ENSP00000354813_81:CAAAAAAAACTC_CAAAAAAACTC ← Variant proteoform (optional Phase 3)
```

If Claude generates a "Top DE Proteins" section without knowing what these are, it
will either invent a plausible-sounding gene function (hallucination risk) or
silently treat them as ordinary UniProt entries (interpretation risk). Neither is
acceptable for a core facility deliverable.

This addendum adds:

1. **Detection** — flag the session as proteogenomic at load time
2. **Classification** — tag every protein with its `source` class
3. **A new export file** `Proteogenomics_Glossary.txt` shipped in the ZIP
4. **Conditional prompt sections** in all three templates explaining how to
   handle non-canonical identifiers
5. **Inline summaries** of proteogenomic discoveries appended to `{ctx}`

No changes are required to the three top-level template structures from the parent
spec — only new conditional blocks slotted into existing placeholders.

---

## Detection: when does this addendum activate?

Add to `app.R` reactive values:

```r
is_proteogenomics      = FALSE,  # TRUE if FASTA contains source= tags or variant prefixes
protein_classification = NULL,   # data.frame: Protein.Group, source, orf_type, parent_gene
```

Set after `readDIANN()` runs:

```r
values$protein_classification <- classify_proteins(values$diann_report)
values$is_proteogenomics <- any(
  values$protein_classification$source %in% c("REF", "NOVEL_GENE", "NOVEL_ISOFORM")
) || any(grepl("^(INDEL|SNV)_ENSP", values$protein_classification$Protein.Group))
```

Helper (new file `R/helpers_proteogenomics.R`):

```r
classify_proteins <- function(diann_report) {
  # DIA-NN's report.parquet preserves FASTA descriptions in the Protein.Names
  # or Protein.Group.Description column depending on version. Parse source=tags.
  
  ids <- unique(diann_report$Protein.Group)
  descriptions <- diann_report$Protein.Group.Description[
    match(ids, diann_report$Protein.Group)
  ]
  
  # Extract source= tag from description (validated format)
  source_tag <- gsub(".*source=([A-Z_]+).*", "\\1", descriptions)
  source_tag[!grepl("source=", descriptions)] <- "UNIPROT"
  
  # Extract ORF_type (NA for UNIPROT entries)
  orf_type <- gsub(".*ORF_type=(\\w+).*", "\\1", descriptions)
  orf_type[!grepl("ORF_type=", descriptions)] <- NA_character_
  
  # Extract parent_gene (NA for UNIPROT entries)
  parent_gene <- gsub(".*parent_gene=(\\S+).*", "\\1", descriptions)
  parent_gene[!grepl("parent_gene=", descriptions)] <- NA_character_
  
  # Detect variant proteoforms by accession prefix (Phase 3 output)
  is_variant <- grepl("^(INDEL|SNV)_ENSP", ids)
  source_tag[is_variant] <- "VARIANT"
  
  data.frame(
    Protein.Group = ids,
    source        = source_tag,
    orf_type      = orf_type,
    parent_gene   = parent_gene,
    stringsAsFactors = FALSE
  )
}
```

The classifier reads explicit metadata from the FASTA — no regex on accession prefix
required (except for the optional VARIANT class). This is more robust than v1.0's
prefix-based inference.

If `is_proteogenomics` is FALSE, this addendum is a no-op — all conditional blocks
below resolve to empty strings.

---

## New shipped file: `Proteogenomics_Glossary.txt`

Written into the ZIP only when `values$is_proteogenomics` is TRUE. Plain-text
reference Claude consults while generating the report.

```
PROTEOGENOMICS GLOSSARY
============================================================================

This experiment used an expanded "proteogenomics" search database that combines
standard reference proteins (UniProt) with sample-specific predicted proteins
derived from matched RNA-seq. Every protein entry has a `source` tag that
identifies which class it belongs to.

----------------------------------------------------------------------------
HOW TO READ A PROTEOGENOMICS IDENTIFIER
----------------------------------------------------------------------------

The accession follows UniProt format:  sp|<protein_id>|<symbol>_<PROJECT_TAG>

The description contains metadata as key=value pairs:
   source=...           one of REF, NOVEL_GENE, NOVEL_ISOFORM, UNPARSED, VARIANT
   ORF_type=...         complete | 5prime_partial | 3prime_partial | internal
   strand=...           + or -
   len=...              ORF length in amino acids
   coords=...           transcript-relative coordinates of the ORF
   parent_gene=...      gene ID this ORF belongs to (ENSMUSG/ENSG or MSTRG)
   transcript=...       transcript ID this ORF was predicted from

Example header:
   sp|MSTRG.10029.5.p2|MSTRG.10029_MM39TEST source=NOVEL_GENE
      ORF_type=5prime_partial strand=- len=112
      coords=MSTRG.10029.5:344-682(-)
      parent_gene=MSTRG.10029 transcript=MSTRG.10029.5

PROJECT_TAG (e.g., "MM39TEST") identifies which proteogenomics build this entry
came from. Useful for tracking when multiple builds are loaded.

----------------------------------------------------------------------------
THE FOUR PROTEIN CLASSES
----------------------------------------------------------------------------

source=REF
   Standard reference protein. The protein_id is an Ensembl transcript ID
   (e.g., ENSMUST00000000001.5.p1). The symbol (e.g., "Gnai3") is the
   standard gene symbol. These are well-characterized proteins from the
   reference proteome and should be interpreted normally — same as a hit
   in any standard UniProt-only search.

source=NOVEL_GENE
   A protein predicted from a StringTie transcript that does NOT overlap
   any annotated gene — a putative novel protein-coding region discovered
   from sample-specific RNA-seq. The protein_id begins with "MSTRG.".
   
   These are CANDIDATE novel ORFs. They require orthogonal validation
   before biological claims can be made. Single-peptide identifications
   are the norm, not the exception.

source=NOVEL_ISOFORM
   A protein predicted from a transcript that overlaps a known gene but
   represents a novel splicing pattern (alternative isoform not in the
   reference annotation). These inherit the parent gene's symbol but
   carry an alternative protein sequence.
   
   For proteogenomics, novel isoforms are often more biologically
   interesting than novel genes — alternative splicing is a well-
   characterized mechanism for regulating protein function. A novel
   isoform with differential expression vs the canonical form is a
   strong biological signal.

source=UNPARSED
   Header could not be parsed by the rewriter — should not appear in a
   clean pipeline run. If present, treat as a data-quality flag.

source=VARIANT  (only if optional Phase 3 ran)
   Variant proteoform of a canonical reference protein. The protein_id
   begins with INDEL_ENSP* or SNV_ENSP*, followed by the parent Ensembl
   protein ID, codon position, and reference→alternate sequence.
   
   Example: INDEL_ENSP00000354813_81:CAAAAAAAACTC_CAAAAAAACTC
   = ENSP00000354813 with a 1-nucleotide deletion at codon 81.

----------------------------------------------------------------------------
ORF_TYPE CONFIDENCE TIERING
----------------------------------------------------------------------------

Among NOVEL_GENE and NOVEL_ISOFORM entries, ORF_type indicates structural
confidence:

   complete           — Full ORF with start (ATG) and stop codon present.
                        Highest confidence for a real protein. Treat as
                        the strongest candidates for biological discussion.
   
   5prime_partial     — Missing N-terminus (no start codon). Could be an
                        alternative start site OR a truncated transcript
                        assembly. Moderate confidence.
   
   3prime_partial     — Missing C-terminus (no stop codon). Usually
                        indicates truncated assembly. Lower confidence.
   
   internal           — Missing both ends. Lowest confidence — typically
                        a fragment of a longer transcript that wasn't
                        fully assembled.

When generating "top discoveries" lists, weight by ORF_type confidence.
A `complete` ORF with modest fold-change is more interpretable than an
`internal` ORF with high fold-change.

----------------------------------------------------------------------------
GENOMIC COORDINATES
----------------------------------------------------------------------------

The coords= field gives TRANSCRIPT-relative coordinates only. GENOMIC
coordinates (which chromosome, which position) can be looked up from
the merged.gtf file that ships alongside the FASTA, using the parent_gene
field as the lookup key.

DE-LIMP's Proteogenomics tab provides this lookup automatically. When
discussing a NOVEL_GENE in the report, you may reference the parent_gene
ID but do NOT invent a chromosome or position — defer to the user's tool
to surface the actual coordinates.

----------------------------------------------------------------------------
INTERPRETATION GUIDANCE FOR THE REPORT
----------------------------------------------------------------------------

For REF and UNIPROT entries:
   - Interpret normally. These are annotated proteins with known function.

For NOVEL_GENE entries:
   - Do NOT invent biological function. By definition these have no
     prior annotation.
   - Refer to them as "candidate novel proteins" or "novel ORF
     candidates" — never as established proteins.
   - When listing top hits, prefer `complete` ORFs over partials.
   - Recommend orthogonal validation: targeted MS, de novo sequencing,
     matched RNA-seq quantification, or western blot using a peptide-
     specific antibody.

For NOVEL_ISOFORM entries:
   - The parent gene's known function is relevant context — you may
     reference it.
   - But emphasize that this is a NOVEL ISOFORM with potentially
     ALTERED function, stability, or localization compared to the
     canonical form.
   - When the canonical reference protein is ALSO in DE_Results_Full.csv,
     comparing fold-changes of canonical vs novel isoform is the most
     informative framing.

For VARIANT entries (if present):
   - Parse the parent ENSP and translate to gene symbol.
   - Report as "sequence-variant proteoform of {GENE} carrying {INDEL/SNV}
     at codon {N}."
   - If both canonical and variant are present, comparing their
     fold-changes is the killer plot — a variant that behaves
     differently from its canonical parent suggests altered function.
   - For tumor samples, variant proteoforms are candidate neoantigens.

GROUPING IN REPORT STRUCTURE:
   Always discuss proteogenomic findings (NOVEL_GENE, NOVEL_ISOFORM,
   VARIANT) in a SEPARATE section from canonical findings. Do not mix
   them into the same "top proteins" lists.

PEPTIDE-LEVEL EVIDENCE:
   For NOVEL_GENE and NOVEL_ISOFORM hits, the supporting peptide MUST
   contain sequence that distinguishes the novel form from any canonical
   protein. DE-LIMP enforces this filter automatically; if you see a
   novel-class entry, trust that it has been peptide-of-discovery
   validated.

----------------------------------------------------------------------------
GENERAL CAUTIONS
----------------------------------------------------------------------------

- The expanded database is roughly 3-8× larger than UniProt alone.
  Multiple-testing burden is higher; some non-canonical hits may be
  marginal even at FDR < 0.05.
- Peptide count matters MORE for proteogenomic hits than canonical:
  a 1-peptide canonical hit is routine; a 1-peptide novel-class hit
  is provisional and should be flagged as such.
- When generating biological interpretation, group non-canonical hits
  separately from canonical hits.
- Do not include MSTRG.* or INDEL_ENSP* identifiers in main-text
  result claims of a manuscript — these belong in supplementary
  tables with explicit "candidate" or "putative" framing.

============================================================================
End of glossary. For questions about specific entries, consult the
Proteogenomics tab in DE-LIMP.
```

R helper to ship the glossary:

```r
write_proteogenomics_glossary <- function(zip_dir, app_dir = NULL) {
  glossary_path <- file.path(zip_dir, "Proteogenomics_Glossary.txt")
  base <- if (!is.null(app_dir)) app_dir else getwd()
  src <- file.path(base, "scripts", "proteogenomics_glossary.txt")
  if (!file.exists(src)) {
    warning("Proteogenomics glossary not found at ", src,
            " — skipping glossary in ZIP")
    return(invisible(NULL))
  }
  file.copy(src, glossary_path, overwrite = TRUE)
  invisible(glossary_path)
}
```

DE-LIMP is run via `shiny::runApp()` and is NOT an installed R package, so
`system.file("extdata", ..., package = "delimp")` does not resolve. The
glossary lives at `scripts/proteogenomics_glossary.txt` at the repo root,
resolved relative to `app_dir`. This mirrors the pattern of
`get_contaminant_fasta(library_name, app_dir = NULL)` in `R/helpers_search.R`.
Update only when classification rules or pipeline output format changes.

---

## Conditional prompt blocks

These slot into all three templates wherever `{phospho_section_*}` and
`{gsea_section_*}` slot. Same conditional gating pattern as those.

### Detection note (top of every template, after `{gsea_note}`)

```r
proteog_note <- if (isTRUE(values$is_proteogenomics)) {
  pc <- values$protein_classification
  n_ref     <- sum(pc$source == "REF")
  n_novel_g <- sum(pc$source == "NOVEL_GENE")
  n_novel_i <- sum(pc$source == "NOVEL_ISOFORM")
  n_variant <- sum(pc$source == "VARIANT")
  n_uniprot <- sum(pc$source == "UNIPROT")
  
  parts <- c(
    sprintf("%s canonical UniProt", format(n_uniprot, big.mark = ",")),
    sprintf("%s reference (Ensembl)", format(n_ref, big.mark = ","))
  )
  if (n_novel_g > 0) parts <- c(parts, sprintf("%s novel genes", format(n_novel_g, big.mark = ",")))
  if (n_novel_i > 0) parts <- c(parts, sprintf("%s novel isoforms", format(n_novel_i, big.mark = ",")))
  if (n_variant > 0) parts <- c(parts, sprintf("%s variant proteoforms", format(n_variant, big.mark = ",")))
  
  paste0(
    "Database type: PROTEOGENOMICS-EXPANDED (",
    paste(parts, collapse = " + "),
    ")\n",
    "→ See Proteogenomics_Glossary.txt for identifier decoding.\n"
  )
} else ""
```

Inject `{proteog_note}` into EXPERIMENT OVERVIEW after `{gsea_note}`.

### File-list line (top of every template)

```r
proteog_file_note <- if (isTRUE(values$is_proteogenomics)) {
  "- Proteogenomics_Glossary.txt — How to decode REF/NOVEL_GENE/NOVEL_ISOFORM identifiers\n"
} else ""
```

Insert `{proteog_file_note}` after `{phospho_file_note}` in each template's file list.

### Brief template section

```r
proteog_section_brief <- if (isTRUE(values$is_proteogenomics)) {
  paste0(
    "## 8. Proteogenomic Discoveries\n",
    "This experiment used an expanded search database — see ",
    "Proteogenomics_Glossary.txt for full guidance.\n\n",
    "Generate a SEPARATE, clearly-labeled subsection for non-canonical hits. ",
    "Do NOT mix these into the canonical \"Key Findings\" section above.\n\n",
    "**Identify proteogenomic hits**: in DE_Results_Full.csv, look at the ",
    "`source` column (or parse from Protein.Group.Description). Non-canonical ",
    "classes are NOVEL_GENE, NOVEL_ISOFORM, and (if present) VARIANT.\n\n",
    "**For the Brief, produce only**:\n",
    "- A 2-3 sentence introduction explaining that this experiment combined ",
    "standard proteomics with sample-specific RNA-seq-derived ORFs.\n",
    "- A compact table of up to 10 significant non-canonical hits, ranked by ",
    "ORF_type confidence (`complete` first), then by |log2FC|:\n",
    "  | Identifier | Class | ORF_type | log2FC | FDR | Notes |\n",
    "- The Notes column should briefly describe each: e.g. ",
    "\"complete ORF, novel intergenic locus on chromosome from parent_gene MSTRG.10029\" ",
    "or \"novel isoform of Trim25, 5'-truncated\".\n",
    "- One closing sentence: \"These hits are candidate novel/variant proteins ",
    "and warrant orthogonal validation before biological interpretation.\"\n\n",
    "**Do NOT**:\n",
    "- speculate on biological function of NOVEL_GENE identifiers\n",
    "- claim a non-canonical hit is a known protein\n",
    "- include non-canonical hits in the volcano plot top-label list — label ",
    "only canonical (REF/UNIPROT) proteins on the volcano. Mention the ",
    "non-canonical count in the volcano caption as: \"An additional N non-canonical ",
    "(proteogenomic) hits are tabulated separately.\"\n"
  )
} else ""
```

Slot:
```
{gsea_section_brief}
{phospho_section_brief}
{proteog_section_brief}            ← NEW
```

### Full template section

```r
proteog_section_full <- if (isTRUE(values$is_proteogenomics)) {
  paste0(
    "## 5b. Proteogenomic Discoveries\n\n",
    "This experiment used a proteogenomics-expanded database. See ",
    "Proteogenomics_Glossary.txt for identifier decoding rules; treat that ",
    "file as authoritative.\n\n",
    
    "### 5b-i. Database Composition Summary\n",
    "Report the breakdown of identified proteins by source class. Use the ",
    "counts from the experiment overview as a starting point, and note how ",
    "many of each class were significantly regulated in the contrasts.\n\n",
    
    "### 5b-ii. Significant Novel Genes (NOVEL_GENE)\n",
    "For each contrast, list all NOVEL_GENE hits passing significance. Format:\n",
    "| Identifier | ORF_type | log2FC | FDR | Peptides | Notes |\n\n",
    "Sort within the table by ORF_type confidence: ",
    "`complete` first, then `5prime_partial`, then `3prime_partial`, then ",
    "`internal`. Within each tier, sort by |log2FC|.\n\n",
    "For each hit, the Notes column should briefly describe the ORF ",
    "(e.g., \"novel intergenic locus, 285 aa, on minus strand of chromosome [parse from parent_gene if possible]\") ",
    "but do NOT invent a function.\n\n",
    
    "### 5b-iii. Significant Novel Isoforms (NOVEL_ISOFORM)\n",
    "For each contrast, list all NOVEL_ISOFORM hits passing significance. ",
    "Critically, for each one CHECK whether the canonical reference parent ",
    "protein is ALSO present in DE_Results_Full.csv. If both are present:\n",
    "- Generate a small comparison: canonical fold-change vs novel-isoform fold-change\n",
    "- Note whether they move concordantly (likely co-regulated) or differentially ",
    "  (novel isoform may have altered stability/regulation)\n",
    "- Differential behavior is the most biologically interesting outcome\n\n",
    "Table format:\n",
    "| Novel ID | Parent Gene Symbol | log2FC (novel) | log2FC (canonical, if detected) | FDR | Interpretation |\n\n",
    
    "### 5b-iv. Variant Proteoforms (if present)\n",
    "Only generate this subsection if any VARIANT-class entries are in the ",
    "results. For each:\n",
    "- Parse parent ENSP, translate to gene symbol\n",
    "- Report codon position and ref→alt change\n",
    "- If canonical parent protein is ALSO present, compare fold-changes\n",
    "- For tumor samples, mention candidate neoantigen status with appropriate ",
    "  caution (requires explicit clinical context to claim)\n\n",
    
    "### 5b-v. Interpretation Caveats\n",
    "Write 1 paragraph noting:\n",
    "- Expanded database increases multiple-testing burden\n",
    "- Single-peptide proteogenomic hits should be flagged as provisional\n",
    "- Orthogonal validation (targeted MS, de novo sequencing, matched RNA-seq) ",
    "is recommended for any biologically significant finding\n",
    "- For NOVEL_GENE entries specifically, the lack of prior annotation is ",
    "by definition — these are discoveries that warrant follow-up, not ",
    "established proteins\n\n",
    
    "### 5b-vi. Figure — Class-Stratified Volcano\n",
    "Generate a volcano plot where points are colored by `source` class:\n",
    "- canonical (REF/UNIPROT): grey, smaller points\n",
    "- NOVEL_GENE: orange, slightly larger points\n",
    "- NOVEL_ISOFORM: purple, slightly larger points\n",
    "- VARIANT (if present): teal, slightly larger points\n\n",
    "Use the same axes and significance lines as the main volcano. Label only ",
    "the top 5 significant hits PER CLASS (not just by overall |logFC|), so that ",
    "non-canonical discoveries are visible even when their fold-changes are ",
    "moderate. Use ORF_type as a label-decoration if space allows (e.g., ",
    "asterisk for `complete`, dash for partial).\n\n",
    "Caption: \"Figure X. Class-stratified volcano plot. Points colored by ",
    "protein source: canonical reference (grey, n=...), StringTie novel genes ",
    "(orange, n=...), novel isoforms of annotated genes (purple, n=...). ",
    "Top 5 per non-canonical class labeled, with ORF_type indicated.\"\n"
  )
} else ""
```

Slot:
```
## 5. Phosphoproteomics Results
{phospho_section_full}

{proteog_section_full}            ← NEW

## 6. Biological Synthesis
```

Also modify the Biological Synthesis prompt (Section 6) with a conditional append:

```r
biosynth_proteog_note <- if (isTRUE(values$is_proteogenomics)) {
  paste0(
    "\n- Discuss canonical and non-canonical findings SEPARATELY. ",
    "When suggesting follow-up experiments, distinguish:\n",
    "  (a) experiments targeting well-annotated proteins (functional validation)\n",
    "  (b) experiments targeting NOVEL_GENE candidates (identity validation — ",
    "      e.g., targeted MS for the specific peptide, matched RNA-seq, ",
    "      antibody generation if peptide is long enough)\n",
    "  (c) experiments targeting NOVEL_ISOFORM candidates (functional ",
    "      comparison to canonical isoform — alternative localization, ",
    "      stability, binding partner studies)\n",
    "- If both canonical and novel-isoform forms of the same gene appear in ",
    "  results, prioritize this finding in the discussion — it's the most ",
    "  interpretable proteogenomic signal."
  )
} else ""
```

Inject before the closing `\n` of section 6.

### Manuscript template section

```r
proteog_section_manuscript <- if (isTRUE(values$is_proteogenomics)) {
  paste0(
    "### Proteogenomic Identifications\n\n",
    "[Write a Results paragraph reporting non-canonical findings. Use the ",
    "following structure verbatim, filling in numbers from DE_Results_Full.csv:]\n\n",
    
    "\"To detect sample-specific proteins absent from standard reference ",
    "proteomes, mass spectra were searched against a proteogenomics-expanded ",
    "database. The database combined the canonical UniProt reference proteome ",
    "with [N] StringTie-predicted open reading frames derived from sample-",
    "matched RNA-seq, comprising [N] putative novel intergenic ORFs and [N] ",
    "novel splicing isoforms of annotated genes. At an FDR threshold of ",
    "{q_cutoff}, we identified [N] significantly regulated proteogenomic hits ",
    "(Supplementary Table SX): [N] novel-gene candidates and [N] novel-isoform ",
    "candidates. Of the novel-gene candidates, [N] were structurally complete ",
    "ORFs (containing both start and stop codons). [If VARIANT entries: ",
    "Additionally, [N] sequence-variant proteoforms of canonical proteins ",
    "were identified.] These proteogenomic identifications are provisional ",
    "and would benefit from orthogonal validation by targeted mass spectrometry ",
    "or matched transcript-level evidence.\"\n\n",
    
    "[Write a Methods paragraph:]\n\n",
    
    "\"For proteogenomics analysis, the spectral search database was constructed ",
    "by concatenating [reference proteome name and version] with predicted ORFs ",
    "from StringTie (Pertea et al., 2015) assembly of matched RNA-seq data. ",
    "RNA-seq reads were quality-trimmed with fastp (Chen et al., 2018), filtered ",
    "against organism-specific ribosomal RNA with bowtie2 (Langmead and Salzberg, ",
    "2012), aligned to the reference genome with STAR (Dobin et al., 2013), ",
    "assembled with StringTie into per-sample then merged transcript models, ",
    "and classified against the reference annotation using gffcompare ",
    "(Pertea and Pertea, 2020) to distinguish novel intergenic loci from novel ",
    "isoforms of annotated genes. Open reading frames were predicted with ",
    "TransDecoder, retaining the single best ORF per transcript. Differential ",
    "expression statistics for proteogenomic identifications were computed using ",
    "the same limpa/limma framework as for canonical proteins; FDR thresholds ",
    "for novel-ORF and isoform hits should be interpreted in the context of the ",
    "expanded search space.\"\n\n",
    
    "**Critical formatting rules**:\n",
    "- Do NOT include MSTRG.* identifiers in main manuscript text body\n",
    "- Do NOT include MSTRG.* identifiers in top-protein lists in Results figures\n",
    "- Do NOT include MSTRG.* identifiers in main volcano plot labels\n",
    "- All non-canonical hits belong in SUPPLEMENTARY tables and figures\n",
    "- Reference them as \"candidate\" identifications\n",
    "- Do NOT propose biological functions for NOVEL_GENE candidates\n",
    "- For NOVEL_ISOFORM candidates, you MAY reference the parent gene's known ",
    "function but ONLY in the context of \"the canonical form is known to ... ; ",
    "this novel isoform may have altered ... pending functional characterization\"\n",
    "- Do NOT use any non-canonical identifier in a sentence claiming biological causation\n"
  )
} else ""
```

Slot:
```
### Pathway Enrichment
{gsea_section_manuscript}

### Phosphoproteomics
{phospho_section_manuscript}

{proteog_section_manuscript}            ← NEW

## Methods
```

---

## Inline data summaries

Append to INLINE DATA SUMMARIES block of all three templates:

```r
proteog_inline <- if (isTRUE(values$is_proteogenomics)) {
  pc <- values$protein_classification
  noncanon_ids <- pc$Protein.Group[pc$source %in% c("NOVEL_GENE", "NOVEL_ISOFORM", "VARIANT")]
  
  contrasts <- colnames(values$fit$contrasts)
  blocks <- lapply(contrasts, function(ct) {
    tt <- limma::topTable(values$fit, coef = ct, number = Inf)
    tt$Protein.Group <- rownames(tt)
    
    sig_nc <- tt %>%
      dplyr::filter(Protein.Group %in% noncanon_ids,
                    adj.P.Val < input$q_cutoff,
                    abs(logFC) > input$logfc_cutoff) %>%
      dplyr::left_join(pc, by = "Protein.Group") %>%
      # Sort by class, then ORF_type confidence, then |logFC|
      dplyr::mutate(
        confidence_rank = dplyr::case_when(
          orf_type == "complete"        ~ 1,
          orf_type == "5prime_partial"  ~ 2,
          orf_type == "3prime_partial"  ~ 3,
          orf_type == "internal"        ~ 4,
          TRUE                          ~ 5
        )
      ) %>%
      dplyr::arrange(source, confidence_rank, dplyr::desc(abs(logFC))) %>%
      utils::head(15) %>%
      dplyr::transmute(
        Identifier  = Protein.Group,
        Source      = source,
        ORF_type    = orf_type,
        ParentGene  = parent_gene,
        logFC       = round(logFC, 2),
        FDR         = signif(adj.P.Val, 2)
      )
    
    if (nrow(sig_nc) == 0) {
      sprintf("Contrast %s: no significant non-canonical hits.", ct)
    } else {
      paste0(
        sprintf("Contrast %s — top non-canonical significant hits (sorted by class, then ORF_type confidence):\n", ct),
        paste(utils::capture.output(print(sig_nc, row.names = FALSE)),
              collapse = "\n")
      )
    }
  })
  
  paste0(
    "\n--- PROTEOGENOMIC DISCOVERIES (non-canonical significant hits) ---\n",
    paste(unlist(blocks), collapse = "\n\n"),
    "\n\n(See Proteogenomics_Glossary.txt for identifier decoding.)"
  )
} else ""
```

Inject `{proteog_inline}` at end of each template's INLINE DATA SUMMARIES block,
after `{phospho_inline}`.

---

## `build_claude_prompt()` changes

```r
build_claude_prompt <- function(type, ctx, values, input, meta, export_date) {
  # ... existing instrument_block, phospho_*, gsea_* assembly ...
  
  # Proteogenomics conditionals
  proteog_note               <- build_proteog_note(values)
  proteog_file_note          <- build_proteog_file_note(values)
  proteog_section_brief      <- build_proteog_section(values, "brief")
  proteog_section_full       <- build_proteog_section(values, "full")
  proteog_section_manuscript <- build_proteog_section(values, "manuscript")
  biosynth_proteog_note      <- build_biosynth_proteog_note(values)
  proteog_inline             <- build_proteog_inline(values, input)
  
  # ... existing glue::glue() template assembly with new placeholders added ...
}
```

The seven `build_proteog_*` helpers live in `R/helpers_proteogenomics.R`. Each
returns `""` when `!isTRUE(values$is_proteogenomics)`. Strictly additive.

---

## ZIP assembly change

```r
# Existing files
write_de_results_csv(...)
write_expression_matrix_csv(...)
write_qc_metrics_csv(...)
# ...

if (isTRUE(values$is_proteogenomics)) {
  write_proteogenomics_glossary(zip_dir)
}
```

---

## Session save/load

```r
# save
session_data$is_proteogenomics      <- values$is_proteogenomics
session_data$protein_classification <- values$protein_classification

# load
values$is_proteogenomics      <- session_data$is_proteogenomics      %||% FALSE
values$protein_classification <- session_data$protein_classification %||% NULL

# back-compat: if older session has no classification, recompute on load
if (!is.null(values$diann_report) && is.null(values$protein_classification)) {
  values$protein_classification <- classify_proteins(values$diann_report)
  values$is_proteogenomics <- any(
    values$protein_classification$source %in% c("REF", "NOVEL_GENE", "NOVEL_ISOFORM", "VARIANT")
  )
}
```

---

## Testing checklist

### Detection
- [ ] Standard UniProt search → `is_proteogenomics = FALSE`, no proteogenomic sections
- [ ] Pipeline-produced FASTA → `is_proteogenomics = TRUE`
- [ ] `classify_proteins()` correctly extracts source/orf_type/parent_gene from descriptions
- [ ] Old-style UniProt IDs without source= tags default to UNIPROT class
- [ ] VARIANT detection by INDEL_ENSP*/SNV_ENSP* prefix works

### Glossary file
- [ ] `Proteogenomics_Glossary.txt` appears in ZIP when proteogenomic, absent otherwise
- [ ] Glossary <12 KB
- [ ] Renders cleanly in plain-text viewers

### Prompt injection — Brief
- [ ] `{proteog_note}` shows correct counts across all classes
- [ ] Section 8 appears with table of top 10 ranked by confidence then |logFC|
- [ ] Volcano caption mentions non-canonical count separately

### Prompt injection — Full
- [ ] Section 5b appears with all six subsections (i–vi)
- [ ] Class-stratified volcano specification present
- [ ] Section 6 includes follow-up experiment differentiation by class
- [ ] NOVEL_ISOFORM canonical-vs-novel comparison requested

### Prompt injection — Manuscript
- [ ] `{proteog_section_manuscript}` appears in Results
- [ ] Methods paragraph correctly references STAR, bowtie2, gffcompare in addition to StringTie/TransDecoder
- [ ] Critical-warning block forbidding MSTRG.* in main text is prominent

### Inline data
- [ ] Up to 15 sig non-canonical hits per contrast, sorted correctly
- [ ] Sort order is: source class → ORF_type confidence → |logFC|
- [ ] Empty case handled
- [ ] Parent gene column populated where available

### Session round-trip
- [ ] Save/reload preserves `is_proteogenomics` and `protein_classification`
- [ ] Pre-feature session reloads cleanly, recomputes classification

### End-to-end
- [ ] Real pipeline-produced FASTA → Brief → Claude does NOT invent functions for MSTRG.* (manual review)
- [ ] Full → class-stratified volcano produced
- [ ] Manuscript → MSTRG.* identifiers appear ONLY in supplementary references

---

## Files modified

| File | Change |
|------|--------|
| `R/helpers_proteogenomics.R` | **New** — `classify_proteins()` reads source= tags; all `build_proteog_*` helpers |
| `scripts/proteogenomics_glossary.txt` | **New** — Glossary text shown above. Resolved relative to `app_dir` at runtime (DE-LIMP is not an installed R package; `inst/extdata` conventions do not apply). |
| `R/server_ai.R` | Add proteog conditionals to `build_claude_prompt()`; glossary in ZIP |
| `R/server_session.R` | Save/load; back-compat recompute |
| `app.R` | Add `is_proteogenomics`, `protein_classification` to `reactiveValues()`. Call `classify_proteins()` from the existing DIA-NN-load handler in `R/server_data.R` or `R/server_search.R` (search the codebase for the post-`readDIANN()` block — that's the insertion point). |

---

## Why v1.1 (vs. expanding v1.0)

- **Strictly additive**: standard UniProt sessions unaffected
- **Reads explicit metadata, not prefix patterns**: source= tag is reliable;
  prefix matching breaks on edge cases (e.g., Ensembl protein IDs without `.p1` suffix)
- **ORF_type confidence tiering**: a `complete` ORF deserves different treatment
  than an `internal` one; v1.1 surfaces this distinction
- **NOVEL_ISOFORM as its own class**: v1.0 didn't have this because the pipeline
  plan didn't include gffcompare; the validation showed gffcompare is essential
- **Manuscript safety preserved**: explicit guardrails against hallucinating
  function claims for novel-gene candidates in publication-bound text
- **Same architectural pattern as Spectronaut/DIA-NN log addendums**: implementation
  reads naturally for someone familiar with the codebase

---

*Addendum version 1.1 — Brett Phinney / UC Davis Proteomics Core — May 2026 (post-validation)*
