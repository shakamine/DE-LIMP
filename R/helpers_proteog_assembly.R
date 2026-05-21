# helpers_proteog_assembly.R — FASTA composition + registry I/O.
# No Shiny reactivity. Pure functions.
#
# Three responsibilities:
#  1) count_proteog_classes() — parse `source=` tags from FASTA headers and
#     count entries by class. Backward-compatible with v0.1 rewriter output
#     (which has no NOVEL_ISOFORM entries).
#  2) assemble_proteogenomics_fasta() — concatenate predicted ORFs + UniProt
#     reference (+ optional VARIANT + contaminants), dedup by sequence, and
#     register in the proteogenomics registry.
#  3) register_proteogenomics_fasta() / load_proteog_registry() — JSON I/O.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# Single source of truth for the proteogenomics-database registry path.
.proteog_registry_path <- function() {
  Sys.getenv(
    "DELIMP_PROTEOG_REGISTRY",
    unset = "/quobyte/proteomics-grp/de-limp/databases/proteogenomics/registry.json"
  )
}

# =============================================================================
# Composition counting
# =============================================================================

#' Count proteins in a FASTA by `source=` class
#'
#' Reads only header lines (`^>`). For each header:
#'   - if it contains `source=...` → that class
#'   - else if accession matches `^(INDEL|SNV)_ENSP` → VARIANT
#'   - else → UNIPROT (canonical reference, no proteogenomic tag)
#'
#' Always returns a list with EVERY class field present (zero-filled), so
#' downstream consumers don't have to handle missing keys.
#'
#' @param fasta_path character — path to a FASTA file (plain or gzipped)
#' @return list with: total, UNIPROT, REF, NOVEL_GENE, NOVEL_ISOFORM, VARIANT, UNPARSED
count_proteog_classes <- function(fasta_path) {
  if (!file.exists(fasta_path)) {
    stop("count_proteog_classes(): file not found: ", fasta_path)
  }

  con <- if (grepl("\\.gz$", fasta_path, ignore.case = TRUE)) {
    gzfile(fasta_path, "rt")
  } else {
    file(fasta_path, "rt")
  }
  on.exit(close(con), add = TRUE)

  out <- list(
    total         = 0L,
    UNIPROT       = 0L,
    REF           = 0L,
    NOVEL_GENE    = 0L,
    NOVEL_ISOFORM = 0L,
    VARIANT       = 0L,
    UNPARSED      = 0L
  )

  # Stream line-by-line — FASTAs can be large (hundreds of MB)
  while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
    if (!startsWith(line, ">")) next
    out$total <- out$total + 1L

    if (grepl("source=", line, fixed = TRUE)) {
      cls <- sub(".*source=([A-Z_]+).*", "\\1", line)
      if (cls %in% names(out)) {
        out[[cls]] <- out[[cls]] + 1L
      } else {
        out$UNPARSED <- out$UNPARSED + 1L
      }
    } else {
      # No source= tag — could be canonical UniProt or a VARIANT proteoform
      # (whose accession prefix marks it without needing a description tag).
      acc <- sub("^>([^[:space:]]+).*", "\\1", line)
      # Strip a leading sp|/tr| if present and grab the accession field.
      stripped <- sub("^(sp|tr)\\|([^|]+)\\|.*", "\\2", acc)
      if (grepl("^(INDEL|SNV)_ENSP", acc) || grepl("^(INDEL|SNV)_ENSP", stripped)) {
        out$VARIANT <- out$VARIANT + 1L
      } else {
        out$UNIPROT <- out$UNIPROT + 1L
      }
    }
  }
  out
}

# =============================================================================
# Registry I/O
# =============================================================================

#' Load the proteogenomics-database registry as a named list
#'
#' Empty/missing registry returns `list()`. Never throws — callers handle
#' an empty registry gracefully (no databases registered yet).
load_proteog_registry <- function() {
  path <- .proteog_registry_path()
  if (!file.exists(path)) return(list())
  raw <- tryCatch(jsonlite::read_json(path), error = function(e) NULL)
  if (is.null(raw)) {
    warning("load_proteog_registry(): could not parse ", path, " — returning empty registry")
    return(list())
  }
  if (length(raw) == 0) list() else raw
}

#' Write the proteogenomics-database registry from a named list
#'
#' @param registry named list (keys = project_name, values = entry list)
.save_proteog_registry <- function(registry) {
  path <- .proteog_registry_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(registry, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}

#' Register a built proteogenomics FASTA in the registry
#'
#' @param path             character — FASTA path
#' @param merged_gtf_path  character — preserved StringTie merged GTF alongside
#' @param project_name     character — user-supplied project identifier
#' @param composition      list — from count_proteog_classes()
#' @param build_metadata   list — optional extra fields (organism, build, n_samples, …)
#' @return invisible(path)
register_proteogenomics_fasta <- function(path,
                                          merged_gtf_path,
                                          project_name,
                                          composition,
                                          build_metadata = list()) {
  registry <- load_proteog_registry()

  entry <- c(
    list(
      path             = path,
      merged_gtf_path  = merged_gtf_path,
      project_name     = project_name,
      composition      = composition,
      created          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      created_by       = Sys.getenv("USER", unset = NA_character_),
      pipeline_version = "1.1"
    ),
    build_metadata
  )

  registry[[project_name]] <- entry
  .save_proteog_registry(registry)
  invisible(path)
}

# =============================================================================
# Assembly
# =============================================================================

#' Build a proteogenomics search database from component FASTAs
#'
#' Concatenates predicted ORFs + UniProt reference (+ optional variants +
#' optional cRAP contaminants), optionally dedups by sequence with `seqkit
#' rmdup -s`, copies the merged GTF alongside for coordinate lookups, and
#' registers the result.
#'
#' Self-describing: returned object carries `pipeline_id`, `methods_paragraph`,
#' `composition` — downstream consumers must read these rather than hardcoding
#' a description of "what we did" (CLAUDE.md architectural rule #1).
#'
#' @param project_name        character
#' @param uniprot_fasta       character path
#' @param predicted_orfs_fasta character path (from header rewriter)
#' @param merged_gtf          character path — StringTie merged.gtf
#' @param variant_fasta       character path or NULL
#' @param contaminants_fasta  character path or NULL
#' @param output_dir          character — where the final FASTA + GTF land
#' @param dedupe              logical — run `seqkit rmdup -s`
#' @param build_metadata      list — optional metadata for registry (organism, etc.)
#' @return list with $path, $merged_gtf_path, $composition, $pipeline_id, $methods_paragraph
assemble_proteogenomics_fasta <- function(project_name,
                                          uniprot_fasta,
                                          predicted_orfs_fasta,
                                          merged_gtf,
                                          variant_fasta = NULL,
                                          contaminants_fasta = NULL,
                                          output_dir,
                                          dedupe = TRUE,
                                          build_metadata = list()) {
  # Validate inputs (boundary check — CLAUDE.md rule on input validation).
  required_files <- c(uniprot_fasta, predicted_orfs_fasta, merged_gtf)
  missing <- required_files[!file.exists(required_files)]
  if (length(missing) > 0) {
    stop("assemble_proteogenomics_fasta(): missing required input(s): ",
         paste(missing, collapse = ", "))
  }
  if (!nzchar(project_name) || !grepl("^[A-Za-z0-9_.-]+$", project_name)) {
    stop("assemble_proteogenomics_fasta(): project_name must be non-empty and contain only [A-Za-z0-9_.-]; got: ", project_name)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out_path <- file.path(
    output_dir,
    sprintf("%s_proteogenomics_%s.fasta",
            project_name,
            format(Sys.Date(), "%Y_%m"))
  )

  # Compose the component list in deterministic order: predicted ORFs first
  # so their headers come up early in the FASTA (useful for spot-checking),
  # then UniProt, then optional VARIANT, then optional contaminants.
  components <- c(predicted_orfs_fasta, uniprot_fasta)
  if (!is.null(variant_fasta) && file.exists(variant_fasta)) {
    components <- c(components, variant_fasta)
  }
  if (!is.null(contaminants_fasta) && file.exists(contaminants_fasta)) {
    components <- c(components, contaminants_fasta)
  }

  if (dedupe) {
    tmp_concat <- tempfile(pattern = "proteog_concat_", fileext = ".fasta")
    on.exit(if (file.exists(tmp_concat)) file.remove(tmp_concat), add = TRUE)

    cat_status <- system2(
      "cat", args = shQuote(components),
      stdout = tmp_concat, stderr = FALSE
    )
    if (cat_status != 0) {
      stop("assemble_proteogenomics_fasta(): cat failed concatenating components")
    }

    rmdup_status <- system2(
      "seqkit",
      args = c("rmdup", "-s", "-o", shQuote(out_path), shQuote(tmp_concat))
    )
    if (rmdup_status != 0) {
      stop("assemble_proteogenomics_fasta(): seqkit rmdup failed (exit ", rmdup_status, ")")
    }
  } else {
    cat_status <- system2(
      "cat", args = shQuote(components),
      stdout = out_path, stderr = FALSE
    )
    if (cat_status != 0) {
      stop("assemble_proteogenomics_fasta(): cat failed concatenating components (no dedupe)")
    }
  }

  composition <- count_proteog_classes(out_path)

  # Preserve merged GTF alongside (coordinate-lookup reference) — registry
  # records the path so the Proteogenomics tab can resolve genomic coords on demand.
  gtf_dest <- sub("\\.fasta$", "_merged.gtf", out_path)
  file.copy(merged_gtf, gtf_dest, overwrite = TRUE)

  register_proteogenomics_fasta(
    path             = out_path,
    merged_gtf_path  = gtf_dest,
    project_name     = project_name,
    composition      = composition,
    build_metadata   = build_metadata
  )

  list(
    pipeline_id       = "proteogenomics_v1.1",
    path              = out_path,
    merged_gtf_path   = gtf_dest,
    composition       = composition,
    methods_paragraph = .proteog_methods_paragraph(composition, build_metadata)
  )
}

# Self-describing methods text used by downstream Claude export / methods readme.
# CLAUDE.md rule #1: pipeline objects must self-describe; downstream consumers
# read $methods_paragraph rather than hardcoding what the pipeline did.
.proteog_methods_paragraph <- function(composition, meta = list()) {
  fmt <- function(n) format(n %||% 0L, big.mark = ",")
  paste0(
    "The spectral search database was constructed by concatenating ",
    (meta$uniprot_release %||% "the UniProt reference proteome"),
    " with StringTie-predicted ORFs from sample-matched RNA-seq (",
    fmt(composition$REF), " reference-derived predicted, ",
    fmt(composition$NOVEL_GENE), " novel-gene candidates, ",
    fmt(composition$NOVEL_ISOFORM), " novel-isoform candidates",
    if ((composition$VARIANT %||% 0L) > 0) {
      paste0(", ", fmt(composition$VARIANT), " variant proteoforms")
    } else "",
    "). Total entries after deduplication: ", fmt(composition$total),
    " (", fmt(composition$UNIPROT), " canonical UniProt + ",
    fmt(composition$REF + composition$NOVEL_GENE + composition$NOVEL_ISOFORM),
    " proteogenomic).",
    if (!is.null(meta$rrna_pct_mean)) {
      sprintf(" Mean rRNA filter rate: %.1f%%.", meta$rrna_pct_mean)
    } else "",
    if (!is.null(meta$uniquely_mapped_pct_mean)) {
      sprintf(" Mean uniquely-mapped rate: %.1f%%.", meta$uniquely_mapped_pct_mean)
    } else "",
    if (!is.null(meta$read_length_tier)) {
      sprintf(" STAR threshold tier: %s.", meta$read_length_tier)
    } else ""
  )
}
