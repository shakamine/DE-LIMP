# helpers_proteog_qc.R — Log parsers and quality-gate enforcement for the
# proteogenomics pipeline. Pure functions, no Shiny reactivity.
#
# Two responsibility groups:
#   1) parse_*() — extract metrics from STAR / fastp / bowtie2 logs
#   2) check_*() — apply quality gates and return structured pass/fail
#
# The R-side check_*() helpers are used by:
#   - The orchestrator to pre-validate logs after each stage
#   - The UI status panel to render gate diagnostics
#   - Tests, against captured logs from validation
#
# The SLURM-level gate (qc_gate sbatch from helpers_rnaseq.R) enforces the same
# rules at compute-time via exit code; these R helpers re-check at submission
# time and at result-display time for richer diagnostics.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# =============================================================================
# Log parsers
# =============================================================================

#' Parse STAR's Log.final.out into a list of numeric metrics
#'
#' @param log_path character — path to a STAR `*_Log.final.out`
#' @return list with $unique_pct, $multi_pct, $toomany_pct, $unmapped_pct,
#'   $n_input_reads, $avg_input_read_length, $avg_mapped_length,
#'   $mismatch_rate_pct, $n_splices_total, $n_splices_annotated, $pct_annotated_splices,
#'   or NULL if the file cannot be parsed.
parse_star_log <- function(log_path) {
  if (!file.exists(log_path)) {
    warning("parse_star_log(): not found: ", log_path)
    return(NULL)
  }
  lines <- readLines(log_path, warn = FALSE)
  if (length(lines) == 0) return(NULL)

  # Helper: pull the "value" after the `|` separator on a STAR log line.
  # STAR's format is "                key |\tvalue", so split on "|" and trim.
  pull <- function(needle, as_pct = FALSE) {
    hit <- grep(needle, lines, fixed = TRUE, value = TRUE)
    if (length(hit) == 0) return(NA_real_)
    raw <- sub(".*\\|", "", hit[1])
    raw <- trimws(raw)
    if (as_pct) raw <- sub("%$", "", raw)
    suppressWarnings(as.numeric(raw))
  }

  unique_pct           <- pull("Uniquely mapped reads %", as_pct = TRUE)
  multi_pct            <- pull("% of reads mapped to multiple loci", as_pct = TRUE)
  toomany_pct          <- pull("% of reads mapped to too many loci", as_pct = TRUE)
  unmapped_short_pct   <- pull("% of reads unmapped: too short", as_pct = TRUE)
  unmapped_other_pct   <- pull("% of reads unmapped: other", as_pct = TRUE)
  unmapped_mm_pct      <- pull("% of reads unmapped: too many mismatches", as_pct = TRUE)
  n_input              <- pull("Number of input reads")
  avg_input_len        <- pull("Average input read length")
  avg_mapped_len       <- pull("Average mapped length")
  mismatch_rate        <- pull("Mismatch rate per base, %", as_pct = TRUE)
  n_splices_total      <- pull("Number of splices: Total")
  n_splices_annotated  <- pull("Number of splices: Annotated (sjdb)")

  unmapped_pct <- sum(
    c(unmapped_short_pct, unmapped_other_pct, unmapped_mm_pct),
    na.rm = TRUE
  )
  pct_annotated_splices <- if (!is.na(n_splices_total) && n_splices_total > 0) {
    100 * (n_splices_annotated %||% 0) / n_splices_total
  } else NA_real_

  list(
    unique_pct             = unique_pct,
    multi_pct              = multi_pct,
    toomany_pct            = toomany_pct,
    unmapped_pct           = unmapped_pct,
    n_input_reads          = n_input,
    avg_input_read_length  = avg_input_len,
    avg_mapped_length      = avg_mapped_len,
    mismatch_rate_pct      = mismatch_rate,
    n_splices_total        = n_splices_total,
    n_splices_annotated    = n_splices_annotated,
    pct_annotated_splices  = pct_annotated_splices
  )
}

#' Parse fastp's JSON report for read-length and pass-rate metrics
#'
#' @param json_path character
#' @return list with $median_read_length_pre, $median_read_length_post,
#'   $total_reads_pre, $total_reads_post, $pct_passed, or NULL.
parse_fastp_json <- function(json_path) {
  if (!file.exists(json_path)) {
    warning("parse_fastp_json(): not found: ", json_path)
    return(NULL)
  }
  raw <- tryCatch(jsonlite::read_json(json_path), error = function(e) NULL)
  if (is.null(raw)) return(NULL)

  total_pre  <- raw$summary$before_filtering$total_reads  %||% NA_real_
  total_post <- raw$summary$after_filtering$total_reads   %||% NA_real_
  median_pre <- raw$summary$before_filtering$read1_mean_length %||% NA_real_
  median_post <- raw$summary$after_filtering$read1_mean_length %||% NA_real_

  pct_passed <- if (is.finite(total_pre) && total_pre > 0) {
    100 * total_post / total_pre
  } else NA_real_

  list(
    median_read_length_pre  = median_pre,
    median_read_length_post = median_post,
    total_reads_pre         = total_pre,
    total_reads_post        = total_post,
    pct_passed              = pct_passed
  )
}

#' Parse bowtie2's rRNA-filter log for overall alignment rate
#'
#' bowtie2 writes a textual summary like "10.57% overall alignment rate".
#'
#' @param log_path character
#' @return list with $overall_alignment_pct, $concordant_pairs_aligned,
#'   $n_input_pairs, or NULL.
parse_rrna_log <- function(log_path) {
  if (!file.exists(log_path)) {
    warning("parse_rrna_log(): not found: ", log_path)
    return(NULL)
  }
  lines <- readLines(log_path, warn = FALSE)
  if (length(lines) == 0) return(NULL)

  overall_pct <- NA_real_
  hit <- grep("overall alignment rate", lines, value = TRUE)
  if (length(hit) > 0) {
    overall_pct <- suppressWarnings(as.numeric(
      sub("([0-9.]+)%.*", "\\1", hit[1])
    ))
  }

  n_input <- NA_real_
  hit2 <- grep("reads; of these:", lines, value = TRUE)
  if (length(hit2) > 0) {
    n_input <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", hit2[1])))
  }

  concordant <- NA_real_
  hit3 <- grep("aligned concordantly >1 times", lines, value = TRUE)
  if (length(hit3) > 0) {
    concordant <- suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*", "\\1", hit3[1])))
  }

  list(
    overall_alignment_pct    = overall_pct,
    n_input_pairs            = n_input,
    concordant_pairs_aligned = concordant
  )
}

# =============================================================================
# Quality gates
# =============================================================================

# Likely-causes catalog for low-mapping diagnostics. Surfaced verbatim to the
# user when a gate fails (CLAUDE.md rule #4 — no silent halts).
.LOW_MAPPING_CAUSES <- c(
  "Wrong reference genome selected (e.g., human reference applied to mouse sample)",
  "Heavy non-rRNA contamination (mitochondrial overload, host cell line, bacterial)",
  "Library type unsuited to this pipeline (Ribo-Seq, CLIP-Seq, 3'-Tag-Seq)",
  "Severe sample degradation (RIN < 5)"
)

#' Apply the uniquely-mapped quality gate to a parsed STAR log
#'
#' @param star_log_path character — path to STAR Log.final.out (or pre-parsed list)
#' @param tier_params list — from select_star_params(); needs $qc_gate_unique_pct
#' @return list with $pass, $unique_pct, $gate, $tier, and on failure $message + $causes
check_alignment_quality <- function(star_log_path, tier_params) {
  parsed <- if (is.list(star_log_path) && "unique_pct" %in% names(star_log_path)) {
    star_log_path
  } else {
    parse_star_log(star_log_path)
  }

  if (is.null(parsed) || !is.finite(parsed$unique_pct %||% NA_real_)) {
    return(list(
      pass       = FALSE,
      unique_pct = NA_real_,
      gate       = tier_params$qc_gate_unique_pct,
      tier       = tier_params$tier,
      message    = "STAR log unreadable or missing 'Uniquely mapped reads %' field — cannot evaluate gate.",
      causes     = .LOW_MAPPING_CAUSES
    ))
  }

  gate <- tier_params$qc_gate_unique_pct
  if (parsed$unique_pct < gate) {
    list(
      pass       = FALSE,
      unique_pct = parsed$unique_pct,
      gate       = gate,
      tier       = tier_params$tier,
      message    = sprintf(
        "Alignment quality below threshold for tier '%s'. Uniquely mapped: %.1f%%, required: %d%%. Pipeline halted.",
        tier_params$tier, parsed$unique_pct, gate
      ),
      causes     = .LOW_MAPPING_CAUSES
    )
  } else {
    list(
      pass       = TRUE,
      unique_pct = parsed$unique_pct,
      gate       = gate,
      tier       = tier_params$tier
    )
  }
}

#' Aggregate per-sample STAR gate results
#'
#' @param star_logs character vector of paths
#' @param tier_params list — from select_star_params()
#' @return list with $pass (all samples passed?), $sample_results (per-sample list)
check_pipeline_gates <- function(star_logs, tier_params) {
  if (length(star_logs) == 0) {
    return(list(pass = FALSE, sample_results = list(),
                message = "No STAR logs provided"))
  }
  results <- lapply(star_logs, check_alignment_quality, tier_params = tier_params)
  names(results) <- vapply(star_logs, basename, character(1))
  all_pass <- all(vapply(results, function(r) isTRUE(r$pass), logical(1)))
  list(
    pass           = all_pass,
    sample_results = results
  )
}

#' Render a human-readable diagnostic block for a failing gate
#'
#' @param gate_result list — from check_alignment_quality()
#' @return character — multi-line text suitable for showing in a modal or log
render_gate_failure <- function(gate_result) {
  if (isTRUE(gate_result$pass)) return("Gate passed.")
  parts <- c(
    sprintf("QC GATE FAILED — alignment quality below threshold"),
    sprintf("  Tier:                  %s", gate_result$tier %||% "unknown"),
    sprintf("  Uniquely mapped:       %.1f%%",
            gate_result$unique_pct %||% NA_real_),
    sprintf("  Required:              %d%%",
            gate_result$gate %||% NA_integer_),
    "",
    if (!is.null(gate_result$message)) gate_result$message else NULL,
    "",
    "Likely causes (in order of frequency):"
  )
  for (i in seq_along(gate_result$causes %||% character())) {
    parts <- c(parts, sprintf("  %d. %s", i, gate_result$causes[i]))
  }
  parts <- c(parts, "",
             "DE-LIMP does NOT auto-fix this — please review the diagnostic and re-submit.")
  paste(parts, collapse = "\n")
}
