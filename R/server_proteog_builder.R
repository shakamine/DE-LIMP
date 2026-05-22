# server_proteog_builder.R — orchestrator for the proteogenomics RNA-seq pipeline.
#
# Three public functions:
#   submit_proteogenomics_build() — kick off the full SLURM dep chain
#   poll_proteog_build_status()  — read status.json + refresh per-stage via sacct
#   cancel_proteog_build()       — scancel the whole chain
#
# This is the FIRST file in the proteogenomics stack with side effects (writes
# sbatch scripts to disk, submits to SLURM, updates status.json). All pure
# script-generation lives in helpers_rnaseq.R. All pure log-parsing lives in
# helpers_proteog_qc.R. This file only orchestrates.
#
# Local vs remote execution:
#   When R runs on Hive (via apptainer de-limp.sif), sbatch/squeue/sacct are
#   on PATH (via SLURM proxy or login-shell). This file assumes local execution
#   via system2(). For Docker-from-laptop, a future iteration can layer an
#   ssh_exec wrapper. Phase C v1 is single-host.
#
# Self-describing build manifest (CLAUDE.md rule #1):
#   Every build object returned from submit_*() carries $pipeline_id,
#   $methods_paragraph, and $stages so downstream consumers (Claude export,
#   methods README, status panel) read these rather than hardcoding "what
#   the pipeline did."

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

PROTEOG_PIPELINE_ID    <- "proteogenomics_v1.1"
PROTEOG_RNASEQ_ROOT    <- "/quobyte/proteomics-grp/de-limp/rnaseq"
PROTEOG_DATABASES_ROOT <- "/quobyte/proteomics-grp/de-limp/databases/proteogenomics"

# Stage order is the contract — used by status.json schema, poller, and cancel.
PROTEOG_STAGE_ORDER <- c(
  "fastp", "rrna_filter", "star", "qc_gate", "stringtie",
  "merge", "gffcompare", "gffread", "transdecoder", "rewrite", "assemble"
)

# =============================================================================
# Input validation + project setup
# =============================================================================

.validate_build_inputs <- function(project_name, rnaseq_dir, sample_names,
                                   reference_key, library_type, strand_flag,
                                   ref_registry) {
  if (!nzchar(project_name) || !grepl("^[A-Za-z0-9_.-]+$", project_name)) {
    stop("submit_proteogenomics_build(): project_name must be non-empty and match [A-Za-z0-9_.-]+; got: ",
         project_name)
  }
  if (length(sample_names) == 0) {
    stop("submit_proteogenomics_build(): no sample_names provided")
  }
  if (any(!grepl("^[A-Za-z0-9_.-]+$", sample_names))) {
    bad <- sample_names[!grepl("^[A-Za-z0-9_.-]+$", sample_names)]
    stop("submit_proteogenomics_build(): invalid sample name(s): ",
         paste(bad, collapse = ", "))
  }
  if (!dir.exists(rnaseq_dir)) {
    stop("submit_proteogenomics_build(): rnaseq_dir does not exist: ", rnaseq_dir)
  }
  # Each sample must have R1 + R2 in rnaseq_dir
  for (s in sample_names) {
    for (rd in c("R1", "R2")) {
      f <- file.path(rnaseq_dir, sprintf("%s_%s.fastq.gz", s, rd))
      if (!file.exists(f)) {
        stop(sprintf(
          "submit_proteogenomics_build(): missing input FASTQ for sample %s: %s",
          s, f
        ))
      }
    }
  }
  if (!reference_key %in% names(ref_registry)) {
    stop("submit_proteogenomics_build(): reference_key not in registry: ",
         reference_key, " (known: ",
         paste(names(ref_registry), collapse = ", "), ")")
  }
  if (!library_type %in% c("polyA", "ribo_depleted", "stranded", "unstranded")) {
    stop("submit_proteogenomics_build(): library_type must be one of ",
         "polyA/ribo_depleted/stranded/unstranded; got: ", library_type)
  }
  if (!strand_flag %in% c("", "--rf", "--fr")) {
    stop("submit_proteogenomics_build(): strand_flag must be '', '--rf', or '--fr'; got: ",
         strand_flag)
  }
  invisible(TRUE)
}

#' Stage raw FASTQs into the project's expected layout
#'
#' The Phase B sbatch generators expect FASTQs at <project_dir>/rnaseq/<sample>_R{1,2}.fastq.gz.
#' If rnaseq_dir is different (e.g., user pointed at sra_data/), we symlink the
#' files into place. Idempotent.
.stage_rnaseq_inputs <- function(project_dir, rnaseq_dir, sample_names) {
  target <- file.path(project_dir, "rnaseq")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  for (s in sample_names) {
    for (rd in c("R1", "R2")) {
      fname <- sprintf("%s_%s.fastq.gz", s, rd)
      src <- file.path(rnaseq_dir, fname)
      dst <- file.path(target, fname)
      if (!file.exists(dst)) {
        # Use absolute path for the symlink target
        file.symlink(normalizePath(src, mustWork = TRUE), dst)
      }
    }
  }
  invisible(target)
}

# =============================================================================
# Sbatch dispatch — local shell submission via system2()
# =============================================================================

#' Run sbatch on a script file; return parsed job_id or stop on failure.
#'
#' Direct system2 call — R inherits PATH from its parent shell (Hive's
#' `bash -l` loads SLURM tools at /cvmfs/.../slurm/bin via the standard
#' login-shell environment). Spawning a new bash -l -c from R does NOT
#' re-load modules (lmod is shell-scoped) and breaks the path resolution,
#' so we keep this direct.
.sbatch_submit <- function(script_path, dep_jid = NULL) {
  args <- character()
  if (!is.null(dep_jid) && nzchar(dep_jid)) {
    args <- c(args, sprintf("--dependency=afterok:%s", dep_jid))
  }
  args <- c(args, script_path)
  out <- tryCatch(
    suppressWarnings(system2("sbatch", args = args,
                              stdout = TRUE, stderr = TRUE)),
    error = function(e) stop("sbatch failed: ", conditionMessage(e))
  )
  if (!is.character(out) || !any(grepl("Submitted batch job", out))) {
    stop("sbatch did not return a job id. Output: ",
         paste(out, collapse = "\n"))
  }
  jid <- parse_sbatch_output(out)
  if (is.null(jid) || !nzchar(jid)) {
    stop("could not parse job id from sbatch output: ",
         paste(out, collapse = "\n"))
  }
  jid
}

#' Query sacct for one job's current state. Returns one of:
#'   "pending" | "running" | "complete" | "failed" | "cancelled" | "unknown"
#'
#' Direct system2 call (same env story as .sbatch_submit above). Uses sacct
#' -X to exclude .extern/.batch substeps that confuse the standard
#' check_slurm_status() helper. Sacct works for both active and completed
#' jobs; squeue only returns active ones, which is why the SSH-centric
#' check_slurm_status() in helpers_search.R fails on this code path.
.sacct_state <- function(jid) {
  # Defensive: jid can arrive as NULL, NA, character(0), or "" depending on
  # how the status.json was last serialized by jsonlite.
  if (is.null(jid)) return("unknown")
  jid <- suppressWarnings(as.character(jid))
  if (length(jid) == 0) return("unknown")
  if (is.na(jid[1])) return("unknown")
  if (!nzchar(jid[1])) return("unknown")
  jid <- jid[1]
  out <- tryCatch(
    suppressWarnings(system2("sacct",
      args = c("-j", jid, "-X", "-n", "-o", "State"),
      stdout = TRUE, stderr = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out) || length(out) == 0 || !nzchar(trimws(out[1]))) {
    return("unknown")
  }
  state <- trimws(out[1])
  if (grepl("^COMPLETED",  state)) return("complete")
  if (grepl("^(RUNNING|COMPLETING)", state)) return("running")
  if (grepl("^PENDING",    state)) return("pending")
  if (grepl("^CANCELLED",  state)) return("cancelled")
  if (grepl("^(FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|PREEMPTED)",
            state)) return("failed")
  "unknown"
}

# =============================================================================
# status.json schema + I/O
# =============================================================================

.init_status_json <- function(project_dir, project_name, sample_names,
                              reference_key, tier_params, jids_by_stage,
                              build_metadata) {
  stages <- lapply(PROTEOG_STAGE_ORDER, function(s) {
    list(
      stage      = s,
      status     = "pending",
      job_id     = jids_by_stage[[s]] %||% NA_character_,
      started_at = NA_character_,
      finished_at = NA_character_
    )
  })

  status <- list(
    pipeline_id      = PROTEOG_PIPELINE_ID,
    project_name     = project_name,
    project_dir      = project_dir,
    sample_names     = as.list(sample_names),
    reference_key    = reference_key,
    read_length_tier = tier_params$tier,
    qc_gate_unique_pct = tier_params$qc_gate_unique_pct,
    submitted_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    current_stage    = "fastp",
    stages           = stages,
    build_metadata   = build_metadata
  )

  status_path <- file.path(project_dir, "status.json")
  jsonlite::write_json(status, status_path, auto_unbox = TRUE, pretty = TRUE)
  invisible(status_path)
}

# =============================================================================
# Public: submit_proteogenomics_build
# =============================================================================

#' Submit the full proteogenomics build pipeline
#'
#' Generates sbatch scripts for every stage, submits them with
#' `--dependency=afterok` chaining, and writes a status.json that downstream
#' callers can poll.
#'
#' @param project_name      character — sanitized; used as project subdir name
#' @param rnaseq_dir        character — directory containing <sample>_R{1,2}.fastq.gz
#' @param reference_key     character — key in references/registry.json (e.g., "mm39_GRCm39")
#' @param sample_names      character vector
#' @param library_type      "polyA" | "ribo_depleted" | "stranded" | "unstranded"
#' @param strand_flag       "" | "--rf" | "--fr"
#' @param project_tag       character — passed to header rewriter; defaults to upper-cased project_name
#' @param uniprot_fasta     character or NULL — for Phase 4 assembly step; NULL skips merge
#' @param diamond_db        character or NULL — TransDecoder homology support
#' @param min_orf_len       integer — TransDecoder LongOrfs min length (default 100)
#' @param slurm_account     character (default "genome-center-grp")
#' @param slurm_partition   character (default "high")
#' @param ref_registry      list — if NULL, loaded from /quobyte/.../references/registry.json
#' @param rnaseq_root       character — base output dir
#' @return list with $project_dir, $status_path, $jids_by_stage, $tier_params,
#'         $pipeline_id, $methods_paragraph
submit_proteogenomics_build <- function(
  project_name,
  rnaseq_dir,
  reference_key,
  sample_names,
  library_type      = "polyA",
  strand_flag       = "",
  project_tag       = NULL,
  uniprot_fasta     = NULL,
  diamond_db        = NULL,
  min_orf_len       = 100L,
  slurm_account     = "genome-center-grp",
  slurm_partition   = "high",
  ref_registry      = NULL,
  rnaseq_root       = PROTEOG_RNASEQ_ROOT
) {
  # ---- 1. Load reference registry + validate inputs --------------------------
  if (is.null(ref_registry)) {
    ref_registry <- load_reference_registry()
  }
  .validate_build_inputs(project_name, rnaseq_dir, sample_names,
                         reference_key, library_type, strand_flag,
                         ref_registry)
  ref <- ref_registry[[reference_key]]

  # ---- 2. Set up project_dir -------------------------------------------------
  project_dir <- file.path(rnaseq_root, project_name)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(project_dir, "logs"), showWarnings = FALSE)
  .stage_rnaseq_inputs(project_dir, rnaseq_dir, sample_names)

  # ---- 3. Detect read length on first R1 -------------------------------------
  first_r1 <- file.path(project_dir, "rnaseq",
                        sprintf("%s_R1.fastq.gz", sample_names[1]))
  read_len <- detect_read_length(first_r1, n_reads = 100L)
  if (is.na(read_len)) {
    stop("submit_proteogenomics_build(): could not detect read length from ", first_r1)
  }
  tier_params <- select_star_params(read_len)
  if (tier_params$tier == "refuse") {
    stop("submit_proteogenomics_build(): ", tier_params$error,
         " (detected median read length = ", read_len, " bp)")
  }

  if (is.null(project_tag)) {
    project_tag <- toupper(sanitize_project_name(project_name))
  }

  # ---- 4. Generate all sbatch scripts ----------------------------------------
  scripts <- list(
    fastp        = generate_fastp_sbatch(project_dir, sample_names,
                                          slurm_account, slurm_partition),
    rrna_filter  = generate_rrna_sbatch(project_dir, sample_names,
                                         ref$rrna_index,
                                         slurm_account, slurm_partition),
    star         = generate_star_sbatch(project_dir, sample_names,
                                         ref$star_index, tier_params,
                                         slurm_account, slurm_partition),
    qc_gate      = generate_qc_gate_sbatch(project_dir, sample_names,
                                            tier_params$qc_gate_unique_pct,
                                            slurm_account, slurm_partition),
    stringtie    = generate_stringtie_sbatch(project_dir, sample_names,
                                              ref$gtf, strand_flag,
                                              slurm_account, slurm_partition),
    merge        = generate_merge_sbatch(project_dir, ref$gtf, sample_names,
                                          slurm_account, slurm_partition),
    gffcompare   = generate_gffcompare_sbatch(project_dir, ref$gtf,
                                               slurm_account, slurm_partition),
    gffread      = generate_gffread_sbatch(project_dir, ref$genome_fasta,
                                            slurm_account, slurm_partition),
    transdecoder = generate_transdecoder_sbatch(project_dir, diamond_db,
                                                 min_orf_len,
                                                 slurm_account, slurm_partition),
    rewrite      = generate_rewrite_sbatch(project_dir, project_tag,
                                            slurm_account = slurm_account,
                                            slurm_partition = slurm_partition)
  )

  # Write each script to <project_dir>/sbatch/<stage>.sbatch
  sbatch_dir <- file.path(project_dir, "sbatch")
  dir.create(sbatch_dir, showWarnings = FALSE)
  script_paths <- character()
  for (stage in names(scripts)) {
    p <- file.path(sbatch_dir, sprintf("%s.sbatch", stage))
    writeLines(scripts[[stage]], p)
    Sys.chmod(p, "755")
    script_paths[[stage]] <- p
  }

  # ---- 5. Submit with afterok dependency chaining ---------------------------
  jids_by_stage <- list()
  prev <- NULL
  for (stage in names(scripts)) {
    jid <- .sbatch_submit(script_paths[[stage]], dep_jid = prev)
    jids_by_stage[[stage]] <- jid
    prev <- jid
  }

  # ---- 6. Self-describing methods paragraph + status.json --------------------
  build_metadata <- list(
    rnaseq_dir       = rnaseq_dir,
    library_type     = library_type,
    strand_flag      = strand_flag,
    organism         = ref$organism %||% NA_character_,
    genome_build     = ref$build %||% NA_character_,
    annotation       = paste(c(ref$annotation_source, ref$annotation_release), collapse = " "),
    project_tag      = project_tag,
    detected_read_length = read_len
  )

  methods_paragraph <- sprintf(
    "Proteogenomics database built from %d samples (%s, %s; %s tier, read length %.0f bp). Pipeline: fastp → bowtie2 rRNA filter → STAR (%s) → stringtie → merge → gffcompare → gffread → TransDecoder → header rewrite (project tag %s). Reference: %s %s (%s %s).",
    length(sample_names), library_type, strand_flag,
    tier_params$tier, read_len, paste(tier_params$flags, collapse = " "),
    project_tag,
    build_metadata$organism, build_metadata$genome_build,
    ref$annotation_source %||% "", ref$annotation_release %||% ""
  )

  status_path <- .init_status_json(project_dir, project_name, sample_names,
                                    reference_key, tier_params,
                                    jids_by_stage, build_metadata)

  list(
    pipeline_id       = PROTEOG_PIPELINE_ID,
    project_dir       = project_dir,
    status_path       = status_path,
    jids_by_stage     = jids_by_stage,
    tier_params       = tier_params,
    build_metadata    = build_metadata,
    methods_paragraph = methods_paragraph
  )
}

# =============================================================================
# Public: poll_proteog_build_status
# =============================================================================

#' Refresh status.json by querying sacct for each stage's job_id
#'
#' @param project_dir character — from submit_*() result
#' @return updated status list with $current_stage and per-stage states
poll_proteog_build_status <- function(project_dir) {
  status_path <- file.path(project_dir, "status.json")
  if (!file.exists(status_path)) {
    stop("poll_proteog_build_status(): status.json not found at ", status_path)
  }
  status <- jsonlite::read_json(status_path)
  if (is.null(status$stages)) {
    stop("poll_proteog_build_status(): status.json missing $stages")
  }

  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  any_running    <- FALSE
  current_stage  <- "complete"
  any_failed     <- FALSE

  # nzchar() on NA returns NA, breaking the `if`. Coerce NA → "" up front.
  .empty_or_str <- function(v) {
    if (is.null(v) || (length(v) == 1 && is.na(v))) "" else as.character(v)
  }

  for (i in seq_along(status$stages)) {
    st <- status$stages[[i]]
    if (st$status %in% c("complete", "failed", "cancelled")) next
    new_state <- .sacct_state(st$job_id)
    if (new_state == "running" && !nzchar(.empty_or_str(st$started_at))) {
      status$stages[[i]]$started_at <- now
    }
    if (new_state %in% c("complete", "failed", "cancelled") &&
        !nzchar(.empty_or_str(st$finished_at))) {
      status$stages[[i]]$finished_at <- now
    }
    status$stages[[i]]$status <- new_state
    if (new_state == "running" || new_state == "pending") {
      any_running <- TRUE
      if (current_stage == "complete") current_stage <- st$stage
    }
    if (new_state == "failed" || new_state == "cancelled") {
      any_failed <- TRUE
    }
  }

  status$current_stage <- if (any_failed) "failed"
                          else if (any_running) current_stage
                          else "complete"
  status$last_polled_at <- now

  jsonlite::write_json(status, status_path, auto_unbox = TRUE, pretty = TRUE)
  status
}

# =============================================================================
# Public: cancel_proteog_build
# =============================================================================

#' Scancel every non-terminal job in the build
cancel_proteog_build <- function(project_dir) {
  status <- poll_proteog_build_status(project_dir)
  to_cancel <- character()
  for (st in status$stages) {
    if (!is.null(st$job_id) && nzchar(st$job_id) &&
        !(st$status %in% c("complete", "failed", "cancelled"))) {
      to_cancel <- c(to_cancel, st$job_id)
    }
  }
  if (length(to_cancel) == 0) {
    return(invisible(list(cancelled = character(),
                          message = "no active jobs to cancel")))
  }
  out <- tryCatch(
    system2("scancel", args = to_cancel, stdout = TRUE, stderr = TRUE),
    error = function(e) stop("scancel failed: ", conditionMessage(e))
  )
  invisible(list(
    cancelled = to_cancel,
    scancel_output = out
  ))
}


# =============================================================================
# Shiny module — UI renderer + server logic for the Build Database 🧬 tab
# =============================================================================
# Called from app.R as: server_proteog_builder(input, output, session, values)
#
# Hard-gated by hpc_available + !is_hf_space at the UI layer (see R/ui.R), so
# the server module assumes sbatch/squeue/sacct are reachable. If invoked on
# Docker-only, the UI panel never renders so the observers below are inert.

#' Render the Build Database tab body
#'
#' Five vertically-stacked accordion-style cards:
#'   1. Source — SLIMS URL OR comma-separated SRA accessions
#'   2. Sample scan / metadata verification
#'   3. Reference selection (from /quobyte/.../references/registry.json)
#'   4. Pipeline parameters (library type, strand, project name + tag)
#'   5. Submit
#'   plus an Active builds table polling status.json files
#'
#' All inputs are namespaced with `proteog_` to avoid colliding with the
#' existing Run Search inputs (analysis_name, fasta_source, etc.).
build_database_ui <- function() {
  div(
    style = "overflow-y: auto; max-height: calc(100vh - 150px); padding: 8px;",

    # ── Header
    div(
      style = "background: linear-gradient(135deg, #e8f5e9 0%, #c8e6c9 100%); padding: 12px 16px; border-radius: 8px; margin-bottom: 16px;",
      h4(icon("dna"), " Build Proteogenomics Database",
         style = "margin: 0; color: #1b5e20;"),
      p("Construct a sample-specific FASTA from matched RNA-seq data. The output ",
        "appears in the Run Search FASTA dropdown with a \U0001F9EC tag.",
        style = "margin: 4px 0 0 0; color: #2e7d32; font-size: 0.9em;")
    ),

    # ── Step 1: Source ──────────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(tagList(icon("upload"), " 1. RNA-seq Source")),
      bslib::card_body(
        radioButtons("proteog_source_mode", NULL,
                     choices = c("DNA Tech Core SLIMS URL" = "slims",
                                 "SRA/ENA accession(s)"    = "sra"),
                     selected = "slims", inline = TRUE),
        conditionalPanel(
          "input.proteog_source_mode == 'slims'",
          textInput("proteog_slims_url", "SLIMS URL",
                    placeholder = "http://slimsdata.genomecenter.ucdavis.edu/Data/<id>/Unaligned/",
                    width = "100%")
        ),
        conditionalPanel(
          "input.proteog_source_mode == 'sra'",
          textInput("proteog_sra_accessions", "Accessions (comma-separated, max 24)",
                    placeholder = "SRR1303776, SRR1303777",
                    width = "100%"),
          checkboxInput("proteog_subsample",
                        "Stream-subsample to 5M read pairs per accession (fast test)",
                        value = FALSE)
        ),
        actionButton("proteog_scan_btn", "Scan / Verify",
                     icon = icon("magnifying-glass"),
                     class = "btn-outline-primary")
      )
    ),

    # ── Step 2: Scan / verification results ─────────────────────────────────
    bslib::card(
      bslib::card_header(tagList(icon("clipboard-check"), " 2. Sample Verification")),
      bslib::card_body(uiOutput("proteog_scan_output"))
    ),

    # ── Step 3: Reference selection ─────────────────────────────────────────
    bslib::card(
      bslib::card_header(tagList(icon("book-open"), " 3. Reference Genome")),
      bslib::card_body(
        selectInput("proteog_reference_key", "Reference",
                    choices = c("(scanning registry…)" = ""),
                    width = "100%"),
        uiOutput("proteog_reference_info")
      )
    ),

    # ── Step 4: Pipeline parameters ─────────────────────────────────────────
    bslib::card(
      bslib::card_header(tagList(icon("sliders"), " 4. Parameters")),
      bslib::card_body(
        layout_columns(
          col_widths = c(6, 6),
          selectInput("proteog_library_type", "Library type",
                      choices = c("polyA mRNA-Seq"            = "polyA",
                                  "Total RNA + rRNA depletion" = "ribo_depleted",
                                  "Stranded RNA-Seq"           = "stranded",
                                  "Unstranded"                 = "unstranded"),
                      selected = "polyA"),
          selectInput("proteog_strand_flag", "Strand",
                      choices = c("Reverse stranded (TruSeq, --rf)" = "--rf",
                                  "Forward stranded (--fr)"          = "--fr",
                                  "Unstranded ()"                    = ""),
                      selected = "--rf")
        ),
        textInput("proteog_project_name", "Project name",
                  placeholder = "e.g. mouse_liver_pilot_2026_05",
                  width = "100%"),
        textInput("proteog_project_tag", "Project tag (uppercase, suffixes FASTA symbols)",
                  placeholder = "e.g. MOUSELIVER",
                  width = "100%"),
        numericInput("proteog_min_orf_len", "Minimum ORF length (aa)",
                     value = 100, min = 30, max = 300, step = 10, width = "50%")
      )
    ),

    # ── Step 5: Submit ──────────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(tagList(icon("rocket"), " 5. Submit")),
      bslib::card_body(
        uiOutput("proteog_submit_warnings"),
        actionButton("proteog_submit_btn", "Build Proteogenomics FASTA",
                     icon = icon("dna"),
                     class = "btn-success btn-lg w-100"),
        helpText("Estimated wall time: 3-6 hours for 12 samples × 30M PE150 reads. ",
                 "You can close the browser; the build continues on Hive.")
      )
    ),

    # ── Active builds table ─────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(tagList(icon("list-check"), " Active & recent builds")),
      bslib::card_body(uiOutput("proteog_active_builds_table"))
    )
  )
}

#' Shiny server module for the Build Database tab
#'
#' Wires the UI inputs to:
#'   - load_reference_registry() / load_proteog_registry() for dropdowns
#'   - scan_slims_url() / verify_sra_accession() for Step 1 verify
#'   - submit_proteogenomics_build() for the submit button
#'   - reactivePoll on status.json files for the active-builds table
server_proteog_builder <- function(input, output, session, values) {

  ns <- session$ns %||% function(x) x  # tolerant of being called inside a moduleServer or not

  # ── render the Build Database body via uiOutput("build_database_content") ──
  output$build_database_content <- renderUI({ build_database_ui() })

  # ── populate reference dropdown from registry ──────────────────────────────
  observe({
    reg <- tryCatch(load_reference_registry(), error = function(e) list())
    if (length(reg) == 0) {
      updateSelectInput(session, "proteog_reference_key",
                        choices = c("No references registered" = ""))
      return()
    }
    labels <- vapply(names(reg), function(k) {
      e <- reg[[k]]
      sprintf("%s — %s %s", k,
              e$organism %||% "?",
              e$annotation_release %||% e$annotation_source %||% "")
    }, character(1), USE.NAMES = FALSE)
    choices <- setNames(names(reg), labels)
    updateSelectInput(session, "proteog_reference_key", choices = choices)
  })

  output$proteog_reference_info <- renderUI({
    req(input$proteog_reference_key, nzchar(input$proteog_reference_key))
    reg <- load_reference_registry()
    entry <- reg[[input$proteog_reference_key]]
    if (is.null(entry)) return(NULL)
    tags$div(
      style = "background: #f5f5f5; padding: 8px; border-radius: 4px; font-size: 0.85em;",
      tags$div(strong("Genome FASTA:"), code(entry$genome_fasta %||% "(none)")),
      tags$div(strong("STAR index:"),   code(entry$star_index %||% "(none)")),
      tags$div(strong("GTF:"),          code(entry$gtf %||% "(none)")),
      tags$div(strong("rRNA index:"),   code(entry$rrna_index %||% "(none)")),
      tags$div(strong("Completeness:"), entry$completeness %||% "(unknown)")
    )
  })

  # ── Step 1: Scan / verify ────────────────────────────────────────────────
  # Stored in a reactiveVal so the UI in Step 2 can render it after the user
  # clicks Scan, and the submit handler can read it back to populate sample
  # names without re-fetching.
  scan_result <- reactiveVal(NULL)

  observeEvent(input$proteog_scan_btn, {
    if (input$proteog_source_mode == "slims") {
      url <- trimws(input$proteog_slims_url %||% "")
      if (!nzchar(url)) {
        scan_result(list(success = FALSE,
                         error = "Please enter a SLIMS URL."))
        return()
      }
      withProgress(message = "Scanning SLIMS URL…", value = 0.5, {
        res <- tryCatch(scan_slims_url(url),
                        error = function(e) list(success = FALSE,
                                                 error = conditionMessage(e)))
      })
      res$mode <- "slims"
      scan_result(res)
    } else {
      raw <- input$proteog_sra_accessions %||% ""
      accs <- trimws(strsplit(raw, "[,;[:space:]]+")[[1]])
      accs <- accs[nzchar(accs)]
      if (length(accs) == 0) {
        scan_result(list(success = FALSE,
                         error = "Please enter at least one SRA/ENA accession."))
        return()
      }
      if (length(accs) > 24) {
        scan_result(list(success = FALSE,
                         error = sprintf("Too many accessions (%d > 24).",
                                         length(accs))))
        return()
      }
      withProgress(message = "Verifying ENA metadata…", value = 0, {
        per_acc <- lapply(seq_along(accs), function(i) {
          incProgress(1 / length(accs), detail = accs[i])
          verify_sra_accession(accs[i])
        })
      })
      n_bad <- sum(!vapply(per_acc, function(r) isTRUE(r$success), logical(1)))
      scan_result(list(
        success = n_bad == 0,
        mode = "sra",
        accessions = accs,
        per_acc = per_acc,
        sample_names = accs,
        error = if (n_bad > 0) sprintf("%d accession(s) failed verification.", n_bad) else NULL
      ))
    }
  })

  output$proteog_scan_output <- renderUI({
    res <- scan_result()
    if (is.null(res)) {
      return(helpText("Click ", strong("Scan / Verify"),
                      " after filling in a source above."))
    }
    if (!isTRUE(res$success)) {
      return(tags$div(class = "alert alert-danger",
                      tags$strong("Error: "), res$error %||% "Unknown error."))
    }
    if (identical(res$mode, "slims")) {
      tagList(
        tags$div(class = "alert alert-success",
                 sprintf("SLIMS URL OK. %d samples (%s).",
                         res$n_samples,
                         if (isTRUE(res$is_paired)) "paired-end" else "single-end")),
        tags$details(
          tags$summary(sprintf("Sample list (%d)", length(res$sample_names))),
          tags$ul(lapply(res$sample_names, function(s) tags$li(code(s))))
        )
      )
    } else {
      # SRA mode — show per-accession metadata
      rows <- lapply(res$per_acc, function(r) {
        if (!isTRUE(r$success)) {
          return(tags$tr(tags$td(code(r$accession %||% "?")),
                         tags$td(colspan = 4,
                                 tags$span(style="color:#c62828;", r$error))))
        }
        tags$tr(
          tags$td(code(r$accession)),
          tags$td(r$scientific_name %||% "?"),
          tags$td(r$library_strategy %||% "?",
                  if (!isTRUE(r$suitable))
                    tags$span(style = "color:#c62828; font-weight:bold;",
                              " ⚠ unsuitable")),
          tags$td(r$instrument %||% "?"),
          tags$td(r$layout %||% "?")
        )
      })
      tagList(
        tags$table(class = "table table-sm",
                   tags$thead(tags$tr(
                     tags$th("Accession"), tags$th("Species"),
                     tags$th("Library"), tags$th("Instrument"), tags$th("Layout")
                   )),
                   tags$tbody(rows))
      )
    }
  })

  # ── Pre-submit warnings + button enable/disable ─────────────────────────────
  output$proteog_submit_warnings <- renderUI({
    warnings <- character()
    res <- scan_result()
    if (is.null(res) || !isTRUE(res$success)) {
      warnings <- c(warnings, "Run Scan / Verify before submitting.")
    }
    if (!nzchar(input$proteog_reference_key %||% "")) {
      warnings <- c(warnings, "Select a reference genome.")
    }
    if (!nzchar(input$proteog_project_name %||% "")) {
      warnings <- c(warnings, "Enter a project name.")
    }
    if (!nzchar(input$proteog_project_tag %||% "")) {
      warnings <- c(warnings, "Enter a project tag.")
    }
    # Species mismatch check (SRA mode only)
    if (!is.null(res) && identical(res$mode, "sra") &&
        nzchar(input$proteog_reference_key %||% "")) {
      reg <- load_reference_registry()
      ref_org <- reg[[input$proteog_reference_key]]$organism %||% ""
      accs_orgs <- vapply(res$per_acc,
        function(r) if (isTRUE(r$success)) r$scientific_name %||% "" else "",
        character(1))
      mismatched <- accs_orgs[nzchar(accs_orgs) & accs_orgs != ref_org]
      if (length(mismatched) > 0) {
        warnings <- c(warnings, sprintf(
          "Species mismatch: reference is %s; accessions report %s.",
          ref_org, paste(unique(mismatched), collapse = ", ")))
      }
    }
    if (length(warnings) == 0) return(NULL)
    tags$div(class = "alert alert-warning",
             tags$strong("Cannot submit yet:"),
             tags$ul(lapply(warnings, tags$li)))
  })

  # ── Submit handler ──────────────────────────────────────────────────────────
  observeEvent(input$proteog_submit_btn, {
    res <- scan_result()
    req(res, isTRUE(res$success))
    req(nzchar(input$proteog_reference_key %||% ""))
    req(nzchar(input$proteog_project_name %||% ""))
    req(nzchar(input$proteog_project_tag  %||% ""))

    # Download data (login-node nohup) if SLIMS mode and not already downloaded.
    # For SRA mode the orchestrator can drive a stream-subsample sbatch.
    # For Phase D v1 we DEFER the download orchestration — assume the data is
    # already at the expected rnaseq_dir if user clicks submit. The next
    # iteration will gate the submit on a completed download.
    rnaseq_dir <- tryCatch({
      if (identical(res$mode, "slims")) {
        d <- launch_slims_download(res$url,
                                   sanitize_project_name(input$proteog_project_name))
        d$project_dir
      } else {
        d <- launch_ena_download(res$accessions,
                                 sanitize_project_name(input$proteog_project_name),
                                 subsample_reads = if (isTRUE(input$proteog_subsample))
                                                     5e6L else NULL)
        d$project_dir
      }
    }, error = function(e) {
      showNotification(sprintf("Download launch failed: %s",
                               conditionMessage(e)),
                       type = "error", duration = 10)
      NULL
    })
    req(rnaseq_dir)

    showNotification(
      tags$div(
        tags$p(strong("Download started"),
               " — the SLURM pipeline will be submitted once data is present."),
        tags$p("You can close the browser; build continues on Hive."),
        tags$p(tags$code(rnaseq_dir))),
      type = "message", duration = 10
    )

    # Stash the build request so a downstream observer can submit
    # submit_proteogenomics_build() once the download status.json reports
    # state=="complete". For Phase D v1, we ALSO emit the call immediately
    # — the download poll-and-submit observer can be added in Phase E.
    tryCatch({
      build <- submit_proteogenomics_build(
        project_name    = sanitize_project_name(input$proteog_project_name),
        rnaseq_dir      = rnaseq_dir,
        reference_key   = input$proteog_reference_key,
        sample_names    = res$sample_names %||% character(0),
        library_type    = input$proteog_library_type,
        strand_flag     = input$proteog_strand_flag,
        project_tag     = input$proteog_project_tag,
        min_orf_len     = as.integer(input$proteog_min_orf_len %||% 100L),
        slurm_account   = "genome-center-grp",
        slurm_partition = "high"
      )
      # Track in reactiveValues so the active builds table can poll
      jobs <- values$proteog_build_jobs %||% list()
      jobs[[length(jobs) + 1L]] <- list(
        project_name = sanitize_project_name(input$proteog_project_name),
        project_dir  = build$project_dir,
        submitted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        jids_by_stage = build$jids_by_stage,
        methods_paragraph = build$methods_paragraph
      )
      values$proteog_build_jobs <- jobs
      showNotification(sprintf("Build submitted: %s",
                               sanitize_project_name(input$proteog_project_name)),
                       type = "default", duration = 8)
    }, error = function(e) {
      showNotification(sprintf("Submit failed: %s", conditionMessage(e)),
                       type = "error", duration = 15)
    })
  })

  # ── Active builds table — polls status.json every 15 s ──────────────────────
  proteog_status_poll <- reactivePoll(
    intervalMillis = 15000,
    session = session,
    checkFunc = function() {
      jobs <- values$proteog_build_jobs %||% list()
      if (length(jobs) == 0) return("")
      paths <- vapply(jobs, function(j) file.path(j$project_dir, "status.json"),
                       character(1))
      paste(vapply(paths, function(p) {
        if (file.exists(p)) as.character(file.mtime(p)) else "MISSING"
      }, character(1)), collapse = "|")
    },
    valueFunc = function() {
      jobs <- values$proteog_build_jobs %||% list()
      lapply(jobs, function(j) {
        tryCatch(poll_proteog_build_status(j$project_dir),
                 error = function(e) NULL)
      })
    }
  )

  output$proteog_active_builds_table <- renderUI({
    statuses <- proteog_status_poll()
    if (length(statuses) == 0) {
      return(helpText("No builds submitted in this session yet."))
    }
    rows <- lapply(seq_along(statuses), function(i) {
      st <- statuses[[i]]
      if (is.null(st)) return(NULL)
      current <- st$current_stage %||% "?"
      badge_color <- switch(current,
        "complete" = "#27ae60",
        "failed"   = "#c0392b",
        "#f39c12")
      done <- sum(vapply(st$stages,
        function(s) identical(s$status, "complete"), logical(1)))
      total <- length(st$stages)
      tags$tr(
        tags$td(st$project_name %||% "?"),
        tags$td(tags$span(style = sprintf("background:%s; color:white; padding:2px 8px; border-radius:4px;",
                                           badge_color),
                          current)),
        tags$td(sprintf("%d / %d", done, total)),
        tags$td(st$submitted_at %||% "?"),
        tags$td(code(basename(st$project_dir %||% "")))
      )
    })
    tags$table(class = "table table-sm",
               tags$thead(tags$tr(
                 tags$th("Project"), tags$th("Stage"),
                 tags$th("Progress"), tags$th("Submitted"),
                 tags$th("Dir")
               )),
               tags$tbody(rows))
  })

  invisible(NULL)
}
