# =============================================================================
# helpers_site.R — Site-configurable defaults for non-UC-Davis deployments
# =============================================================================
#
# Background: DE-LIMP was originally built for the UC Davis Proteomics core
# facility's HPC layout. Several structural assumptions about paths, SLURM
# accounts, and shared resources got hardcoded across the codebase. UI defaults
# are user-overrideable (typed over in textboxes), but these structural values
# were buried in helpers and silently broke for non-UCD users.
#
# v3.10.15 — single source of truth for site-specific values. Defaults
# preserve historic UCD behavior. Non-UCD sites override via:
#
#   1. ~/.delimp_site.yaml (preferred — version this in your lab repo)
#   2. Environment variables (good for one-off overrides):
#         DELIMP_STORAGE_LOCAL, DELIMP_STORAGE_HPC,
#         DELIMP_SHARED_FASTA_LIB_LOCAL, DELIMP_SHARED_FASTA_LIB_HPC,
#         DELIMP_SHARED_ACTIVITY_LOG,
#         DELIMP_SLURM_ACCOUNT, DELIMP_SLURM_PARTITION,
#         DELIMP_SLURM_FALLBACK_ACCOUNT, DELIMP_SLURM_FALLBACK_PARTITION,
#         DELIMP_FASTA_DIR_HPC, DELIMP_FASTA_DIR_LOCAL,
#         DELIMP_GENE_MAP_DIRS (colon-separated)
#   3. Function arguments — every callsite still accepts an explicit override.
#
# Cached in globalenv so we don't re-read the YAML / env on every call.
# Invalidate via `delimp_site_invalidate()` (e.g. after editing settings).
# =============================================================================

#' Site-specific configuration as a single list.
#'
#' Defaults preserve UC Davis Proteomics core facility behavior — Brett's
#' lab keeps working unchanged. Other sites override via env vars or
#' `~/.delimp_site.yaml`. Cached on first call.
#'
#' @return list with named entries — see top-of-file comment for keys.
#' @export
delimp_site <- function() {
  cached <- get0(".delimp_site_cache", envir = globalenv(), ifnotfound = NULL)
  if (!is.null(cached)) return(cached)

  # Try ~/.delimp_site.yaml first (if `yaml` package available).
  yaml_path <- Sys.getenv("DELIMP_SITE_YAML",
    file.path(Sys.getenv("HOME"), ".delimp_site.yaml"))
  yaml_cfg <- list()
  if (file.exists(yaml_path) && requireNamespace("yaml", quietly = TRUE)) {
    yaml_cfg <- tryCatch(yaml::read_yaml(yaml_path),
      error = function(e) {
        message("[DE-LIMP] Could not parse ", yaml_path, ": ", e$message)
        list()
      })
  }

  # Helper: look up a key with priority env > yaml > UCD-default.
  pick <- function(env_key, yaml_key, ucd_default) {
    v <- Sys.getenv(env_key, "")
    if (nzchar(v)) return(v)
    if (!is.null(yaml_cfg[[yaml_key]])) return(as.character(yaml_cfg[[yaml_key]]))
    ucd_default
  }

  cfg <- list(
    # --- Storage layout: local mount <-> HPC native ---
    # translate_storage_path() rewrites between these two prefixes.
    storage_local = pick("DELIMP_STORAGE_LOCAL", "storage_local",
      "/Volumes/proteomics-grp/"),
    storage_hpc   = pick("DELIMP_STORAGE_HPC",   "storage_hpc",
      "/quobyte/proteomics-grp/"),

    # --- Shared FASTA library cache (lab-wide predicted spectral libraries) ---
    shared_fasta_lib_local = pick("DELIMP_SHARED_FASTA_LIB_LOCAL",
      "shared_fasta_lib_local", "/Volumes/proteomics-grp/dia-nn/fasta_library"),
    shared_fasta_lib_hpc   = pick("DELIMP_SHARED_FASTA_LIB_HPC",
      "shared_fasta_lib_hpc",   "/quobyte/proteomics-grp/dia-nn/fasta_library"),

    # --- DIA-NN shared install dir (for SIF / binaries) ---
    shared_diann_local = pick("DELIMP_SHARED_DIANN_LOCAL", "shared_diann_local",
      "/Volumes/proteomics-grp/dia-nn"),
    shared_diann_hpc   = pick("DELIMP_SHARED_DIANN_HPC",   "shared_diann_hpc",
      "/quobyte/proteomics-grp/dia-nn"),

    # --- Pre-staged FASTA download target (UniProt / NCBI saves here) ---
    fasta_dir_hpc   = pick("DELIMP_FASTA_DIR_HPC",   "fasta_dir_hpc",
      "/quobyte/proteomics-grp/de-limp/fasta"),
    fasta_dir_local = pick("DELIMP_FASTA_DIR_LOCAL", "fasta_dir_local",
      "/Volumes/proteomics-grp/de-limp/fasta"),

    # --- Cross-machine activity log ---
    shared_activity_log = pick("DELIMP_SHARED_ACTIVITY_LOG",
      "shared_activity_log", "/quobyte/proteomics-grp/de-limp/activity_log.csv"),

    # --- SLURM primary partition (where queue-submit lands by default) ---
    slurm_account   = pick("DELIMP_SLURM_ACCOUNT",   "slurm_account",
      "genome-center-grp"),
    slurm_partition = pick("DELIMP_SLURM_PARTITION", "slurm_partition",
      "high"),

    # --- SLURM fallback partition (auto-queue-switch destination) ---
    slurm_fallback_account   = pick("DELIMP_SLURM_FALLBACK_ACCOUNT",
      "slurm_fallback_account",   "publicgrp"),
    slurm_fallback_partition = pick("DELIMP_SLURM_FALLBACK_PARTITION",
      "slurm_fallback_partition", "low"),

    # --- gene_map.tsv search dirs for NCBI ortholog lookup ---
    # Colon-separated. These get tried in order during gene-symbol resolution.
    gene_map_dirs = strsplit(pick("DELIMP_GENE_MAP_DIRS", "gene_map_dirs",
      "/data/fasta:/quobyte/proteomics-grp/de-limp/fasta"), ":",
      fixed = TRUE)[[1]],

    # --- Source of truth for "is this a UCD site?" ---
    # Used to gate UCD-only features (auto-queue-switch, shared FASTA cache,
    # cluster usage history). Auto-detect via mount existence; can override
    # with DELIMP_UCD_MODE=true|false.
    is_ucd_site = {
      override <- Sys.getenv("DELIMP_UCD_MODE", "")
      if (nzchar(override)) {
        isTRUE(as.logical(override))
      } else {
        dir.exists("/quobyte/proteomics-grp") ||
          dir.exists("/Volumes/proteomics-grp")
      }
    },

    yaml_path = yaml_path
  )

  assign(".delimp_site_cache", cfg, envir = globalenv())
  cfg
}

#' Resolve the FASTA download/cache directory at runtime.
#'
#' Returns (in priority order) the first existing path among:
#'   1. `getOption("delimp.fasta_dir", NULL)` — explicit programmatic override
#'   2. `delimp_site()$fasta_dir_local` (default `/Volumes/proteomics-grp/de-limp/fasta`)
#'   3. `delimp_site()$fasta_dir_hpc`   (default `/quobyte/proteomics-grp/de-limp/fasta`)
#'   4. `~/.delimp_fasta` — created if missing
#'
#' Replaces the v3.10.x pattern of `getOption("delimp.fasta_dir",
#' default = "/quobyte/proteomics-grp/de-limp/fasta")` which silently
#' created `/quobyte/...` on non-UCD users' filesystems.
#'
#' @return Character path that exists (or was just created).
#' @export
resolve_fasta_dir <- function() {
  cfg <- delimp_site()
  candidates <- c(
    getOption("delimp.fasta_dir", NULL),
    cfg$fasta_dir_local,
    cfg$fasta_dir_hpc
  )
  candidates <- candidates[!vapply(candidates, is.null, logical(1)) & nzchar(candidates)]
  for (p in candidates) {
    if (dir.exists(p)) return(p)
  }
  # No configured path exists — fall back to a sane local default
  fallback <- file.path(Sys.getenv("HOME"), ".delimp_fasta")
  if (!dir.exists(fallback)) dir.create(fallback, recursive = TRUE, showWarnings = FALSE)
  fallback
}

#' Invalidate the site-config cache.
#' Call after editing ~/.delimp_site.yaml or env vars at runtime.
#' @export
delimp_site_invalidate <- function() {
  if (exists(".delimp_site_cache", envir = globalenv())) {
    rm(".delimp_site_cache", envir = globalenv())
  }
  invisible(NULL)
}
