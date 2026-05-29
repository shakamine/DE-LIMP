# =============================================================================
# helpers_coreomics.R — Integration with Adam Olshen's coreomics submission
# system (UC Davis core facility submissions portal).
#
# Feature is GATED on presence of a coreomics token. If neither
# `COREOMICS_TOKEN` env nor `~/.coreomics_token` exists, every helper here
# becomes a no-op or returns NULL — so UI conditionals like
# `if (coreomics_enabled()) { ... }` hide the feature entirely for
# deployments outside UCD Proteomics Core (Flinders, HF, etc.).
#
# Endpoints:
#   GET  <base>/submissions/?lab=PROTEOMICS&page=N&page_size=K
#                                            → DRF-paginated list
#   GET  <base>/submissions/<id>/            → submission detail
#                                              (.submission_data.samples)
#
# Auth: bearer token, `Authorization: Token <hex>` header.
# =============================================================================


# ── Gate / config ────────────────────────────────────────────────────────────

#' Is the coreomics integration available in this deployment?
#'
#' Returns TRUE iff a token is reachable via env var `COREOMICS_TOKEN` OR
#' a readable file at `~/.coreomics_token`. UI panels should be wrapped in
#' `if (coreomics_enabled()) { ... }` so non-UCD installs never see them.
#'
#' @return logical(1)
coreomics_enabled <- function() {
  nzchar(Sys.getenv("COREOMICS_TOKEN", "")) ||
    file.exists(path.expand("~/.coreomics_token"))
}

#' Base URL for the coreomics REST API. Defaults to UCD's deployment.
#' Override via env `COREOMICS_BASE_URL` for any other site that ever
#' deploys coreomics.
coreomics_base_url <- function() {
  url <- Sys.getenv("COREOMICS_BASE_URL", "https://ucdavis.coreomics.com/server/api")
  sub("/+$", "", url)
}

#' Lab filter passed to `?lab=<X>` on the submissions list endpoint.
coreomics_lab_filter <- function() {
  Sys.getenv("COREOMICS_LAB", "PROTEOMICS")
}

#' Load the bearer token from env or the file at `~/.coreomics_token`.
#' Trims whitespace. Returns NULL if neither source has a token (callers
#' should already have gated on `coreomics_enabled()` before reaching here).
coreomics_token_load <- function() {
  env <- Sys.getenv("COREOMICS_TOKEN", "")
  if (nzchar(env)) return(trimws(env))
  f <- path.expand("~/.coreomics_token")
  if (!file.exists(f)) return(NULL)
  tok <- trimws(readLines(f, warn = FALSE, n = 1))
  if (!nzchar(tok)) return(NULL)
  tok
}


# ── HTTP plumbing ────────────────────────────────────────────────────────────

#' Internal: GET a coreomics endpoint with auth + JSON parse.
#'
#' Uses `curl::curl_fetch_memory()` rather than a text-mode connection.
#' The latter trips `can only read from a binary connection` on R 4.x
#' when SSL+content-encoding is involved.
#'
#' @param path URL path relative to the base, e.g. "submissions/abc123/"
#' @param query Named list of query params (NULL → none)
#' @param timeout Seconds before bailing
#' @return Parsed list (from `jsonlite::fromJSON`), or NULL on auth / network failure.
coreomics_get <- function(path, query = NULL, timeout = 15) {
  if (!coreomics_enabled()) return(NULL)
  tok <- coreomics_token_load()
  if (is.null(tok)) return(NULL)

  url <- paste0(coreomics_base_url(), "/", sub("^/", "", path))
  if (!is.null(query) && length(query) > 0) {
    qs <- paste(sprintf("%s=%s", names(query),
                        vapply(query, function(v) utils::URLencode(as.character(v), reserved = TRUE),
                               character(1))),
                collapse = "&")
    url <- paste0(url, "?", qs)
  }

  tryCatch({
    h <- curl::new_handle(timeout = timeout)
    curl::handle_setheaders(h,
      "Authorization" = paste("Token", tok),
      "Accept" = "application/json")
    resp <- curl::curl_fetch_memory(url, handle = h)
    if (resp$status_code >= 400) {
      message("[coreomics] HTTP ", resp$status_code, " from ", url)
      return(NULL)
    }
    body <- rawToChar(resp$content)
    Encoding(body) <- "UTF-8"
    jsonlite::fromJSON(body, simplifyVector = FALSE)
  }, error = function(e) {
    message("[coreomics] GET ", url, " failed: ", conditionMessage(e))
    NULL
  })
}


# ── Top-level API ────────────────────────────────────────────────────────────

#' Path to the on-disk submission cache. One per user.
coreomics_cache_path <- function() {
  file.path(Sys.getenv("HOME"), ".delimp_coreomics_cache.rds")
}

#' Construct the canonical HIVE project directory for a coreomics
#' submission, per Adam Schaal's convention (2026-05-19):
#'   /quobyte/proteomics-grp/coreomics/projects/<YYYY>/<MM>/<submission_id>/
#'
#' The path is fixed for the lifetime of the submission (date + UUID).
#' Outputs that DE-LIMP places under this dir get picked up automatically
#' by coreomics' bioshare symlinks and customer-facing views.
#'
#' @param submission The submission dict (from coreomics_get_submission /
#'   coreomics_list_submissions). Must have $id and $submitted.
#' @return Character path, or NULL if either field is missing/malformed.
coreomics_project_dir <- function(submission) {
  if (is.null(submission)) return(NULL)
  sub_id <- submission$id %||% ""
  submitted <- submission$submitted %||% ""
  if (!nzchar(sub_id) || !nzchar(submitted)) return(NULL)
  ts <- tryCatch(as.POSIXct(submitted), error = function(e) NULL)
  if (is.null(ts) || is.na(ts)) return(NULL)
  sprintf("/quobyte/proteomics-grp/coreomics/projects/%s/%s/%s",
          format(ts, "%Y"), format(ts, "%m"), sub_id)
}

#' Compose a short, filesystem-safe analysis_name from a coreomics submission.
#' Prefers `internal_id` (PROT_0701); falls back to a sanitized snippet of
#' the description. Used to pre-fill the Analysis Name field.
coreomics_suggest_analysis_name <- function(submission) {
  if (is.null(submission)) return("")
  iid <- submission$internal_id %||% ""
  if (nzchar(iid)) return(iid)
  desc <- (submission$submission_data %||% list())$description %||% ""
  if (!nzchar(desc)) return("")
  # Sanitize: keep only alnum + _, cap at 40 chars
  clean <- gsub("[^A-Za-z0-9_]+", "_", desc)
  clean <- gsub("^_+|_+$", "", clean)
  substr(clean, 1, 40)
}

#' Build a human-readable label for a single submission, suitable for use
#' as a selectizeInput choice. Format:
#'   "PROT_0701 — Erik Chow (Paszek Lab) — 2026-05-19 — Mouse"
coreomics_submission_label <- function(s) {
  iid <- s$internal_id %||% s$id %||% ""
  submitter <- trimws(paste(s$first_name %||% "", s$last_name %||% ""))
  pi <- trimws(paste(s$pi_first_name %||% "", s$pi_last_name %||% ""))
  pi_part <- if (nzchar(pi)) sprintf(" (%s Lab)", pi) else ""
  date_part <- if (nzchar(s$submitted %||% "")) substr(s$submitted, 1, 10) else ""
  organism <- (s$submission_data %||% list())$organism %||% ""
  org_part <- if (nzchar(organism)) sprintf(" — %s", organism) else ""
  paste0(iid, " — ", submitter, pi_part,
         if (nzchar(date_part)) sprintf(" — %s", date_part) else "",
         org_part)
}

#' List all proteomics submissions, paginating internally. Returns a list
#' of submission summary objects (the `results` array from the API,
#' concatenated across pages).
#'
#' @param page_size How many per page. The API gets SLOWER per row at higher
#'   page sizes (server-side query time scales) — 100 is empirically optimal.
#' @param max_pages Safety cap so we don't accidentally fetch a runaway loop
#' @return List of submission dicts. Includes `submission_data` (with the
#'   `samples` array nested inside) — so per-submission detail calls aren't
#'   needed for the common match-raw-files use case.
coreomics_list_submissions <- function(page_size = 100, max_pages = 100) {
  if (!coreomics_enabled()) return(list())
  all <- list()
  for (page in seq_len(max_pages)) {
    d <- coreomics_get("submissions/",
                       query = list(lab = coreomics_lab_filter(),
                                    page = page, page_size = page_size))
    if (is.null(d) || length(d$results %||% list()) == 0) break
    all <- c(all, d$results)
    if (is.null(d$`next`) || !nzchar(d$`next`)) break
  }
  all
}

#' Load the submission list — prefers an on-disk cache (`~/.delimp_coreomics_cache.rds`)
#' when fresh enough; otherwise fetches from the API and writes the cache.
#'
#' Wall time savings: API fetch is ~30s for ~4400 submissions. Cache load is
#' < 1s. New sessions are effectively instant after the first warm.
#'
#' @param ttl_hours Cache lifetime. Defaults to 1 hour — captures new
#'   submissions on most working timescales without re-fetching on every
#'   new tab. Users can force a fresh fetch via the Refresh button.
#' @param force If TRUE, ignore the cache and fetch from the API. Used by
#'   the Refresh button.
#' @return List of submissions; also writes the cache file as a side effect
#'   when fetching from API.
coreomics_load_submissions <- function(ttl_hours = 1, force = FALSE) {
  if (!coreomics_enabled()) return(list())
  cache_file <- coreomics_cache_path()

  # Try the cache first
  if (!force && file.exists(cache_file)) {
    age_hours <- as.numeric(difftime(Sys.time(),
                                      file.info(cache_file)$mtime,
                                      units = "hours"))
    if (age_hours < ttl_hours) {
      cached <- tryCatch(readRDS(cache_file),
                         error = function(e) {
                           message("[coreomics] cache read failed: ",
                                   conditionMessage(e))
                           NULL
                         })
      if (!is.null(cached) && is.list(cached$submissions) &&
          length(cached$submissions) > 0) {
        return(cached$submissions)
      }
    }
  }

  # Fetch from API
  subs <- coreomics_list_submissions()
  if (length(subs) > 0) {
    tryCatch({
      saveRDS(list(submissions = subs, fetched_at = Sys.time()),
              file = cache_file)
    }, error = function(e) {
      message("[coreomics] cache write failed: ", conditionMessage(e))
    })
  }
  subs
}

#' Fetch one submission's full detail, including `submission_data.samples`.
#'
#' @param submission_id Coreomics submission ID, e.g. "5a8c37df3bd1"
#' @return Parsed submission dict, or NULL on failure.
coreomics_get_submission <- function(submission_id) {
  if (!coreomics_enabled() || !nzchar(submission_id %||% "")) return(NULL)
  coreomics_get(sprintf("submissions/%s/", submission_id))
}

#' Pull the samples list from a submission's `submission_data.samples`.
#' Each sample is a list with at minimum (`unique_id`, `sample_name`,
#' `condition_name`); additional fields if the submission template
#' included them.
#'
#' @param submission_id Coreomics submission ID
#' @return List of sample dicts (possibly empty), or NULL if submission
#'   not found.
coreomics_get_samples <- function(submission_id) {
  sub <- coreomics_get_submission(submission_id)
  if (is.null(sub)) return(NULL)
  sd <- sub$submission_data %||% list()
  samples <- sd$samples
  if (is.null(samples)) return(list())
  samples
}

#' Heuristic: given a raw file basename, find the coreomics submission +
#' sample that most likely produced it. Uses fuzzy prefix match on
#' `sample_name` (separators stripped).
#'
#' This is intended for LIVE lookups from the UI (Assign Groups & Run,
#' History detail panel, etc.) — for bulk archive matching, prefer the
#' PG side (see scripts/import_coreomics_api.py --auto-link-raws).
#'
#' @param raw_basename Basename of the raw file, e.g. "MCP1.d" or "20260522_293F-WT_Control_R01.d"
#' @param submissions Optional pre-fetched submission list (avoids re-pulling
#'   when called multiple times in a loop). If NULL, fetches the full list.
#' @return List with named elements: $submission_id, $unique_id, $sample_name,
#'   $condition_name, $organism, $description. NULL if no match.
coreomics_match_raw_file <- function(raw_basename, submissions = NULL) {
  if (!coreomics_enabled() || !nzchar(raw_basename %||% "")) return(NULL)
  norm <- function(x) gsub("[._-]", "", x, fixed = FALSE)
  rb_norm <- norm(raw_basename)

  if (is.null(submissions)) submissions <- coreomics_list_submissions()
  if (length(submissions) == 0) return(NULL)

  best <- NULL
  best_match_len <- 0L
  for (sub in submissions) {
    sub_id <- sub$id %||% ""
    samples <- (sub$submission_data %||% list())$samples
    if (!is.list(samples) || length(samples) == 0) next
    for (s in samples) {
      sn <- s$sample_name %||% ""
      if (!nzchar(sn)) next
      sn_norm <- norm(sn)
      if (startsWith(rb_norm, sn_norm) || grepl(sn_norm, rb_norm, fixed = TRUE)) {
        L <- nchar(sn_norm)
        if (L > best_match_len) {
          best_match_len <- L
          best <- list(
            submission_id = sub_id,
            internal_id = sub$internal_id,
            unique_id = s$unique_id,
            sample_name = sn,
            condition_name = s$condition_name,
            organism = (sub$submission_data %||% list())$organism,
            description = (sub$submission_data %||% list())$description,
            submitter = paste(trimws(sub$first_name %||% ""),
                              trimws(sub$last_name %||% ""))
          )
        }
      }
    }
  }
  best
}

#' Convenience for the Group-Assignment UI: given a vector of raw
#' basenames (the files currently loaded), return a data.frame with
#' coreomics-derived condition assignments. Suitable for pre-filling
#' the group selector.
#'
#' @param raw_basenames Character vector of basenames
#' @return data.frame with columns: raw_basename, condition_name,
#'   coreomics_unique_id, coreomics_submission_id, sample_name. Rows
#'   with no match get NA in the coreomics columns.
coreomics_conditions_for_raws <- function(raw_basenames) {
  if (!coreomics_enabled() || length(raw_basenames) == 0) {
    return(data.frame(raw_basename = raw_basenames,
                      condition_name = NA_character_,
                      coreomics_unique_id = NA_character_,
                      coreomics_submission_id = NA_character_,
                      sample_name = NA_character_,
                      stringsAsFactors = FALSE))
  }
  # Fetch once, match many
  subs <- coreomics_list_submissions()
  out <- lapply(raw_basenames, function(rb) {
    m <- coreomics_match_raw_file(rb, submissions = subs)
    list(
      raw_basename = rb,
      condition_name = m$condition_name %||% NA_character_,
      coreomics_unique_id = m$unique_id %||% NA_character_,
      coreomics_submission_id = m$submission_id %||% NA_character_,
      sample_name = m$sample_name %||% NA_character_
    )
  })
  do.call(rbind.data.frame, c(out, list(stringsAsFactors = FALSE)))
}
