# DE-LIMP History DB — Design v2 (target v3.11.0)

> **Status:** Draft, 2026-05-18. Supersedes `HISTORY_DB_DESIGN.md` (v1,
> 2026-05-06) which targeted local SQLite. The v1 design is preserved for
> historical reference but should not be implemented as written.
>
> **Why a v2:** between v1 and now we (a) committed to PG Farm rather than
> SQLite, (b) decided the database should serve as a community AI training
> corpus in addition to internal history, (c) discovered HUPO PSI standards
> we need to bake in from day one (USI, ProForma 2.0, SDRF-Proteomics).
> These together change enough of the design that a fresh write is cleaner
> than incremental edits.

---

## 1. Vision

DE-LIMP's history database is **two things at once**:

1. **An operational tool** for the UC Davis Proteomics Core — staff can find
   every search ever run, know who it was for, see what state it's in, and
   load the results into DE-LIMP without hunting through `/quobyte/` and
   `/nfs/` for paths that may have moved.

2. **A community AI training corpus** — every confidently-identified
   peptide, with its modifications / charge state / RT / IM / spectrum,
   accumulated across thousands of searches, indexed for fast lookup. Future
   de novo / spectral-library / RT-prediction model training can pull labels
   directly from this corpus instead of chasing raw files across cluster
   storage.

These two needs share most of the same data. The difference is at the
**access boundary**: the operational tool wants customer names, project
codes, free-text notes, sample sheets. The community corpus wants none of
that — just the biophysical measurements. So the design is a **two-tier
schema**:

- **Public layer** — searchable, exportable, AI-ready, contains only data
  that's safe to share. Adopts PSI standards (USI, ProForma 2.0,
  SDRF-Proteomics) so the corpus is interoperable with ProteomeXchange and
  community ML tooling (MassIVE-KB, ProteomicsML).
- **Internal layer** — joined to the public layer by ID, contains customer
  / project / staff-notes data. Access-controlled at the PG role level so
  the export pipeline physically cannot include it.

## 2. Non-goals (what this is NOT)

These will save grief later if we're clear about them now:

- **Not a raw spectrum archive.** Top-150 peaks per confident-ID MS2 go in
  PG. The full raw .d/.raw files stay on cluster storage (or in PRIDE for
  published data). Storing full raw spectra would blow the storage budget
  by ~100x.
- **Not a real-time ingest pipeline.** Searches are ingested after they
  complete, from `report.parquet`. No coupling to a running DIA-NN process.
- **Not the source of truth for raw acquisition metadata** — that lives on
  the instrument PCs and in STAN. DE-LIMP just stores a cached snapshot at
  ingest time.
- **Not a replacement for ProteomeXchange.** PRIDE / MassIVE are the right
  homes for published datasets. The DE-LIMP DB is the **internal corpus**,
  with a sanitized public-export pipeline for releasing curated subsets.
- **Not a sample LIMS.** Sample sheets, customer billing, and consent
  management belong elsewhere (existing core facility systems). The
  internal layer can REFERENCE these but shouldn't try to own them.

## 3. Architecture

```
                            DE-LIMP App (R/Shiny)
                                    |
                                    | psycopg/RPostgres
                                    v
                ┌─────────── PG Farm (pgfarm.library.ucdavis.edu) ───────────┐
                |                                                            |
                |   Database: uc-davis-genome-center-proteomics-core/delimp  |
                |                                                            |
                |   Public layer (role: delimp_public, SELECT-only):         |
                |   ├─ delimp_searches              (one row per search)     |
                |   ├─ raw_files                    (per-raw metadata)       |
                |   ├─ search_raw_files             (junction)               |
                |   ├─ delimp_sample_metadata       (SDRF-Proteomics-style)  |
                |   ├─ delimp_proteins              (per-search per-protein) |
                |   ├─ delimp_precursors            (per-search per-peptide  |
                |   |                                with top-150 peaks)     |
                |   ├─ delimp_cohorts               (AI training cohorts)    |
                |   ├─ delimp_cohort_members        (junction)               |
                |   ├─ delimp_consensus_ids         (cross-engine consensus) |
                |   └─ delimp_schema_version        (migration history)      |
                |                                                            |
                |   Internal layer (role: delimp_internal, SELECT/INSERT):   |
                |   ├─ delimp_searches_internal     (customer, project,      |
                |   |                                staff notes)            |
                |   └─ delimp_raw_files_internal    (sample labels, custom   |
                |                                    FASTAs, anonymization) |
                └────────────────────────────────────────────────────────────┘
                                    ^
                                    |
                                    | Discovery walker (every N min)
                                    |
                  Cluster storage: /quobyte/, /nfs/, /Volumes/
                  ├─ search output dirs (each has search_info.md)
                  └─ raw files (.d, .raw)
```

**Connection pattern** mirrors STAN's `db_pg.py`:

- Token at `/Volumes/proteomics-grp/brett/.pgfarm_delimp_token` (Mac) and
  `/quobyte/proteomics-grp/brett/.pgfarm_delimp_token` (HPC). Note:
  **separate token file from STAN's** so rotation of one doesn't accidentally
  clobber the other.
- 7-day CAS bearer; refresh weekly via `pgfarm auth login`. The existing
  STAN Wednesday-9am-PT reminder cron should be extended (or duplicated)
  for the DE-LIMP token.
- `sslmode=require` (not `verify-full` — Mac cert path broken; see STAN doc)
- **Module-level cached connection** in R (analogous to STAN's `_connect()`).
  Without it, each insert is ~3.5 s of SSL handshake — bulk ingest of a
  282-file search at 30k precursors each becomes hours instead of minutes.
- `DELIMP_DB_BACKEND=pg` env var to enable. Off → falls back to current CSV
  activity log. Gradual rollout: dev → core staff → all.

## 4. Schema — public layer

All DDL below is intended to be runnable against a fresh `delimp` database.
A migration script (`scripts/migrate_pg_v1.sql`) will execute these in
order with `delimp_schema_version` rows appended.

### 4.1 `delimp_schema_version`

Tracks every migration. The very first thing that runs.

```sql
CREATE TABLE delimp_schema_version (
    version TEXT PRIMARY KEY,                -- semver: '1.0.0'
    applied_at TIMESTAMPTZ DEFAULT NOW(),
    migration_script TEXT,                   -- path to the .sql file
    breaking_change BOOLEAN DEFAULT FALSE,
    notes TEXT
);

INSERT INTO delimp_schema_version (version, migration_script, notes)
VALUES ('1.0.0', 'scripts/migrate_pg_v1.sql',
        'Initial schema — public + internal layers, PSI standards baked in');
```

### 4.2 `delimp_searches`

One row per DE-LIMP search run. Public columns only — customer / project
data lives in `delimp_searches_internal`.

```sql
CREATE TABLE delimp_searches (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),

    -- Identifying
    search_name TEXT NOT NULL,               -- 'Taha_Big_Dog_VER216' (anonymized form for public)
    output_dir TEXT NOT NULL,                -- '/nfs/lssc0/.../...20260518'
    completed_at TIMESTAMPTZ,
    submitted_at TIMESTAMPTZ NOT NULL,

    -- Lineage / reproducibility
    delimp_version TEXT,                     -- '3.10.38' — the app version that ran it
    resubmit_of_search_id UUID REFERENCES delimp_searches(id),
    parent_chain_depth INT DEFAULT 0,        -- 0 = original, 1 = first re-run, etc.

    -- Cached speclib reuse: when s1 (library prediction) was skipped
    -- because a previously-built predicted speclib matched. The link points
    -- back to the search that built it — so the Log Viewer can resolve
    -- "show me the s1 log" even when s1 didn't run here.
    speclib_origin_search_id UUID REFERENCES delimp_searches(id),
    speclib_origin_speclib_md5 TEXT,         -- hash of the .predicted.speclib file
    speclib_origin_built_at TIMESTAMPTZ,     -- when the origin's s1 finished

    -- Search engine + pipeline
    search_engine TEXT NOT NULL DEFAULT 'diann'    -- 'diann','sage','spectronaut','fragpipe'
        CHECK (search_engine IN ('diann','sage','spectronaut','fragpipe','other')),
    search_engine_version TEXT,
    pipeline_id TEXT NOT NULL,               -- 'dpc_quant_limpa', 'maxlfq_limma'
    pipeline_version TEXT,

    -- Search parameters (raw — JSONB so we don't have to schema-evolve every flag)
    search_params_json JSONB,                -- full DIA-NN flag set as submitted
    fasta_path TEXT,                         -- primary FASTA used
    fasta_md5 TEXT,                          -- hash for reproducibility
    fasta_n_proteins INT,
    contaminant_lib TEXT,                    -- 'universal', 'crap', NULL
    contaminant_lib_version TEXT,

    -- Outcome
    n_raw_files INT NOT NULL,
    n_precursors_total INT,                  -- denormalised count, refreshed at ingest
    n_proteins_total INT,                    -- denormalised
    status TEXT NOT NULL                     -- 'queued','running','completed','failed','cancelled'
        CHECK (status IN ('queued','running','completed','failed','cancelled')),
    failure_reason TEXT,

    -- Public-facing FAIR / sharing metadata
    sharing_status TEXT NOT NULL DEFAULT 'private'
        CHECK (sharing_status IN ('private','collaborator','public_pending','public')),
    pride_accession TEXT,                    -- PXD123456 if published
    doi TEXT,                                -- paper DOI if cited
    license TEXT,                            -- 'CC-BY-4.0', etc. NULL for private
    embargo_until DATE,
    citation TEXT,                           -- 'Phinney et al. 2026, doi:...'

    -- Schema versioning for this row
    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version),
    ingested_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_searches_completed_at ON delimp_searches(completed_at DESC);
CREATE INDEX idx_searches_search_name_trgm
    ON delimp_searches USING gin (search_name gin_trgm_ops);
CREATE INDEX idx_searches_pipeline ON delimp_searches(pipeline_id, search_engine);
CREATE INDEX idx_searches_sharing ON delimp_searches(sharing_status);
CREATE INDEX idx_searches_pride ON delimp_searches(pride_accession)
    WHERE pride_accession IS NOT NULL;
```

### 4.3 `raw_files`

Normalised per-raw metadata. Many searches share the same raws; we want
instrument/acquisition info stored once.

```sql
CREATE TABLE raw_files (
    raw_path TEXT PRIMARY KEY,               -- '/quobyte/.../sample.d'
    raw_basename TEXT NOT NULL,              -- 'sample.d' (for filename queries)
    raw_name_anonymized TEXT,                -- 'run_a1b2c3d4.d' for public export

    -- Platform classification (KEY axis for AI training cohorts)
    platform TEXT NOT NULL                   -- 'orbitrap','timstof','tof','other'
        CHECK (platform IN ('orbitrap','timstof','tof','other')),

    -- Instrument (PSI CV terms where possible)
    instrument_model TEXT,                   -- 'timsTOF HT', 'Orbitrap Exploris 480'
    instrument_serial TEXT,
    instrument_cv_accession TEXT,            -- 'MS:1003123' if known
    instrument_cv_name TEXT,                 -- canonical PSI CV term

    -- Acquisition
    acquisition_method TEXT,                 -- 'DDA','DIA','diaPASEF','ddaPASEF'
    acquisition_date TIMESTAMPTZ,
    gradient_minutes REAL,
    lc_method TEXT,                          -- 'EvoSep 60SPD', 'Bruker nanoElute 22min', etc.
    samples_per_day REAL,                    -- SPD if applicable
    sample_amount_ng REAL,

    -- Orbitrap-specific (NULL for timsTOF)
    ms1_resolution INT,                      -- 60000, 120000, etc.
    ms2_resolution INT,
    agc_target INT,
    max_inject_time_ms REAL,
    activation_method TEXT,                  -- 'HCD','CID','EThcD' — affects fragment patterns
    nce REAL,                                -- normalised collision energy (Prosit input)

    -- timsTOF-specific (NULL for Orbitrap)
    mobility_min REAL,                       -- 1/K0 range
    mobility_max REAL,
    n_ms1_frames INT,
    n_ms2_frames INT,
    cycle_time_sec REAL,

    -- Common
    mass_range_min REAL,
    mass_range_max REAL,

    -- File integrity
    file_size_bytes BIGINT,
    md5 TEXT,                                -- expensive to compute; populate lazily

    -- Pointers to derived artifacts (XIC parquets, MGF archives — see § 6)
    xic_parquet_path TEXT,
    xic_parquet_md5 TEXT,
    labeled_mgf_path TEXT,                   -- /quobyte/.../search_X/labeled_spectra.mgf.zst
    labeled_mgf_md5 TEXT,
    labeled_mgf_n_spectra INT,
    labeled_mgf_q_cutoff REAL,

    -- TIC + raw metadata extracted at ingest
    tic_metrics_json JSONB,                  -- AUC, peak RT, gradient width, etc.
    instrument_metadata_json JSONB,          -- raw HyStarMetadata.xml dump, Thermo header

    -- Schema versioning
    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version),
    first_seen_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_raw_files_platform ON raw_files(platform);
CREATE INDEX idx_raw_files_instrument ON raw_files(instrument_model);
CREATE INDEX idx_raw_files_acquisition ON raw_files(acquisition_method, platform);
CREATE INDEX idx_raw_files_date ON raw_files(acquisition_date DESC);
```

### 4.4 `search_raw_files` (junction)

```sql
CREATE TABLE search_raw_files (
    search_id UUID NOT NULL REFERENCES delimp_searches(id) ON DELETE CASCADE,
    raw_path TEXT NOT NULL REFERENCES raw_files(raw_path),
    -- Per-search-per-raw stats (cached for fast filter)
    n_precursors INT,                        -- from report.stats.tsv
    n_proteins INT,
    PRIMARY KEY (search_id, raw_path)
);

CREATE INDEX idx_srf_raw ON search_raw_files(raw_path);
```

### 4.5 `delimp_sample_metadata` (SDRF-Proteomics-compatible)

One row per raw file with ontology-pinned sample-level metadata. This is
the data that PRIDE / ProteomeXchange would want in an SDRF-Proteomics TSV.
Storing it natively means we can emit SDRF on demand.

```sql
CREATE TABLE delimp_sample_metadata (
    raw_path TEXT PRIMARY KEY REFERENCES raw_files(raw_path),

    -- Sample type tagging (KEY for AI training cohort building — see § 7)
    sample_type TEXT                         -- 'study_sample','hela_qc','lysate_std',
                                             -- 'plasma_std','blank','wash',
                                             -- 'standard_other','unknown'
        CHECK (sample_type IN ('study_sample','hela_qc','lysate_std','plasma_std',
                                'blank','wash','standard_other','unknown')),

    -- SDRF-Proteomics fields with ontology accessions
    organism_taxon_id INT,                   -- NCBI Taxonomy ID (9606 = H. sapiens)
    organism_name TEXT,                      -- 'Homo sapiens'
    tissue_efo_accession TEXT,               -- 'EFO_0001185' (kidney), etc.
    tissue_name TEXT,
    cell_line_clo_accession TEXT,            -- 'CL_0000540' (neuron), optional
    cell_line_name TEXT,
    disease_doid_accession TEXT,             -- 'DOID:14330' (Parkinson's), etc.
    disease_name TEXT,

    -- Predicted organism (Phase 4 — populated by the species scanner § 4.11)
    -- These are ADVISORY; organism_taxon_id above remains the authoritative
    -- value when set manually or from SDRF. The predicted columns let the
    -- scanner backfill the species for archived data where organism is
    -- unknown / uncertain.
    predicted_organism_taxon_id INT,
    predicted_organism_name TEXT,
    predicted_organism_confidence REAL,      -- 0.0-1.0, fraction of unique peptides supporting top species
    predicted_organism_method TEXT
        CHECK (predicted_organism_method IS NULL
               OR predicted_organism_method IN (
                   'diamond_blast_diann_peps',
                   'diamond_blast_casanovo_denovo',
                   'manual_override')),
    predicted_organism_top3_json JSONB,      -- runner-up species + counts for sanity review
    predicted_organism_n_peptides_scored INT,
    predicted_organism_at TIMESTAMPTZ,

    -- Experimental design
    biological_replicate INT,
    technical_replicate INT,
    fraction INT,
    label_type TEXT,                         -- 'label_free','TMT16','SILAC_heavy', etc.
    enrichment TEXT,                         -- 'phospho','glyco','plasma_immunodepleted','none'

    -- SDRF emission helper
    sdrf_row_json JSONB,                     -- Pre-built SDRF row for cheap export

    -- Free-text fallback (things SDRF doesn't cover)
    custom_metadata_json JSONB,

    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version)
);

CREATE INDEX idx_sample_type ON delimp_sample_metadata(sample_type);
CREATE INDEX idx_sample_organism ON delimp_sample_metadata(organism_taxon_id);
CREATE INDEX idx_sample_tissue ON delimp_sample_metadata(tissue_efo_accession);
```

### 4.6 `delimp_proteins`

Per-search per-protein-group summary. Mirrors DIA-NN's `report.pg_matrix.tsv`.

```sql
CREATE TABLE delimp_proteins (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    search_id UUID NOT NULL REFERENCES delimp_searches(id) ON DELETE CASCADE,
    raw_path TEXT NOT NULL,                  -- which run; denormalised for query speed
    protein_group TEXT NOT NULL,             -- 'P12345', 'P12345;P12346' for groups
    gene TEXT,
    n_unique_peptides INT,
    n_precursors INT,
    intensity DOUBLE PRECISION,
    normalized_intensity DOUBLE PRECISION,
    pg_q_value REAL,
    is_contaminant BOOLEAN DEFAULT FALSE,    -- Cont_ prefix or in contaminant FASTA
    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version)
);

CREATE INDEX idx_proteins_search ON delimp_proteins(search_id);
CREATE INDEX idx_proteins_group ON delimp_proteins(protein_group);
CREATE INDEX idx_proteins_gene ON delimp_proteins(gene) WHERE gene IS NOT NULL;
```

### 4.7 `delimp_precursors` (the heart of the AI training corpus)

One row per (search × peptide × charge × raw file). For 30k precursors × 2000
searches/year, this grows ~60M rows/year. Postgres handles this fine if
indexes are right.

```sql
CREATE TABLE delimp_precursors (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    search_id UUID NOT NULL REFERENCES delimp_searches(id) ON DELETE CASCADE,
    raw_path TEXT NOT NULL,                  -- denormalised — common filter

    -- Identity (peptide + mods + charge)
    stripped_seq TEXT NOT NULL,              -- 'LLPGFMCQGGDFTR' (no mods)
    modified_seq_diann TEXT,                 -- 'LLPGFMC+57.021QGGDFTR' (DIA-NN inline)
    modified_seq_proforma TEXT,              -- 'LLPGFMC[Carbamidomethyl]QGGDFTR' (PSI)
    charge SMALLINT NOT NULL,
    precursor_id_diann TEXT,                 -- DIA-NN's Precursor.Id
    usi TEXT,                                -- Universal Spectrum Identifier (PSI standard)

    -- Modification breakdown (queryable)
    mods JSONB,                              -- [{"accession":"UNIMOD:35","name":"Oxidation",
                                             --   "aa":"M","pos":7}]
    n_mods INT DEFAULT 0,

    -- Mass / position
    precursor_mz REAL NOT NULL,              -- observed m/z
    predicted_mz REAL,                       -- theoretical m/z from peptide+mods
    precursor_mass DOUBLE PRECISION,         -- neutral monoisotopic mass

    -- Spectrum coordinates (so training can extract spectrum from raw if needed)
    rt REAL,                                 -- retention time, minutes
    im REAL,                                 -- 1/K0 (NULL for Orbitrap)
    ms1_apex_scan INT,                       -- frame/scan number
    ms2_apex_scan INT,

    -- Quality / confidence (essential for cohort filtering)
    q_value REAL,                            -- precursor-level FDR
    global_q_value REAL,                     -- global FDR
    pg_q_value REAL,                         -- protein-group FDR
    pep REAL,                                -- posterior error probability
    empirical_quality REAL,                  -- DIA-NN confidence
    site_localization_probability REAL,      -- for phospho/PTM (per InstaNovo-P)

    -- Quantification
    intensity DOUBLE PRECISION,              -- raw
    normalized_intensity DOUBLE PRECISION,   -- normalised
    intensity_log2 REAL,                     -- log2(intensity), pre-computed for speed

    -- Peak shape summary (cheap, ~30 B; for QC / library quality)
    peak_fwhm REAL,
    peak_asymmetry REAL,
    peak_n_points INT,
    peak_snr REAL,

    -- Spectrum payload (top-150 peaks — the AI training core)
    peak_mz DOUBLE PRECISION[],              -- ordered by intensity descending
    peak_intensity DOUBLE PRECISION[],
    peak_annotation TEXT[],                  -- 'y4+','b3+','' (unannotated) — same length
    n_peaks_total INT,                       -- before top-150 filter
    ms2_spectrum_md5 TEXT,                   -- so cohorts hash-pin exact spectrum content

    -- Cross-engine consensus (denormalised — populated from delimp_consensus_ids)
    n_engines_confirming INT DEFAULT 1,      -- 1 if only this search engine ID'd it

    -- Library context
    library_match TEXT                       -- 'predicted','empirical','dda_confirmed'
        CHECK (library_match IN ('predicted','empirical','dda_confirmed','unknown')),

    -- Supersession (for re-classification without losing the original row)
    superseded_by UUID REFERENCES delimp_precursors(id),
    superseded_reason TEXT,
    superseded_at TIMESTAMPTZ,

    -- Schema versioning
    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version)
);

-- Indexes (sized carefully — too many slows ingest, too few slows queries)
CREATE INDEX idx_prec_search ON delimp_precursors(search_id);
CREATE INDEX idx_prec_raw_path ON delimp_precursors(raw_path);
CREATE INDEX idx_prec_stripped_charge ON delimp_precursors(stripped_seq, charge);
CREATE INDEX idx_prec_protein_via_proteins ON delimp_precursors(search_id, raw_path);
CREATE INDEX idx_prec_qvalue ON delimp_precursors(q_value) WHERE q_value < 0.05;
CREATE INDEX idx_prec_mods_gin ON delimp_precursors USING gin (mods);
CREATE INDEX idx_prec_stripped_trgm
    ON delimp_precursors USING gin (stripped_seq gin_trgm_ops);
CREATE INDEX idx_prec_modified_trgm
    ON delimp_precursors USING gin (modified_seq_proforma gin_trgm_ops);
```

### 4.8 `delimp_consensus_ids`

For AI training, "peptide X identified by ≥2 search engines independently"
is a higher-confidence label than single-engine. This table tracks that.

```sql
CREATE TABLE delimp_consensus_ids (
    raw_path TEXT NOT NULL REFERENCES raw_files(raw_path),
    stripped_seq TEXT NOT NULL,
    modified_seq_proforma TEXT NOT NULL,
    charge SMALLINT NOT NULL,
    engines TEXT[] NOT NULL,                 -- {'diann','sage','spectronaut'}
    best_q_value REAL,
    n_engines INT NOT NULL,                  -- denormalised count
    last_updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (raw_path, modified_seq_proforma, charge)
);

CREATE INDEX idx_consensus_stripped ON delimp_consensus_ids(stripped_seq);
CREATE INDEX idx_consensus_n_engines ON delimp_consensus_ids(n_engines);
```

### 4.9 `delimp_cohorts` + `delimp_cohort_members`

The reproducibility layer for AI training. When you train Casanovo Track H,
the cohort row captures *exactly* which precursor rows it was trained on,
hashed so future you can verify.

```sql
CREATE TABLE delimp_cohorts (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL UNIQUE,                       -- 'Casanovo_TrackH_train_v1'
    description TEXT,
    sql_query TEXT NOT NULL,                          -- WHERE clause used to build it
    parent_cohort_id UUID REFERENCES delimp_cohorts(id),  -- for derived cohorts

    created_at TIMESTAMPTZ DEFAULT NOW(),
    created_by TEXT,
    frozen BOOLEAN DEFAULT FALSE,                     -- false = dynamic, true = locked
    frozen_at TIMESTAMPTZ,
    n_rows INT,
    cohort_md5 TEXT,                                  -- hash of sorted member IDs

    -- Publication / citation metadata
    license TEXT,                                     -- 'CC-BY-4.0'
    citation TEXT,
    doi TEXT,                                         -- Zenodo / DataCite DOI on release
    benchmark_category TEXT                           -- 'in_distribution', 'cross_instrument',
        CHECK (benchmark_category IN (                -- 'cross_species', 'cross_modality',
            'in_distribution','cross_instrument',     -- 'difficult_pep', 'rare_mod', 'general'
            'cross_species','cross_modality',
            'difficult_pep','rare_mod','general')),
    intended_use TEXT                                 -- 'training','validation','test',
        CHECK (intended_use IN (                      -- 'comparison'
            'training','validation','test','comparison')),

    notes TEXT,
    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version)
);

CREATE TABLE delimp_cohort_members (
    cohort_id UUID NOT NULL REFERENCES delimp_cohorts(id) ON DELETE CASCADE,
    precursor_id UUID NOT NULL REFERENCES delimp_precursors(id),
    PRIMARY KEY (cohort_id, precursor_id)
);

CREATE INDEX idx_cohort_members_precursor ON delimp_cohort_members(precursor_id);
```

### 4.10 `delimp_search_step_logs`

DIA-NN emits per-step `.out` and `.err` files (`<output_dir>/logs/diann_s<N>_<step_label>_<jobid>.{out,err}`).
Currently these live only on cluster storage and disappear when the search
dir gets archived — at which point we lose the record of what library was
predicted, which params were used, and why a step failed. Storing log
content directly in PG gives us:

- **Survival past disk archival.** Output dirs get cleaned up; PG rows don't.
- **Cross-search queryability.** *"Find all s1 runs that took >30 min last year"* or *"all s2 array tasks with an OOM message"* becomes pure SQL.
- **Origin-log resolution for cached speclib reuse** (the v3.10.38 Log Viewer feature, but DB-resident so it works without SSH and without the original disk path being valid). When the current search has `speclib_origin_search_id`, the s1 log row lives under THAT search's `search_id`.

```sql
CREATE TABLE delimp_search_step_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    search_id UUID NOT NULL REFERENCES delimp_searches(id) ON DELETE CASCADE,
    step TEXT NOT NULL                           -- 'step1','step2','step3','step4','step5'
        CHECK (step IN ('step1','step2','step3','step4','step5')),
    slurm_job_id TEXT,                           -- '14170165' (per-job ID; .out filename derives from this)
    array_task_id INT,                           -- NULL for single jobs; 0..N for array tasks
    log_type TEXT NOT NULL                       -- 'stdout','stderr'
        CHECK (log_type IN ('stdout','stderr')),

    -- Log payload. PG's TOAST mechanism auto-compresses TEXT > 2 KB, so
    -- there's no need for an explicit bytea/lz4 wrapper. Soft cap of ~5 MB
    -- per row (anything bigger usually means a runaway log we don't want
    -- to keep verbatim — store a summary instead).
    log_content TEXT,
    log_size_bytes INT,
    truncated BOOLEAN DEFAULT FALSE,             -- TRUE if we hit the soft cap

    -- SLURM accounting captured at ingest
    runtime_seconds INT,
    exit_code INT,
    state TEXT,                                  -- 'COMPLETED','FAILED','TIMEOUT', etc.
    cpus_used INT,
    mem_max_gb REAL,
    nodelist TEXT,

    captured_at TIMESTAMPTZ DEFAULT NOW(),
    ingested_schema_version TEXT NOT NULL REFERENCES delimp_schema_version(version)
);

-- Indexes optimized for the queries the Log Viewer + future Lab QC tab need:
CREATE INDEX idx_logs_search_step ON delimp_search_step_logs(search_id, step, log_type);
CREATE INDEX idx_logs_array ON delimp_search_step_logs(search_id, step, array_task_id);
CREATE INDEX idx_logs_failed ON delimp_search_step_logs(exit_code) WHERE exit_code != 0;
CREATE INDEX idx_logs_state ON delimp_search_step_logs(state)
    WHERE state IN ('FAILED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL');
-- Full-text search on log content for grep-style queries
CREATE INDEX idx_logs_content_trgm ON delimp_search_step_logs
    USING gin (log_content gin_trgm_ops);
```

### Smart-capture policy

Storing every byte of every log for every array task would balloon the
table. The discovery walker applies a per-step policy at ingest:

| Step | What's stored | Rationale |
|---|---|---|
| s1 libpred | **Full stdout + stderr.** One file per stream, ~50 KB. Always valuable: library composition, FASTA stats, peptide count. | One job, high info density. |
| s3 assembly | **Full stdout + stderr.** ~100 KB. Cross-run quant params. | Same. |
| s5 report | **Full stdout + stderr.** ~30 KB. Final report stats. | Same. |
| s2 firstpass array | **Per-task summary row** (state + runtime + exit + first 100 + last 100 lines). **Full content for tasks where state ∉ COMPLETED.** | 282 successful tasks of essentially identical logs = wasteful; the 1 task that failed is what matters. |
| s4 finalpass array | Same as s2 | Same reasoning. |

Storage estimate at 2000 searches/year, ~282 raws/search avg:
- s1/s3/s5 full logs: ~180 KB × 3 × 2000 = **~1.1 GB/year**
- s2/s4 summary rows: ~3 KB × 282 tasks × 2 steps × 2000 = **~3.4 GB/year**
- Failed/timeout full content: ~50 KB × ~10/search × 2000 = **~1 GB/year**
- **Total ~5.5 GB/year** for logs (negligible against the 2 TB target).

### Querying for the origin s1 log

The hardest case the Log Viewer needs to handle: "this search reused a
cached predicted speclib (s1 was skipped), show me the original s1 log."
With the FK linkage:

```sql
-- Resolve to whichever search actually ran s1 — this one, or the origin
WITH origin AS (
    SELECT COALESCE(speclib_origin_search_id, id) AS resolved_id
    FROM delimp_searches WHERE id = $current_search_id
)
SELECT l.log_content,
       s.search_name AS log_from_search,
       s.completed_at AS log_built_at,
       (l.search_id != $current_search_id) AS is_origin_fallback
FROM delimp_search_step_logs l
JOIN delimp_searches s ON l.search_id = s.id
JOIN origin o ON l.search_id = o.resolved_id
WHERE l.step = 'step1' AND l.log_type = 'stdout'
ORDER BY l.captured_at DESC LIMIT 1;
```

One query handles both "s1 ran here" and "s1 ran on the origin search."
The boolean `is_origin_fallback` tells the UI whether to show the
"showing log from origin search X" hint.

### Migration / backfill from existing logs on disk

A one-time backfill script (`scripts/backfill_logs.py`) walks the discovery
roots, finds historical `<output_dir>/logs/`, applies the same capture
policy, and INSERTs. Idempotent: skips rows where `(search_id, step,
slurm_job_id, log_type)` already exists.

For searches that already happened before `speclib_origin_search_id` was
populated, the backfill cross-references the speclib cache RDS files
(local `~/.delimp_speclib_cache.rds`, shared `/quobyte/.../speclib_cache.rds`)
and matches each cache entry's `output_dir` to a `search_id` to populate
the origin link retroactively.

### 4.11 Species scanner (Phase 4)

The core's archives accumulate raw files where the recorded organism is
wrong, missing, or just guessed from the filename. Manually re-annotating
thousands of historical files isn't viable. But we ALREADY have the
machinery to predict species from peptide content — it's what we built
for the Cascadia / Casanovo paleoproteomics pipeline (see
`R/helpers_denovo.R` on the `feature/cascadia-denovo` branch:
`run_diamond_blast()`, `normalize_cascadia_blast()`). Phase 4 adapts that
infrastructure into a batch scanner that walks the archive, predicts
species for each raw file, and writes results into the
`predicted_organism_*` columns of `delimp_sample_metadata`.

#### How it works

```
For each raw file in the archive:

  Mode A (default, fast — ~5 min/file):
    1. Find associated report.parquet via PG (delimp_searches.raw_path)
    2. Extract stripped_seq for peptides with Q.Value < 0.01
    3. Filter: peptide length >= 9 (tryptic peptides shorter than this
       hit too many species in SwissProt and are non-discriminating)
    4. Filter: exclude contaminants (protein_group starts with 'Cont_')
    5. Write FASTA, DIAMOND blastp against SwissProt
    6. Tally hits per species, weighted by UNIQUE peptide count
    7. Top species → predicted_organism_*

  Mode B (fallback, slow — ~30 min/file):
    For raw files that have no report.parquet:
    1. Run Casanovo (DDA) or Cascadia (diaPASEF) for de novo IDs
    2. Take novel peptides (no DB match)
    3. Same DIAMOND BLAST + voting as Mode A
```

DIAMOND DB target: `/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_sprot`
(same as the existing paleoproteomics pipeline uses).

#### R-side helper signature

```r
species_scan_archive(
  raw_paths = vector,                # files to scan
  mode = "diann_peps",               # "diann_peps" (fast) or "casanovo_denovo" (slow)
  diamond_db = "/quobyte/.../uniprot_sprot",
  min_q_value = 0.01,
  min_peptide_length = 9,            # protect against cross-species short-peptide noise
  drop_contaminants = TRUE,
  write_to_pg = TRUE,                # auto-update delimp_sample_metadata
  overwrite_existing = FALSE         # skip files already scored
)
```

Returns a data frame:

| raw_path | predicted_organism | confidence | top3 | n_peps |
|---|---|---|---|---|
| `.../HeLa_3.d` | Homo sapiens (9606) | 0.92 | {Hs: 0.92, Mm: 0.04, Rn: 0.02} | 2147 |
| `.../mystery_S2-B1.d` | Canis lupus familiaris | 0.78 | {Cf: 0.78, Hs: 0.08, Mm: 0.05} | 1832 |

#### Sample queries

```sql
-- "Archived files where the predicted species disagrees with the recorded one"
SELECT raw_path, organism_name, predicted_organism_name, predicted_organism_confidence
FROM delimp_sample_metadata
WHERE organism_taxon_id IS NOT NULL
  AND predicted_organism_taxon_id IS NOT NULL
  AND organism_taxon_id != predicted_organism_taxon_id
  AND predicted_organism_confidence > 0.7
ORDER BY predicted_organism_confidence DESC;

-- "Archived files with no recorded organism — fill in from prediction (manual review first)"
SELECT raw_path, predicted_organism_name, predicted_organism_confidence,
       predicted_organism_top3_json
FROM delimp_sample_metadata
WHERE organism_taxon_id IS NULL
  AND predicted_organism_confidence > 0.8;
```

#### Caveats

These caused trouble in the original Cascadia work and would bite again here:

1. **Tryptic peptides are cross-species sticky.** A 7-AA tryptic peptide hits ~50 species in SwissProt. Use `min_peptide_length >= 9`; ideally `>= 12` if the corpus is rich enough to afford it. ~60% of human tryptic peptides have ≥1 mouse hit even at 9-mer cutoff.
2. **Human keratin everywhere.** Almost every spectrum gets human-keratin / human-albumin hits regardless of species. Voting must weight by **unique-peptide count per species** (not all-hits), otherwise *Homo sapiens* always wins.
3. **SwissProt is biased to common species.** Rare organisms (paleoproteomics, marine, niche plants) mis-predict to a more-represented relative. The `predicted_organism_top3_json` is meant for manual review — verify the right genus shows up in top-3 even when top-1 is wrong.
4. **Contaminant DB hits are non-informative.** `Cont_`-prefixed protein groups must be filtered before BLAST. Otherwise universal contaminants (BSA, trypsin, keratin) dominate the vote.
5. **organism_taxon_id is authoritative, predicted_* is advisory.** A high-confidence prediction does NOT auto-overwrite the recorded organism — the scanner only writes the `predicted_*` columns. Any promotion to the canonical `organism_*` field is a human decision, ideally captured via `predicted_organism_method = 'manual_override'`.

#### Sizing

The core's archive contains ~5000 files with uncertain species. Mode A at ~5 min/file:
- Serial: ~16 days (impractical)
- SLURM array, 50 concurrent: ~5 hours (tractable)
- Result: ~5000 new `delimp_sample_metadata` rows updated, minimal storage impact

DIAMOND DB itself stays on cluster storage; not in PG. PG only stores the predictions.

### 4.12 Coreomics submission linkage (cache + auto-link)

Adam Schaal's coreomics system is the canonical record of every customer
submission to the core — 4,362 submissions and 9,419 individual samples
as of 2026-05-19. It holds the rich data we'd otherwise duplicate (PI,
submitter, sample prep choices, gradient length, biohazard status, billing
account, NIH user status, free-text description, etc.) at 61 columns per
submission.

Rather than mirror all 61 columns into `delimp_searches_internal`, we
**cache** a small subset in PG (for fast cross-join queries) and **link**
via the coreomics submission ID. Everything else stays in coreomics; we
query / join on demand.

These cache tables live in the **public layer** despite mirroring data
that's customer-facing, because:

- The data already exists in coreomics; we're not adding new sensitive info.
- The cache supports both internal admin queries AND general analyses
  (e.g., "what HeLa samples were submitted in the last year").
- Customer-specific PII (email, phone) is filtered out before caching —
  only names + institutional affiliation persist.

```sql
-- Cache of coreomics submission rows (4,362 rows as of the 2026-05-19 export).
-- Source: scripts/import_coreomics_xlsx.py reads Adam's periodic .xlsx export.
-- Eventually: live API/DB connection (when coreomics exposes one).
CREATE TABLE coreomics_submissions_cache (
    submission_id TEXT PRIMARY KEY,            -- '5a8c37df3bd1'
    internal_id TEXT,                          -- 'PROT_0701'
    type TEXT,                                 -- 'Proteomics', 'Genomics', etc.
    submitted_at TIMESTAMPTZ,
    status TEXT,                               -- 'Submitted', 'Delivered', etc.
    send_date DATE,

    -- Submitter (names + institution kept; phone/email filtered out at import)
    submitter_first_name TEXT,
    submitter_last_name TEXT,
    submitter_email TEXT,                      -- core staff need this for follow-up
    pi_first_name TEXT,
    pi_last_name TEXT,
    pi_email TEXT,
    institute TEXT,

    -- Frequently-queried submission metadata
    num_samples INT,
    organism TEXT,                             -- customer-supplied, free text
    species TEXT,                              -- legacy / alternate field in coreomics
    prot_or_pep TEXT,                          -- 'Intact Proteins', 'Peptides', etc.
    proteomics_type TEXT,
    mass_spec_wanted TEXT,
    sample_prep TEXT,
    gradient_length TEXT,
    dia TEXT,
    tmt TEXT,
    description TEXT,
    other_info TEXT,

    -- Compliance / billing (admin queries need these for monthly reports)
    biohazard TEXT,
    pathogenic TEXT,
    nih_s10_user TEXT,
    is_nih_major_user TEXT,
    transgenic TEXT,
    po_account_number TEXT,

    -- Full original payload for anything we didn't pull into a column
    raw_payload JSONB,

    -- Cache provenance
    imported_at TIMESTAMPTZ DEFAULT NOW(),
    source_export TEXT,                        -- filename of the .xlsx import this came from
    source_export_md5 TEXT
);

CREATE INDEX idx_coreomics_submitter_email
    ON coreomics_submissions_cache(submitter_email);
CREATE INDEX idx_coreomics_submitter_name
    ON coreomics_submissions_cache(submitter_last_name, submitter_first_name);
CREATE INDEX idx_coreomics_pi_email
    ON coreomics_submissions_cache(pi_email);
CREATE INDEX idx_coreomics_internal_id
    ON coreomics_submissions_cache(internal_id);
CREATE INDEX idx_coreomics_submitted_at
    ON coreomics_submissions_cache(submitted_at DESC);

-- Per-sample rows (9,419 as of 2026-05-19). FK to submission.
CREATE TABLE coreomics_samples_cache (
    submission_id TEXT NOT NULL
        REFERENCES coreomics_submissions_cache(submission_id) ON DELETE CASCADE,
    unique_id TEXT NOT NULL,                   -- 'EC2'
    sample_name TEXT,                          -- '20260522_293F-WT_Control_R01' (filename match key)
    condition_name TEXT,                       -- 'CTL', 'Treated', etc.
    amt_to_inject TEXT,
    internal_id TEXT,                          -- coreomics' own per-sample internal ID
    internal_notes TEXT,
    imported_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (submission_id, unique_id)
);

-- Filename-match index for auto-linkage from raw_files
CREATE INDEX idx_coreomics_sample_name_trgm
    ON coreomics_samples_cache USING gin (sample_name gin_trgm_ops);
CREATE INDEX idx_coreomics_sample_condition
    ON coreomics_samples_cache(condition_name);
```

#### Auto-linkage at ingest

When the discovery walker (§9) ingests a new raw file, after parsing
instrument metadata it runs one additional lookup:

```sql
-- Find the coreomics sample row that this raw file most likely belongs to.
-- Match heuristic: raw_basename startsWith(sample_name) OR contains it
-- with normalised separators.
WITH normalized AS (
    SELECT regexp_replace(raw_basename, '[._-]', '', 'g') AS norm_raw FROM raw_files
    WHERE raw_path = $ingested_raw_path
)
SELECT cs.submission_id, cs.unique_id, cs.sample_name, cs.condition_name
FROM coreomics_samples_cache cs, normalized n
WHERE cs.sample_name IS NOT NULL
  AND (regexp_replace(cs.sample_name, '[._-]', '', 'g') = LEFT(n.norm_raw, length(cs.sample_name))
       OR n.norm_raw LIKE '%' || regexp_replace(cs.sample_name, '[._-]', '', 'g') || '%')
ORDER BY length(cs.sample_name) DESC                -- prefer longer/more-specific matches
LIMIT 1;
```

If found, the discovery walker `UPDATE`s `delimp_raw_files_internal` with
the four `coreomics_*` columns. For the ~5,000-file archive backfill
this is a one-time scan; thereafter it happens per new raw file.

**Expected match rate**: ~70–80%. The remaining 20-30% either pre-date
coreomics (old archived files), got renamed manually, or correspond to
QC/test files that were never submitted via the portal.

#### Sample cross-join queries

```sql
-- "Erik Chow's searches in the last 6 months, with coreomics context"
SELECT s.search_name, s.completed_at,
       cs.internal_id AS coreomics_id,
       cs.organism AS recorded_organism,
       cs.description
FROM delimp_searches s
JOIN delimp_searches_internal si ON s.id = si.search_id
JOIN coreomics_submissions_cache cs ON si.coreomics_submission_id = cs.submission_id
WHERE cs.submitter_email = 'enchow@ucdavis.edu'
  AND s.completed_at > NOW() - INTERVAL '6 months'
ORDER BY s.completed_at DESC;

-- "All raw files linked to a known coreomics sample, grouped by lab/PI"
SELECT cs.pi_last_name, cs.institute, COUNT(*) AS n_files
FROM delimp_raw_files_internal rfi
JOIN coreomics_submissions_cache cs ON rfi.coreomics_submission_id = cs.submission_id
WHERE rfi.coreomics_sample_unique_id IS NOT NULL
GROUP BY cs.pi_last_name, cs.institute
ORDER BY n_files DESC;

-- "Archived raws that don't yet have a coreomics linkage — candidates for review"
SELECT rf.raw_path, rf.acquisition_date
FROM raw_files rf
LEFT JOIN delimp_raw_files_internal rfi ON rf.raw_path = rfi.raw_path
WHERE rfi.coreomics_sample_unique_id IS NULL
ORDER BY rf.acquisition_date DESC;
```

#### Import workflow

For the v3.11.x rollout phase:

1. **Initial backfill**: take the 2026-05-19 .xlsx Adam provided, run
   `scripts/import_coreomics_xlsx.py` once. Loads 4,362 submissions +
   9,419 samples into the cache tables.
2. **Auto-link the archive**: one-time pass that runs the auto-linkage
   query against every existing row in `raw_files`. Expected ~70–80%
   match rate.
3. **Periodic refresh — API-first (2026-05-19 onward)**: Adam provided a
   lab-automation token tied to coreomics' REST API. Cron runs
   `scripts/import_coreomics_api.py` (no human in the loop). Token lives
   at `~/.coreomics_token` (Mac, mode 600) and
   `/quobyte/proteomics-grp/brett/.coreomics_token` (HIVE, mode 640
   group-readable so any lab member's automation can read it). Endpoints
   used: `GET /api/submissions/?lab=PROTEOMICS&page=N` for the list, then
   `GET /api/submissions/<id>/` for each detail (samples live at
   `submission_data.samples` inside the detail view, NOT in the list).
   The xlsx importer (`import_coreomics_xlsx.py`) remains as an offline
   fallback for one-shot backfills.

#### Storage

- `coreomics_submissions_cache`: 4,362 rows × ~2 KB (with JSONB payload) ≈ **9 MB**
- `coreomics_samples_cache`: 9,419 rows × ~150 B ≈ **1.5 MB**
- Total: **~11 MB**. Trivial against the 2 TB target.

#### Privacy note

The PSI / FAIR / public-export pipeline (§ 8) must **not** join through
these tables when emitting public data. Even though submitter names are
arguably less sensitive than the raw file paths we already filter, the
defensive rule is: anything tied to a real person's name stays
admin-only. The `delimp_public` SQL role does NOT have SELECT on the
coreomics cache tables — same protection as the internal layer.

## 5. Schema — internal layer

These two tables hold the customer / project / staff-notes data. They are
**never** queried by the public-export pipeline. Access is controlled at
the role level (§ 5.3).

### 5.1 `delimp_searches_internal`

```sql
CREATE TABLE delimp_searches_internal (
    search_id UUID PRIMARY KEY REFERENCES delimp_searches(id) ON DELETE CASCADE,

    -- Customer / project context
    customer_id TEXT,                        -- 'taha', 'cunningham', etc. (internal codes)
    customer_full_name TEXT,                 -- 'Ameer Taha'
    customer_email TEXT,
    customer_affiliation TEXT,               -- 'UCD Nutrition'
    pi_name TEXT,                            -- principal investigator
    pi_lab TEXT,                             -- 'UCD Wood Lab'

    -- Project metadata
    project_code TEXT,                       -- 'AD_Oct_2025_research'
    project_description TEXT,                -- 'AD cohort proteomics year 2'
    grant_or_funding_source TEXT,
    expected_deliverables TEXT,              -- 'volcano + GSEA + Methods'

    -- Sample sheet (the actual experimental design)
    sample_sheet_path TEXT,                  -- internal path to customer sample sheet
    sample_sheet_md5 TEXT,
    experimental_design_notes TEXT,

    -- Internal workflow tracking
    submitted_via TEXT,                      -- 'core-portal','email','in-person'
    customer_submitted_at TIMESTAMPTZ,
    invoiced BOOLEAN DEFAULT FALSE,
    invoice_id TEXT,
    delivery_status TEXT,                    -- 'in-progress','delivered','archived','blocked'
    delivered_at TIMESTAMPTZ,

    -- Core staff notes
    staff_notes TEXT,
    quality_concerns TEXT,                   -- 'sample S2-C11 likely carryover'
    customer_concerns TEXT,

    -- Lineage to source experimental context
    related_experiments JSONB,

    -- Delivered artifacts
    delivered_report_path TEXT,              -- internal Quarto/PDF report
    invoice_pdf_path TEXT,

    -- Linkage to Adam Olshen's coreomics submission system (the canonical
    -- record of every customer submission to the core). All the rich
    -- submitter / PI / sample-prep / billing data lives in coreomics; we
    -- store the FK and a small cache (see § 4.12) rather than duplicating.
    coreomics_submission_id TEXT,            -- '5a8c37df3bd1'  (coreomics short UUID)
    coreomics_internal_id TEXT,              -- 'PROT_0701'    (human-readable)
    coreomics_last_synced_at TIMESTAMPTZ,

    last_updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- No public indexes (only delimp_internal role queries this).
CREATE INDEX idx_internal_customer ON delimp_searches_internal(customer_id);
CREATE INDEX idx_internal_project ON delimp_searches_internal(project_code);
CREATE INDEX idx_internal_status ON delimp_searches_internal(delivery_status);
CREATE INDEX idx_internal_coreomics ON delimp_searches_internal(coreomics_submission_id)
    WHERE coreomics_submission_id IS NOT NULL;
```

### 5.2 `delimp_raw_files_internal`

```sql
CREATE TABLE delimp_raw_files_internal (
    raw_path TEXT PRIMARY KEY REFERENCES raw_files(raw_path),

    -- Linkage
    customer_id TEXT,
    project_code TEXT,

    -- Original sample identifiers (customer-supplied, often PII-adjacent)
    sample_label_original TEXT,              -- 'JohnDoe_KinaseTreated_3wk_S1A1'
    sample_label_in_paper TEXT,              -- 'Control_3wk_rep1' (per publication)

    -- Custom FASTAs (proprietary protein sequences)
    custom_fasta_path TEXT,
    custom_fasta_md5 TEXT,
    custom_fasta_n_proteins INT,
    custom_fasta_notes TEXT,                 -- 'customer proprietary, do not export'

    -- Anonymisation
    anonymized_run_name TEXT,                -- 'run_a1b2c3d4.d' — populated at first review
    anonymized_at TIMESTAMPTZ,
    anonymized_by TEXT,

    -- Linkage to coreomics per-sample row (the (submission_id, unique_id) PK
    -- in coreomics_samples_cache — see § 4.12). Populated by:
    --   (a) auto-match at ingest time via sample_name prefix on raw_basename
    --   (b) manual override during analysis
    coreomics_submission_id TEXT,            -- '5a8c37df3bd1'  (parent submission)
    coreomics_sample_unique_id TEXT,         -- 'EC2'           (per-sample short ID)
    coreomics_sample_name TEXT,              -- '20260522_293F-WT_Control_R01' (raw basename match key)
    coreomics_condition_name TEXT,           -- 'CTL', 'Treated', etc.
    coreomics_last_synced_at TIMESTAMPTZ,

    last_updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_raw_internal_coreomics_sub
    ON delimp_raw_files_internal(coreomics_submission_id)
    WHERE coreomics_submission_id IS NOT NULL;
CREATE INDEX idx_raw_internal_coreomics_sample
    ON delimp_raw_files_internal(coreomics_sample_name)
    WHERE coreomics_sample_name IS NOT NULL;
```

### 5.3 Roles and access control

```sql
-- The PUBLIC role: read-only access to public layer + write to cohorts.
-- The export pipeline runs as this role. Cannot see customer data.
CREATE ROLE delimp_public NOLOGIN;
GRANT SELECT ON
    delimp_searches, raw_files, search_raw_files, delimp_sample_metadata,
    delimp_proteins, delimp_precursors, delimp_consensus_ids,
    delimp_cohorts, delimp_cohort_members, delimp_schema_version
    TO delimp_public;
GRANT INSERT, UPDATE ON delimp_cohorts, delimp_cohort_members TO delimp_public;

-- The INTERNAL role: full access. Used by DE-LIMP staff sessions.
CREATE ROLE delimp_internal NOLOGIN;
GRANT delimp_public TO delimp_internal;                       -- inherits public access
GRANT SELECT, INSERT, UPDATE, DELETE ON
    delimp_searches_internal, delimp_raw_files_internal
    TO delimp_internal;

-- The INGEST role: writes new data from the discovery walker. Subset of internal.
CREATE ROLE delimp_ingest NOLOGIN;
GRANT delimp_internal TO delimp_ingest;
GRANT INSERT, UPDATE ON
    delimp_searches, raw_files, search_raw_files, delimp_sample_metadata,
    delimp_proteins, delimp_precursors, delimp_consensus_ids
    TO delimp_ingest;

-- Specific user accounts get one of these roles. brettsp → delimp_internal initially.
-- Public/external collaborators → delimp_public.
```

## 6. Spectrum storage strategy

Per § 4.7, top-150 peaks per confident MS2 go directly into PG as
`DOUBLE PRECISION[]` arrays. Rationale:

- **Top-150 covers Casanovo's input filter** (depthcharge `filter_intensity(max_num_peaks=150)`).
  Storing more wastes space; storing less loses training data.
- **PG arrays are the right primitive** — `peak_mz` and `peak_intensity` are
  aligned by index. `peak_annotation[i]` gives the fragment ion label (e.g.,
  `'y4+'`) for peak `i`, empty string if unannotated.
- **Hash-pinning** via `ms2_spectrum_md5` lets cohorts assert *"this exact
  spectrum content"* — survives schema drift and content rewrites.

For full-spectrum extraction (>150 peaks, e.g., for Cascadia training which
prefers unfiltered), the raw `.d`/`.raw` is the source. The
`raw_files.labeled_mgf_path` column points to a per-search MGF archive
generated at search completion — see § 8.

### XIC parquets

DIA-NN already emits XIC parquet files per sample at `<output_dir>/_xic/`.
These are large (~100 MB per sample) and stay on shared HPC storage — we
just store the pointer (`raw_files.xic_parquet_path`). The
discovery walker should mark `xic_parquet_path` as NULL if the file moves
or is deleted (re-check on each walk).

## 7. AI training corpus design

### 7.1 Cohort lifecycle

A cohort goes through three states:

```
dynamic ──→ frozen ──→ published
   │           │           │
   │           │           └─ DOI minted, license set, citation populated
   │           │              (sharing_status='public' on member rows required)
   │           │
   │           └─ membership locked, cohort_md5 computed and stored,
   │              members can never change for this cohort
   │
   └─ membership is a live SQL query result; row count may change as
      new data is ingested. Useful for "what does the current corpus
      look like under these filters?"
```

### 7.2 Building a cohort (example)

```sql
-- Define a Casanovo-style timsTOF training cohort (dynamic)
INSERT INTO delimp_cohorts (name, description, sql_query, intended_use, benchmark_category)
VALUES (
    'casanovo_timstof_train_v1',
    'Confident timsTOF ddaPASEF IDs across the lab corpus, training set',
    'sample_type = ''study_sample'' AND platform = ''timstof'' AND ' ||
    'acquisition_method IN (''DDA'',''ddaPASEF'') AND ' ||
    'q_value < 0.01 AND n_engines_confirming >= 1',
    'training',
    'general'
);

-- Materialise membership (this is what `frozen` does)
INSERT INTO delimp_cohort_members (cohort_id, precursor_id)
SELECT (SELECT id FROM delimp_cohorts WHERE name = 'casanovo_timstof_train_v1'),
       p.id
FROM delimp_precursors p
JOIN search_raw_files srf ON p.raw_path = srf.raw_path
JOIN raw_files rf ON p.raw_path = rf.raw_path
JOIN delimp_sample_metadata sm ON p.raw_path = sm.raw_path
WHERE sm.sample_type = 'study_sample'
  AND rf.platform = 'timstof'
  AND rf.acquisition_method IN ('DDA','ddaPASEF')
  AND p.q_value < 0.01;

-- Freeze it
UPDATE delimp_cohorts
SET frozen = TRUE,
    frozen_at = NOW(),
    n_rows = (SELECT COUNT(*) FROM delimp_cohort_members WHERE cohort_id = id),
    cohort_md5 = (SELECT md5(string_agg(precursor_id::text, ',' ORDER BY precursor_id))
                  FROM delimp_cohort_members WHERE cohort_id = id)
WHERE name = 'casanovo_timstof_train_v1';
```

### 7.3 Export to MGF (for Casanovo / Cascadia training)

A helper script (lives at `scripts/export_cohort_mgf.py`) takes a cohort
name and emits annotated MGF compatible with Casanovo:

```python
# Pseudo-code
cur.execute("""
    SELECT p.modified_seq_diann, p.charge, p.precursor_mz, p.rt, p.im,
           p.peak_mz, p.peak_intensity, rf.raw_basename, p.ms2_apex_scan
    FROM delimp_cohort_members m
    JOIN delimp_precursors p ON m.precursor_id = p.id
    JOIN raw_files rf ON p.raw_path = rf.raw_path
    WHERE m.cohort_id = %s
""", [cohort_id])

with open(out_path, 'w') as f:
    for row in cur:
        f.write("BEGIN IONS\n")
        f.write(f"TITLE={row.raw_basename}:scan:{row.ms2_apex_scan}\n")
        f.write(f"PEPMASS={row.precursor_mz}\n")
        f.write(f"CHARGE={row.charge}+\n")
        f.write(f"SEQ={row.modified_seq_diann}\n")
        if row.im is not None:
            f.write(f"ION_MOBILITY={row.im}\n")
        for mz, intensity in zip(row.peak_mz, row.peak_intensity):
            f.write(f"{mz} {intensity}\n")
        f.write("END IONS\n\n")
```

The export is **fully reproducible** from `cohort_md5`: anyone with read
access to the public layer can regenerate the exact MGF.

## 8. Sanitization pipeline (export-time anonymisation)

When a row's `sharing_status` is upgraded from `private` to `public_pending`,
a review process runs:

1. **Anonymise `search_name`**. Original is fine internally but:
   - Strip customer codes (`Taha_*` → `srch_<hash>`)
   - Strip date suffixes if they're sub-day granular
2. **Anonymise `output_dir`**. Don't expose the raw path — replace with
   `<corpus_id>/<sanitised_name>` form.
3. **Anonymise `raw_files.raw_basename`**. Use `anonymized_run_name` if
   populated (populated during review).
4. **Drop custom FASTAs entirely**. Public rows must reference only
   standard UniProt + standard contaminant libraries. Searches against
   proprietary FASTAs aren't shareable as published-search rows.
5. **Round timestamps**. `acquisition_date` → month granularity for public.
6. **Drop `instrument_serial`** from public view (instrument model is fine;
   serial can be triangulated to a specific lab).
7. **Strip free-text from `tic_metrics_json` / `instrument_metadata_json`**
   if those JSONB columns ever happen to contain customer-supplied text
   (they shouldn't, but defensive).

Implementation: this is a per-row review step, not a fully-automatic export.
A `public_pending` row needs a human's sign-off (`reviewed_by` field on a
review table — TBD) before flipping to `public`.

## 9. Discovery walker

A periodic process (cron or systemd timer) that finds completed searches
and ingests them. Walks configured roots, identifies new search dirs by
their `search_info.md` (existing DE-LIMP convention), ingests into PG.

```yaml
# ~/.delimp_discovery.yaml — site-configurable
roots:
  - /quobyte/proteomics-grp/de-limp        # UCD HIVE shared storage
  - /quobyte/proteomics-grp/service        # core facility customer dirs
  - /Volumes/proteomics-grp/de-limp        # Mac-mounted same
  # Add NFS / mounts here as needed
ignore_globs:
  - "*_failed_*"
  - "*archive*"
scan_interval_minutes: 30
batch_size: 10                              # max searches per scan
```

For each new search:
1. Parse `search_info.md` → search params
2. Read `report.parquet` → per-precursor + per-protein data
3. Read `report.stats.tsv` → per-file precursor counts
4. Parse raw file metadata (re-use existing `parse_raw_file_metadata()`)
5. **Walk `<output_dir>/logs/` and ingest into `delimp_search_step_logs`** via the smart-capture policy in §4.10 (full content for single-job steps + failed array tasks; summary-only for healthy array tasks).
6. **Resolve `speclib_origin_search_id`** by cross-referencing the speclib cache entry's `output_dir` to an existing `delimp_searches.id`.
7. Insert into PG within a transaction
8. Emit `labeled_spectra.mgf.zst` archive (extracted from raw files)
9. Update `xic_parquet_path` / `labeled_mgf_path` pointers

Failure modes:
- Search dir disappears mid-ingest → transaction rollback, retry on next scan
- `report.parquet` corrupt → log to `delimp_ingest_failures` table, skip
- Customer's path is `/Volumes/...` (Mac-only) → translate via the same
  `PATH_TRANSLATIONS` map STAN uses; see `stan/cli.py::ingest_orphans_cmd`

## 10. Robust "Load this search" UX

The History tab's "Load" button should never break because a path moved.
Resolution order:

1. **Original `output_dir`** — if accessible, fetch `report.parquet` from there.
2. **PG-canonical path** — `raw_files.xic_parquet_path` etc., which the
   discovery walker keeps current.
3. **PG payload** — if neither path works, reconstitute the analysis from
   PG itself: pull all precursors / proteins for `search_id`, build an
   in-memory `report.parquet`-equivalent, hand to DE-LIMP's loader.

That third path is the "DB-as-source-of-truth" mode and is what makes the
DB a true backup of the analysis (not just an index of where files used to
live).

## 11. Migration from CSV (current activity log)

DE-LIMP currently writes per-event rows to
`~/.delimp_activity_log.csv` (per user, per Mac). Migration plan:

**Phase 0 (now):** CSV continues to be the source of truth.

**Phase 1 (v3.11.0 release):** Dual-write enabled when `DELIMP_DB_BACKEND=pg`.
Every event goes to both CSV and PG. CSV is still the read source for the
History tab. PG is shadowed.

**Phase 2 (v3.11.1):** PG becomes the read source. CSV is dual-written for
back-compat. Existing CSV rows are migrated via
`scripts/migrate_csv_to_pg.py` — one-shot import of historical activity log.

**Phase 3 (v3.12.0):** CSV write disabled. PG-only. CSV files are kept
read-only as a fallback in case PG is down.

**Phase 4 (v3.12.x):** AI cohort registry + sanitization pipeline goes
public. First community releases.

## 12. Phased rollout — concrete deliverables

| Version | Deliverables | ETA |
|---|---|---|
| **v3.11.0** | PG dual-write; schema v1.0 deployed; discovery walker MVP; History tab reads from PG when `DELIMP_DB_BACKEND=pg` | 2-3 weeks after admin grants access |
| **v3.11.1** | PG as primary read; CSV migration script; activity log rows backfilled | +2 weeks |
| **v3.11.2** | Lab Proteome tab (cross-search analytics, cumulative stats) | +2 weeks |
| **v3.12.0** | CSV write disabled; PG only | +4 weeks |
| **v3.11.2 (Phase 2.5)** | Coreomics submission linkage (§4.12): cache tables + import script (`scripts/import_coreomics_xlsx.py`) + auto-link discovery walker step. One-time archive backfill (~5,000 raws → estimated 70-80% match rate). | Concurrent with v3.11.2 |
| **v3.12.x** | Cohort registry UI; sanitisation pipeline; first public cohort release; SDRF-Proteomics emission | TBD |
| **v3.13.x (Phase 4)** | Species scanner (§4.11): batch-predict organism for archive files via DIAMOND blastp on existing DIA-NN peptide IDs; auto-populate `predicted_organism_*` columns. Needs Mode B (Casanovo fallback) only if a substantial fraction of archive lacks search results. | After Phase 3 |

## 13. Risks & open questions

### Risks

1. **PG Farm token expiry every 7 days.** If the refresh is missed,
   everything silently breaks (auth failures). Mitigation: extend STAN's
   existing Wednesday-9am-PT reminder cron to cover DE-LIMP too; consider
   adding a "PG token expires in < 24h" warning banner inside DE-LIMP.
2. **Schema drift over years.** Mitigation: `delimp_schema_version` table
   + every row carries `ingested_schema_version`. Migrations are
   forwards-only and documented.
3. **Discovery walker missing data**. If a search dir doesn't have
   `search_info.md` (older runs from before v3.6.0), it's invisible.
   Mitigation: a backfill script that reconstructs `search_info.md` from
   sbatch scripts in the dir.
4. **Customer data leaks into public exports.** Mitigation: role-based
   access (`delimp_public` cannot see `_internal` tables); export pipeline
   runs as `delimp_public`; every public-bound row passes through review.
5. **Funding lapses.** Community resources die without sustained
   engineering. Mitigation: design for low ongoing maintenance (no daemons,
   minimal moving parts beyond the periodic walker); plan a successor
   before applying for the first grant.

### Open questions

1. **Decoys.** DIA-NN drops decoys before `report.parquet`. Do we want a
   parallel ingest path that pulls them in (with `--decoys` enabled)?
   Useful for FDR re-training. Defer to v3.12.x.
2. **mzTab / mzML export.** Should we emit those formats for PRIDE
   submission? Probably yes, eventually. Defer.
3. **Public read-only API.** REST/GraphQL endpoint? Or just "ask the lab
   for a SQL dump"? Probably PostgREST as a first pass. Defer to v3.12.x.
4. **Review workflow for `public_pending` → `public`.** Manual sign-off
   table? Multi-reviewer approval? Out of scope for v3.11.0.
5. **DOI minting.** Zenodo / DataCite integration. Defer until first
   public cohort release.

## 14. Standards reference (for future audits)

| Standard | Used in | Source |
|---|---|---|
| HUPO PSI **USI** | `delimp_precursors.usi` | [psidev.info/usi](https://www.psidev.info/usi) |
| HUPO PSI **ProForma 2.0** | `modified_seq_proforma`, `mods.accession` | [psidev.info/proforma](https://www.psidev.info/proforma) |
| HUPO PSI **CV** (mass-spec ontology) | `instrument_cv_accession`, `instrument_cv_name` | [psidev.info/mass-spectrometry](https://www.psidev.info/mass-spectrometry) |
| **SDRF-Proteomics** | `delimp_sample_metadata` (whole table mirrors SDRF columns) | [github.com/bigbio/sdrf-proteomics](https://github.com/bigbio/sdrf-proteomics) |
| **UniMod** | `mods.accession` (e.g. `UNIMOD:35` for Oxidation) | [unimod.org](https://www.unimod.org) |
| **NCBI Taxonomy** | `delimp_sample_metadata.organism_taxon_id` | [ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) |
| **EFO** (Experimental Factor Ontology) | `tissue_efo_accession` | [ebi.ac.uk/efo](https://www.ebi.ac.uk/efo/) |
| **Cell Line Ontology** | `cell_line_clo_accession` | [obofoundry.org/ontology/clo](https://obofoundry.org/ontology/clo.html) |
| **Disease Ontology** | `disease_doid_accession` | [disease-ontology.org](https://disease-ontology.org) |
| **ProteomeXchange 2026 (FAIR/AI)** | Architecture motivations | [academic.oup.com/nar/article/54/D1/D459](https://academic.oup.com/nar/article/54/D1/D459/8315797) |
| **MassNet** (billion-scale AI corpus) | Schema inspiration | [biorxiv.org/content/10.1101/2025.06.20.660691v1](https://www.biorxiv.org/content/10.1101/2025.06.20.660691v1.full) |

## 15. Decision log (the "why" for each non-obvious choice)

- **PG Farm, not local SQLite.** v1 design picked SQLite for simplicity.
  Switched because: (a) SQLite at the size we're targeting hits write
  contention with multiple users, (b) UCD Library's PG Farm is free for
  research, (c) STAN already pioneered the pattern there.
- **Same org as STAN, separate DB.** Cross-tool joins via SQL aren't
  needed (raw-file path equality is enough for application-side matching);
  separate DBs let DE-LIMP schema evolve independently.
- **Top-150 peaks, not full spectra.** Top-150 covers Casanovo's input
  filter; storing more is wasted space for the dominant training use case.
  Cascadia needs unfiltered → falls back to per-search MGF archive on disk.
- **`raw_files` separate from `delimp_searches`.** Many searches use the
  same raws (especially HeLa QCs); normalisation saves space and keeps
  instrument metadata canonical.
- **Two-tier (`_internal` vs public).** Earlier draft had `internal_only`
  bool columns; that's error-prone because each query has to remember to
  filter. Separate tables + role-based access make leaks impossible by
  construction.
- **UUID primary keys, append-only writes.** Once a row is published, its
  ID must be citable forever. Auto-increment ints renumber on migration.
  Append-only + `superseded_by` preserves history.
- **JSONB for `mods` and `search_params_json`.** Modifications and DIA-NN
  flag sets evolve faster than the schema should. JSONB indexed via GIN
  gives queryability without ALTER TABLE per change.
- **PSI standards (USI, ProForma, SDRF) baked in from day one.**
  Retrofitting these into a populated DB is painful. Better to mint USIs
  at ingest and emit SDRF on demand than try to backfill years later.
- **`sharing_status` defaults to `private`.** Every row enters private;
  promotion to public is an explicit decision. Inverse default (public-
  unless-flagged) is unsafe.

---

**Next steps after admin provisions the database:**

1. Save credentials to `/Volumes/proteomics-grp/brett/.pgfarm_delimp_token`
   and the parallel HPC path.
2. Run `scripts/migrate_pg_v1.sql` to create the schema.
3. Create R-side connection helper (`R/db_pg.R`) modelled on `stan/db_pg.py`.
4. Wire DE-LIMP's activity log to dual-write when `DELIMP_DB_BACKEND=pg`.
5. Run discovery walker against the existing search history; verify a few
   searches roundtrip cleanly.
6. Cut v3.11.0 release.

---

*End of design doc. Open issues filed in the project tracker under
`v3.11.0-history-db` label.*
