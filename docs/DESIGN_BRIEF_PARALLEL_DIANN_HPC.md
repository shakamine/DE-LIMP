# ASMS Poster Content Brief — Parallel Multi-Step DIA-NN Search & HPC Cluster Integration in DE-LIMP

> **Purpose:** Source material for designing an ASMS conference poster.
> **For:** Claude (poster design / layout).
> **From:** Brett Phinney, UC Davis Proteomics Core.
> **Audience:** ASMS attendees — mass spectrometrists, core-facility staff, anyone running large DIA searches.
> **Level:** Highlights — headline story + key results + figure ideas.
> **Date:** 2026-05-22
> **Companion briefs:** `DESIGN_BRIEF_DENOVO_PROTEOGENOMICS.md`, `DESIGN_BRIEF_COMPARATOR_MULTIOMICS.md`.

---

## 0. The one-sentence story

DE-LIMP turns a DIA-NN search on an HPC cluster into a **point-and-click job** — generating a 5-step parallelized SLURM pipeline, submitting it over SSH with self-healing queue management, and auto-loading the results — so a biologist never touches a command line or a sbatch script.

---

## 1. Suggested poster title options

- *"DIA-NN at Scale Without the Command Line: Parallel HPC Search Orchestration in DE-LIMP"*
- *"From the Bench to the Cluster: Automated Parallel DIA-NN Search and SLURM Management in a Shiny GUI"*
- *"Self-Healing HPC Proteomics: Parallel DIA-NN Search with Automatic Queue Switching"*

(Authors / affiliation: Brett Phinney et al., UC Davis Proteomics Core / Genome Center — fill in co-authors.)

---

## 2. The problem (poster intro hook)

Large DIA-NN searches (dozens–hundreds of raw files) are slow and resource-hungry. Running them well on a shared HPC cluster requires writing SLURM scripts, chaining job dependencies, splitting work across nodes, respecting per-user CPU limits, monitoring queues, and recovering from preemption — all command-line bioinformatics that a typical proteomics core user can't and shouldn't have to do. **DE-LIMP hides all of it behind a GUI**, while applying the parallelization strategy recommended by DIA-NN's own developer.

---

## 3. What it does — the 5-step parallel pipeline

Instead of one monolithic DIA-NN job, DE-LIMP generates the **5-step parallelized workflow** (per DIA-NN author Vadim Demichev, Discussion #1414), chained as SLURM dependencies:

```
Step 1: Library prediction (from FASTA)          ── single job
Step 2: First-pass per-file quant   ── ARRAY job (1 task / raw file, parallel)
Step 3: Empirical library assembly               ── single job
Step 4: Final per-file quant        ── ARRAY job (1 task / raw file, parallel)
Step 5: Cross-run report + matrices              ── single job
```

The per-file array steps (2 & 4) are **embarrassingly parallel** — N raw files run as N simultaneous array tasks across the cluster, collapsing wall-clock time. The assembly/report steps (1, 3, 5) are single jobs that depend on the array steps completing.

**Submission is done in just 3 SSH/SCP calls** (mkdir + upload-all-scripts + one launcher that chains every `sbatch --dependency=afterok:$PREV`) to avoid tripping the cluster's connection-rate throttling (`MaxStartups`).

### The hard-won correctness details (credibility points for a methods poster)
These are the things that silently corrupt results if you get them wrong — DE-LIMP encodes the fixes:
- **`--quant-ori-names` on all steps** — preserves original filenames in `.quant` files; without it, container bind-mount path differences cause naming mismatches between steps.
- **`--fasta-search`/`--predictor` in Step 1 only** — leaking them into later steps re-digests the FASTA from scratch.
- **Fixed mass accuracy with `--use-quant`** — auto-optimization + quant reuse gives *different* results; the pipeline forces manual mass-accuracy mode.
- **Empirical library is `.parquet`, not `.speclib`** (DIA-NN 2.0+).
- **Quant-file verification before assembly** — array tasks can fail silently (preemption/OOM); Steps 3 & 5 verify all expected `.quant` files exist before running.
- **Step-2 backup before Step 3** — Step 3 overwrites Step-2 quant files; a backup enables smart resume.

---

## 4. The differentiator — self-healing HPC queue management

This is the most novel piece and the best poster story. DE-LIMP doesn't just submit jobs — it **actively manages them against a real shared-cluster scheduler** (UC Davis "Hive"):

### Automatic queue switching
- Jobs start on the **priority partition** (`genome-center-grp/high`), which has a **per-user 64-CPU cap** (`MaxTRESPU` — the binding constraint, *not* the group limit).
- A monitoring observer polls every 15 s. If a job sits **pending >5 min**, or hits `QOSMaxCpuPerUserLimit`/`InvalidQOS`, DE-LIMP **moves the parallel array steps to the preemptible public partition** (`publicgrp/low`, 1000+ idle CPUs) via `scontrol update`.
- Assembly steps (single jobs that can't restart mid-way) **stay on the priority queue**; only the safely-restartable array steps move.
- Moved jobs get **`Requeue=1`** so SLURM auto-restarts them if preempted; `PREEMPTED` is mapped to "queued," not "failed."
- When priority capacity frees up, queued public-partition jobs are **moved back**.

### Smart resource awareness
- A live **traffic-light cluster indicator** (green/yellow/red) polls `sacctmgr`/`squeue`/`sinfo` every 60 s for group + per-user CPU availability and partition idle capacity.
- Search CPU requests are **auto-capped to the per-user SLURM limit** so jobs aren't rejected.
- `select_best_partition()` picks the starting partition based on live availability.

### Failure recovery & dependency repair
- After retrying failed array tasks (which get a *new* job ID), the pipeline **rewrites the downstream step's `--dependency`** to wait on the retry job — otherwise the assembly step starts before retries finish and the quant-verify fails.

---

## 5. The full GUI experience (end-to-end, no command line)

- **Dual backend:** the same flag-builder (`build_diann_flags()`) drives both **HPC (SLURM over SSH)** and **local Docker** execution; backend auto-detected at startup.
- **SSH file browser** for picking raw files / FASTA on the remote cluster (no scp-by-hand).
- **ControlMaster SSH multiplexing** — one reused TCP connection across all calls (dodges `MaxStartups` throttling; macOS-safe short socket paths).
- **Instrument metadata auto-extraction** at file-scan time (timsTOF `.tdf` / Thermo `.raw`) sets the correct m/z range automatically.
- **Live job monitoring** with per-step / per-array-task progress; **NCBI/UniProt FASTA download** with gene-symbol mapping built in.
- **One-click "Load from HPC"** — SCP the finished `report.parquet` back and auto-run the DE pipeline.
- **`search_info.md`** auto-generated per run: all parameters, job IDs, file list, log paths — a reproducible record (queue-switch events are appended too).

---

## 6. Headline numbers / impact (for large callouts)

- **5-step** parallel pipeline, **N raw files → N simultaneous** array tasks.
- Submitted in **3 SSH calls** (throttle-safe).
- **2 partitions**, automatic switching; jobs survive preemption via auto-requeue.
- **64-CPU per-user limit** auto-respected; **1000+ idle CPUs** tapped on the public queue when priority is full.
- **9 SLURM command paths** proxied so the whole thing works even from inside an Apptainer container with no scheduler access.
- Zero command-line steps for the end user.

> (If you have a concrete wall-clock comparison — e.g. "X files: serial Y h vs parallel Z h" — that single before/after bar would be the strongest possible figure. I can help compute it from a real run's `sacct` timing if you point me at a job ID.)

---

## 7. Suggested figures

- **The 5-step pipeline schematic** — boxes for Steps 1–5, with Steps 2 & 4 shown fanning out into parallel array tasks across nodes; dependency arrows between steps. This is the centerpiece.
- **Queue-switching state diagram** — priority queue → (pending >5 min) → preemptible queue → (capacity returns) → back. The "self-healing" loop is a great visual.
- **Wall-clock before/after bar** (serial vs parallel) — if a real timing is available, this sells the whole poster.
- **Live cluster traffic-light + GUI screenshot** — communicates "no command line."
- **Architecture diagram** — laptop (Docker/GUI) ↔ SSH ↔ HPC (SLURM), showing the dual backend.

---

## 8. Unifying message

DE-LIMP makes **production-scale, correctly-parallelized DIA-NN search on a shared HPC cluster** a GUI operation — encoding the DIA-NN developer's own multi-step recipe *and* the real-world cluster-management logic (per-user limits, preemption, queue switching, dependency repair) that normally requires a bioinformatician. **The biologist clicks "search"; the platform handles the cluster.**

---

## 9. Platform / "by the numbers" sidebar (optional, shared across the poster set)

- Open-source, single-developer Shiny platform: **700 commits / ~16 weeks**, **~47K lines of R**, 39 releases.
- Deploys local / Docker / HPC (Apptainer + SLURM proxy); GitHub + Hugging Face.
- UC Davis Aggie Blue/Gold poster figures already made (`~/Downloads/DE-LIMP_commits_timeline.png`, `DE-LIMP_stat_tiles.png`).

---

## 10. Practical notes for the poster designer (Claude)

- **Lead with the 5-step parallel schematic and the queue-switching loop** — they carry the whole story visually.
- **Punchy numbers to enlarge:** 5 steps · N→N parallel array · 3 SSH calls · 2 auto-switched partitions · 64-CPU cap respected · 9 SLURM paths proxied · 0 command-line steps.
- The "self-healing / survives preemption" framing is the most distinctive angle vs. other GUI search tools — make it prominent.
- UC Davis palette: **Aggie Blue `#022851`**, **Aggie Gold `#FFBF00`**.
- Caption SLURM/HPC jargon (array job, preemption, QOS, partition) — much of the ASMS audience is MS-side, not HPC-side.
- I can generate a real wall-clock parallel-vs-serial figure or a pipeline schematic on request — just point me at a finished job's output dir or `sacct` job ID.
