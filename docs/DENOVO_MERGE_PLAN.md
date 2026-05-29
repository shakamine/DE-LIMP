# De Novo Branch — Merge-to-Main Prep Plan

> **Status as of 2026-05-22.** `feature/cascadia-denovo` is NOT yet ready to merge into `main`.
> Do this as one focused session **after the proteogenomics GUI is finished** (so both
> feature branches merge together and the de novo branch doesn't drift behind again).
> Written for whoever (future Claude session or Brett) does the merge.

## Branch state (measured 2026-05-22)

- `feature/cascadia-denovo` is **67 commits ahead**, **148 commits behind** `main`.
- Common ancestor: `3bb8bb7` ("docs: Add nightly documentation GitHub Action to TODO").
- Branch VERSION is **3.7.0**; `main` is **3.10.33** (it will inherit main's version on merge).
- 67 commits of Cascadia/Casanovo/DDA work: ~25K lines across 52 files
  (cascadia/* python, R/server_dda.R, R/server_denovo*.R, R/helpers_denovo.R, R/helpers_dda.R, etc.).

## The 6 conflict files (from `git merge-tree` preview)

Both branches edited these; they need manual conflict resolution:

| File | Why it conflicts | Resolution guidance |
|------|------------------|--------------------|
| `R/ui.R` | cascadia added the De Novo dropdown + 3-mode (DIA/DDA/XL-MS) switcher; main rebuilt large parts of the navbar/UI | Keep main's current navbar structure; re-insert the De Novo dropdown + mode switcher into it. Highest-judgment conflict — review carefully. |
| `R/server_search.R` | both heavily edited (main: Docker/WSL/queue work; cascadia: de novo submission hooks) | Keep main's version as the base; layer cascadia's de novo SLURM submission additions on top. |
| `R/server_ai.R` | both touched AI prompt/export | Take main's base; re-add any de novo-specific AI context. |
| `R/server_qc.R` | both touched QC | Take main's base; re-add de novo QC bits if any. |
| `R/server_gsea.R` | main has the v3.10.33 GSEA org-db fix + other gsea changes; cascadia has older gsea | **Take main's version** (it has the org-db install fix). cascadia has nothing newer here. |
| `CLAUDE.md` | both edited project docs | Merge by hand — keep main's structure, fold in any de novo doc notes. |

`app.R` does NOT conflict (auto-merges) despite both editing it — good.

## Recommended sequence (safe order)

1. **Finish the proteogenomics GUI first** (on `feature/proteogenomics-builder`).
2. **Update the de novo branch, not main:** check out `feature/cascadia-denovo`,
   then `git merge main`. Resolve the 6 conflicts **on the feature branch** so `main` stays clean.
3. **Launch the app and test** — not just "does it compile":
   - App starts, navbar renders, all existing tabs work.
   - De Novo dropdown appears; the 12 de novo sub-tabs render.
   - Ideally run one real Cascadia or Casanovo GPU job end-to-end (de novo features need the cluster to truly verify).
4. **Bump VERSION + CHANGELOG** on the feature branch for the de novo feature set.
5. **Then merge into main** — now a clean, low-risk merge.
6. Do the same for `feature/proteogenomics-builder` (it's already current with main via the
   v3.10.33 merge, so its merge to main should be much smaller).

## Already done (2026-05-22)

- GSEA org-db install fix committed to `main` (v3.10.33) and merged into
  `feature/proteogenomics-builder`. The de novo branch will inherit it automatically
  at merge time (it never touched those gsea lines beyond the old version → main's wins).

## Stashes — IGNORE THEM

- `stash@{1}` "cascadia wip 2" and `stash@{2}` "cascadia wip + v3.8.0 edits" are
  **identical and obsolete** (v3.8.0-era; they delete PATTERNS.md and set VERSION 3.8.0).
  `main` is long past that. Do **not** apply them — they'd revert newer work. They hold no
  unique unsaved work. (Left in place; dropping is destructive and unnecessary.)
- `stash@{0}` "pre-proteogenomics WIP" is a separate `main` stash — assess separately if relevant.

## Do NOT

- Merge cascadia-denovo directly into `main` without updating+testing it first.
- Apply the old stashes.
- Run heavy de novo verification on HPC login nodes — submit via sbatch.
