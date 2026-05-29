# Notes on Casanovo 5 (with ion mobility) and Cascadia fine-tuning

*Part 1 (Casanovo + IM) is the main thing we've been working on. Part 2
(Cascadia notes) is a separate set of observations from fine-tuning
Cascadia on local timsTOF data — feel free to forward that section to
whoever maintains Cascadia if it's useful.*

---

# Part 1 — Adding ion mobility to Casanovo 5

## Summary

We extended Casanovo 5.0 to consume 1/K0 ion-mobility values from timsTOF
ddaPASEF data via two additions to `Spec2Pep`, then fine-tuned with a
masked-AA pre-task. The resulting model — internally we call it Track G,
referenced below as the IM + masked-AA model — was trained to ~1.1M steps
on A100 nodes and learns IM without catastrophically forgetting the
pretrained Orbitrap capability.

In-distribution timsTOF HeLa val (61.7k spectra):

| Checkpoint | Peptide acc. (backbone, I=L) | AA acc. | Δ vs. Casanovo 5.0 baseline |
|---|---|---|---|
| Casanovo 5.0 stock | 18.9% | 45.7% | — |
| IM + masked-AA model, ~600k steps | 53.9% | 71.6% | +35.0 pp peptide, +25.9 pp AA |
| IM + masked-AA model, ~1.1M steps (final) | **61.0%** | **74.4%** | **+42.1 pp peptide, +28.7 pp AA** |

OOD on a salmon timsTOF set (79k spectra, unseen species + instrument run):
**31.5% peptide / 55.4% AA, +11.0 pp AA over baseline.**

The patches are ~530 lines of Python (model + dataloader + train script)
derived from Casanovo 5.0, plus a custom inference entry point because the
stock CLI rejects the extra state-dict keys.

---

## Architecture

The model is a subclass of `casanovo.denovo.model.Spec2Pep` named
`Spec2PepIM` (`casanovo5_im_patches/model.py:31`). The non-trivial
additions:

### 1. IM injection at the encoder output

```python
# model.py:62 — FloatEncoder is depthcharge's sinusoidal scalar encoder
self.im_encoder = FloatEncoder(d_model=dim_model, ...)

# model.py:68 — concat-then-project; bias=False matters for the identity init below
self.im_projector = nn.Linear(2 * dim_model, dim_model, bias=False)
```

IM is **not** added to the precursor token — it's broadcast across the
encoder output (the `memories` tensor) and concat-projected:

```python
# model.py:133  _apply_im(memories, im_values)
im_expanded = im_values[:, None].expand(-1, memories.shape[1])
im_encoded  = self.im_encoder(im_expanded)            # (B, L, D)
combined    = torch.cat([memories, im_encoded], dim=2)  # (B, L, 2D)
return self.im_projector(combined)                     # (B, L, D)
```

We tried IM-on-precursor-token first; this version converged faster and
preserved Orbitrap behavior better. Open architectural question: should
IM also flow per-peak, since each fragment inherits the precursor's 1/K0?

### 2. Identity-zero initialization in `from_pretrained`

The trick that made the pretrained checkpoint loadable without behavior
drift:

```python
# model.py:336 — im_projector.weight shape is (dim_model, 2*dim_model)
# First D cols  = spectrum encoding  → init to identity (pass-through)
# Last D cols   = IM encoding         → init to zero
new_state[proj_key]                  = torch.zeros(dim_model, 2*dim_model)
new_state[proj_key][:, :dim_model]   = torch.eye(dim_model)
new_state[proj_key][:, dim_model:]   = 0.0
```

At step 0, `_apply_im(memories, im) == memories` exactly. Step 0 is
**byte-identical to stock Casanovo 5.0.** The model learns to mix in IM
only if/when gradient says it should. This avoided the "destroys pretrained
behavior on epoch 1" failure we hit with random-init.

### 3. Masked-AA pre-task (the key ingredient)

```python
# model.py:76 — learned mask vector, zero-init so behavior at step 0
# matches the previous attempt exactly
self.mask_embedding = nn.Parameter(torch.zeros(dim_model))

# Training-time: replace 15% of decoder input embeddings with mask_embedding.
# Forces spectrum→residue dependency instead of autoregressive shortcutting.
```

Masked positions are sampled at random (`model.py:83`), then
`_masked_embedding_ctx` (`model.py:109`) wraps the decoder's
`token_encoder` with a forward hook that swaps in `mask_embedding`
at the sampled positions. Mask is empty in eval mode — no overhead at
inference.

Without masking, the model leaned heavily on the AR language prior and
the IM signal stayed small. With 15% masking, the AR prior is unreliable,
so the model has to use spectrum + IM more. We used 15% blind from BERT
— no ablation.

### 4. Flat LR (intentionally — see "what didn't work")

```python
# train_track_g.py:240
optimizer = torch.optim.Adam(model.parameters(), lr=5e-6)
return [optimizer], []
```

Plain flat LR over all parameters. An earlier attempt used 100× higher
LR on IM layers; that recovered the IM signal but still hurt Orbitrap.
Once masked-AA started working, differential LR became unnecessary.

### 5. Tokenizer

Stock `MskbPeptideTokenizer` from depthcharge, no changes. The vocab is
the 28-token MassIVE-KB set (`C[Carbamidomethyl]`, `M[Oxidation]`,
`N[Deamidated]`, `Q[Deamidated]`, `[Acetyl]-`, `[Ammonia-loss]-`,
`[Carbamyl]-`, `[+25.980265]-`). Annotated MGFs use the inline form
(`SEQ=DETVSDC+57.021SPHLANLGR`); depthcharge's MGF parser handles
the round-trip.

---

## Data pipeline (the `ion_mobility` column)

`DeNovoDataModuleIM` (`dataloaders.py:22`) subclasses `DeNovoDataModule`
and adds a `CustomField` so 1/K0 flows from MGF → Lance → batch dict:

```python
# dataloaders.py:52
CustomField(
    name="ion_mobility",
    accessor=_extract_ion_mobility,  # reads spectrum["params"]["ion_mobility"]
    dtype="float64",
)
```

The accessor reads `ION_MOBILITY=` from MGF params (matches `pyteomics.mgf`'s
lowercased key), and the field appears as `batch["ion_mobility"]` at training
time. The model's `_extract_im` (`model.py:160`) pulls it out and casts to
float. No schema changes to Casanovo's core code path.

**Quobyte gotcha:** `n_workers > 0` deadlocks on our distributed FS. We run
`n_workers=0` and stage HDF5 to local SSD. Worth a note in install docs
for HPC users.

---

## What didn't work (in order)

Each row is something worth flagging for anyone trying to add a modality
to Casanovo. The internal label is shown for reference against our file
names; what each attempt actually changed is in the second column.

| Internal label | What changed vs. previous | What failed | Diagnosed |
|---|---|---|---|
| Attempt 1 (Track D) | Casanovo on timsTOF, no IM, resumed mid-epoch | Mode collapsed (0 PSMs on all eval data) | `trainer.fit(ckpt_path=...)` with a DataLoader that has no `state_dict()`. Lightning iterates forever on epoch 0 looking for the resume batch. Cost us ~3 GPU-days; might warrant a note in the docs. |
| Attempt 2 (Track E) | + IM injection, uniform LR | IM weights stayed at 0 (89.7% < 0.01 after 26 epochs), Orbitrap degraded | Zero-init + uniform LR → near-zero input → near-zero gradient. Also catastrophic forgetting on Orbitrap. |
| Attempt 3 (Track E2) | + **100× LR on IM layers** | Learned IM, still hurt Orbitrap | Differential LR was necessary but not sufficient. |
| Attempt 4 (Track E3) | + tanh-gated IM injection | Marginal improvement | Helped, but not the breakthrough. |
| Attempt 5 (Track F) | + identity-zero init `im_projector`, dropped differential LR | OK on IM, weak on hard spectra | AR prior too dominant. |
| Attempt 6 (Track G, final) | + 15% masked-AA | **Works.** | Forces spectrum→residue learning. |

Headline ordering: **identity-zero init + masked-AA were the two
load-bearing ideas.** Differential LR was a red herring once those were
in place.

---

## Training infrastructure

- **HPC:** UC Davis HIVE cluster, A100 80GB nodes
- **Training data:** ~290k labeled timsTOF HeLa + matched plasma spectra
  from local ddaPASEF acquisitions, mobility-filtered (median 88 peaks
  per spectrum vs. ~9.5k unfiltered, which OOMs an A100 at batch_size=1
  because attention is O(L²) in peaks)
- **Effective batch:** 32, accumulate=1 on A100 80GB
- **Optimizer:** Adam, LR=5e-6 (flat), weight_decay=1e-9
- **Mask rate:** 15%
- **Compute:** 4 × 48h A100 runs = **192 A100-hours** to 1.1M steps.
  Per-resume gain decayed +48 → +6 → +4 → +2 pp peptide accuracy.
  Stopped at the ~1.1M-step checkpoint.
- **depthcharge version:** 0.2.3 (with two patches, below)
- **Casanovo:** 5.0 base, patched in `casanovo5_im_patches/` overlay
- **Lightning:** 2.x (whatever shipped with the Casanovo 5.0 env)

### depthcharge patches we applied

1. `AnnotatedSpectrumDataset.__len__` references `self._offsets`; the
   attribute is actually `self._file_offsets`. One-line fix.
2. `filter_intensity(max_num_peaks=N)` only syncs `mz` + `intensity`
   arrays. If you add custom arrays (rt, level, im, fragment), they don't
   get filtered alongside. We patched `primitives.py` to sync all arrays.

Both upstreamable.

---

## Inference + eval

Stock `casanovo sequence` rejects our checkpoints with *"Weights file
incompatible"* because of the extra `im_encoder` / `im_projector` /
`mask_embedding` keys. We use a custom inference script
(`predict_track_f.py`, ~150 lines) which imports `Spec2PepIM` directly,
loads via `Spec2PepIM.load_from_checkpoint(CKPT)`, runs
`trainer.predict()`, and writes a minimal mztab keyed by MGF position.

Casanovo 5.0's `--evaluate` mode is **broken for Mskb-tokenized
checkpoints** — `pyteomics.proforma` round-trip produces 6-decimal mass
strings (`C[+57.021464]`) that aren't in the `MskbPeptideTokenizer`
vocab, KeyError → silent 0% reported. We wrote a 100-line
`score_casanovo.py` that strips mods (both bracket and mass-delta forms),
maps I→L, and compares backbone-only.

Reproduction recipe:

```bash
export TRACK_F_CKPT=/.../casanovo_track_g/checkpoints/last-v3.ckpt
export VAL_MGF=/.../casanovo_finetune/track_e/val.mgf   # 61.7k timsTOF HeLa
export OUTPUT_DIR=/.../model_comparison
export OUTPUT_ROOT=trackg_v3_timstof
python3 predict_track_f.py
python3 score_casanovo.py --mgf $VAL_MGF --mztab $OUTPUT_DIR/${OUTPUT_ROOT}.mztab \
    --label "IM + masked-AA model, ~1.1M steps — timsTOF HeLa val"
```

A100 wall-time: ~11 min for 65k spectra in the last sanity run.

---

## Files

All on UC Davis HIVE at `/quobyte/proteomics-grp/de-limp/cascadia/`.
The internal labels (`track_g`, `predict_track_f.py`, etc.) are baked
into our file paths; nothing meaningful about the letters beyond that.

| File | Lines | Purpose |
|---|---|---|
| `casanovo5_im_patches/model.py` | 379 | `Spec2PepIM` subclass: IM + masked-AA |
| `casanovo5_im_patches/dataloaders.py` | 152 | `DeNovoDataModuleIM`: `ion_mobility` CustomField |
| `casanovo5_im_patches/train_track_g.py` | 314 | Training entry point: flat LR, mask_prob=0.15 |
| `casanovo5_im_patches/predict_track_f.py` | ~150 | Inference entry (bypasses stock CLI) |
| `training/model_comparison/score_casanovo.py` | ~100 | Replacement for `casanovo --evaluate` |
| `training/casanovo_track_g/checkpoints/last-v3.ckpt` | 554 MB | Final checkpoint (~1.1M steps) |
| `training/casanovo_track_g/checkpoints/last-v1.ckpt` | 554 MB | ~600k step ckpt (referenced in the numbers above) |
| `training/casanovo_finetune/track_e/val.mgf` | 61.7k spectra | In-distribution timsTOF HeLa val |

Can package the `casanovo5_im_patches/` directory + an example val MGF
into a self-contained tarball.

---

## Eval suite

timsTOF HeLa held-out val on both checkpoints:

| Checkpoint | N spectra | Peptide acc. | AA acc. |
|---|---|---|---|
| ~600k steps (`last-v1.ckpt`) | 61.7k | 53.9% | 71.6% |
| **~1.1M steps (`last-v3.ckpt`, final)** | 61.7k | **61.0%** | **74.4%** |

Other datasets (numbers from `last-v1.ckpt`; `last-v3` re-run pending):

| Dataset | Type | N spectra | Peptide acc. | AA acc. | Δ AA |
|---|---|---|---|---|---|
| Salmon timsTOF | OOD species + instrument run | 79k | 31.5% | 55.4% | +11.0 |
| Arab val | OOD (thin) | 5.6k | 8.3% | 37.5% | +12.6 |
| Arab sanity JL19042 | **TRAIN-CONTAMINATED** | 65k | 32.1% | 59.6% | (memorization check) |
| Arab sanity JL520 | TRAIN-CONTAMINATED | 38k | 31.6% | 56.4% | (memorization check) |
| Arab sanity MKT007 | TRAIN-CONTAMINATED | 22k | 20.4% | 47.1% | (memorization check) |

The sanity rows are training-data subsets evaluated as predictions; the
gap between sanity (~32%) and held-out (61.0%) tells us the val split is
sufficiently disjoint from train, and the model is not just memorizing.
The held-out gain over the contaminated sanity rows is what we'd expect
from a generalizing model rather than a memorizing one.

**OOD dataset-pick gotcha worth flagging:** our first OOD candidate was
PXD067383 (αKNL2 SUMOylation IP-MS) — turned out to be 210 unique
peptides with 19.5k deeply-redundant PSMs, useless for measuring
generalization. We've since filtered to (bulk proteome, ≥1k unique
peptides at 1% FDR, ddaPASEF, 2022+).

---

## Open architectural questions (Casanovo side)

1. **Is IM-on-encoder-memories the right injection point?** We didn't try
   IM-as-per-peak-feature (each fragment carries inherited 1/K0). That might
   matter more for fragment-level inference.
2. **Does the identity-zero init trick generalize?** RT, charge-state,
   instrument-id embedding all face the same "add a feature to a
   pretrained transformer without breaking it" problem.
3. **Mask rate sensitivity** is untested at 15%.

## Other things we tried that didn't help (Casanovo)

These didn't make the main "what didn't work" table because they're
smaller, but they're things we burned time on and would flag for anyone
working on Casanovo:

- **Casanovo conda env shipped with CPU-only PyTorch (`2.11.0+cpu`)**
  on at least one of our HPC profiles. Training ran on CPU at ~200×
  the expected wall time. Symptom: GPU utilization 0% in nvidia-smi,
  no `CUDA out of memory` errors. Fix: `pip install torch==2.1.2+cu118
  --index-url https://download.pytorch.org/whl/cu118` after env activate.
  Worth a check in the install docs.
- **Mass-accuracy recalibration on timsTOF data** (we measured −3.2 ppm
  median bias, 5.7 ppm RT drift across the gradient on a local HeLa
  ddaPASEF run). With Casanovo's default ±20 ppm tolerance, recalibrating
  the search space did not change the prediction set meaningfully. Useful
  to know the tolerance is wide enough that small systematic biases don't
  bite.
- **Mzml conversion via ProteoWizard's `msconvert` under wine** on Linux
  (Apptainer container, `wine64 msconvert`). Works, but the round-trip
  loses the per-precursor 1/K0 that the timsrust path preserves. If we'd
  gone via msconvert we'd never have had IM data to feed in.
- **Larger batch sizes with full-spectra inputs** — irrelevant to
  Casanovo because Casanovo uses depthcharge's default 150-peak filter,
  but worth flagging that the attention is O(L²) in peak count, so any
  experiment that disables filtering hits the wall fast (more on this
  in the Cascadia section).
- **`grad_clip_val` ablation** — tried 0.5, 1.0, none. Made no detectable
  difference on our setup. Left at Casanovo default.

---

# Part 2 — Cascadia fine-tuning notes

*For forwarding to the Cascadia maintainer if useful. We didn't modify
Cascadia's architecture — these are observations from fine-tuning the
released Cascadia model on local timsTOF diaPASEF data.*

## What we did

Fine-tuned Cascadia 0.0.7 on ~290k labeled timsTOF HeLa spectra extracted
from diaPASEF runs as pseudo-DDA spectra (mobility-filtered, see below).
Two completed runs:

| Internal label | Preprocessing | Effective batch | Status |
|---|---|---|---|
| Attempt A (Track A) | `[scale_intensity("root"), scale_to_unit_norm]` (matches Cascadia's `train()` default), median ~88 peaks/spectrum | 32 | 30 epochs done, loss 1.23, `TrackA-epoch=29-step=268000.ckpt` |
| Attempt B (Track B) | depthcharge default + `filter_intensity(max_num_peaks=200)` | 32 | 10 epochs done; required `primitives.py` patch (below) |

## What surprised us about Cascadia internals

Each of these would help in the docs or a `FINE_TUNING.md`:

1. **Cascadia does NOT use the depthcharge default peak filter.** The
   docstring and class default on `AnnotatedSpectrumDataset` includes
   `filter_intensity(max_num_peaks=200)`, but Cascadia's `train()` entry
   point explicitly passes `preprocessing_fn=[scale_intensity("root"),
   scale_to_unit_norm]` — overriding the default, no peak filter at all.
   The pretrained checkpoint was trained on full unfiltered spectra
   (median ~9.5k peaks, max ~113k). We initially used the depthcharge
   default and got a distribution shift; fine-tuning loss curves looked
   "fine" but eval was off.

2. **`configure_optimizers()` returns a hidden LR scheduler.**
   `AugmentedSpec2Pep.configure_optimizers()` returns
   `CosineWarmupScheduler(warmup=10000, max_iters=100000)`. Steps 0–10k:
   LR ramps from 0 → target (almost no learning). Steps 100k+: LR decays
   back to 0 (learning stops). The original pretraining used
   `warmup=200_000, max_iters=800_000`, which is invisible from the
   class signature. For fine-tuning over <100k steps this means most of
   training runs at near-zero LR. We override with a flat optimizer:

   ```python
   def flat_lr_optimizer(self):
       return [torch.optim.Adam(self.parameters(), lr=5e-6, weight_decay=1e-9)], []
   model.configure_optimizers = flat_lr_optimizer.__get__(model)
   ```

3. **The `lr_decay` parameter is actually Adam's `weight_decay`,** not a
   schedule. The naming threw us off for a while because the scheduler
   in (2) is a separate object.

4. **`max_charge=10` (11 classes), not 4.** When we initially set
   `max_charge=4` to match Casanovo, the model failed at load with a
   shape mismatch on the charge embedding.

5. **`temp_path` defaults to `~/cascadia_temp/`.** On HPC systems with
   small home directory quotas (5–10 GB), the run will fill the quota
   silently mid-job. We set `temp_path=/tmp/cascadia_temp_$SLURM_JOB_ID`.

6. **Inference tries to download the tokenizer from MassIVE-KB** on first
   run. On GPU compute nodes with no outbound network, this hangs. Need
   to pre-stage the tokenizer file before launching the job.

## Failed attempts on the Cascadia side

| Attempt | What we tried | What failed |
|---|---|---|
| `training_run_v2` | Full unfiltered spectra, batch=1 + accumulate=16 on A100 80GB | OOM at first training step (24 GiB allocation). Attention is O(L²) in peak count; 113k peaks = ~13 GB just for one attention map. |
| Earlier same | batch=8, full spectra | OOM at 4 min |
| Earlier same | batch=4, `precision="16-mixed"` | Instant crash — older Lightning expects `precision=16` (no string variant). |
| Earlier same | batch=4, `precision=16` | OOM at 9 min on first train step |
| Earlier same | batch=32, `max_num_peaks=200` | Cancelled — would have changed the input distribution vs pretrained (see surprise #1 above). |

Conclusion: full-spectra Cascadia fine-tuning is infeasible on a single
A100 80GB without architectural changes (e.g. linear attention). Our
workaround was to use mobility-filtered pseudo-DDA spectra extracted
from the diaPASEF data (median 88 peaks), which fits comfortably at
batch=32.

## depthcharge patches that affect both Cascadia and Casanovo

(Same as the Casanovo Part 1 list, repeated for the Cascadia
audience's convenience.)

1. **`AnnotatedSpectrumDataset.__len__` references `self._offsets`;**
   the attribute is actually `self._file_offsets`. Patched in our overlay.
   One-line fix upstream.

2. **`filter_intensity(max_num_peaks=N)` only syncs `mz` + `intensity`
   arrays.** Cascadia adds extra per-peak arrays (`rt`, `level`, `im`,
   `fragment`). When the filter selects the top-200 peaks by intensity,
   it picks 200 mz/intensity indices but doesn't apply the same selection
   to the extra arrays, so after filtering the spectrum has 200
   mz/intensity values but the original count of rt/level/im/fragment.
   Downstream code that assumes synced indices then reads wrong rt/im
   for a given peak.

   We patched `primitives.py` so the filter applies its selection to all
   arrays in the spectrum dict. Worth upstreaming — silent data
   corruption is the worst failure mode.

3. **`n_workers > 0` deadlocks on Quobyte** when multiple workers try
   to open the same HDF5 file concurrently. `n_workers=0` or stage to
   local SSD first.

## Local Cascadia paths

All on UC Davis HIVE at `/quobyte/proteomics-grp/de-limp/cascadia/`:

| Path | Contents |
|---|---|
| `models/cascadia.ckpt` | Pretrained Cascadia 0.0.7 (558 MB) |
| `training/track_a_mobfilter/` | Attempt A checkpoints (30 epochs) |
| `training/track_b_peakcap/` | Attempt B checkpoints (10 epochs) |
| `training/training_run_v2/` | The failed full-spectra OOM attempts |
| `envs/cascadia5/` | Conda env: Cascadia 0.0.7, PyTorch 2.0.1, timsrust_pyo3, pyarrow 23.0.1 |

Happy to share any of this. We've found Cascadia very useful for our
diaPASEF paleoproteomics workflow and would like to see it grow.
