# qs_framework — quasi-steady FSI analysis for the RC-19 panel

One code base. Data per case. Settings per study. Nothing is ever copied to
run something new.

```
qs_framework/
├── qs_startup.m       run this once per MATLAB session
├── src/               the code. NEVER edited to run a case or a study.
├── cases/<case_id>/   data for one CFD case + its case.json
├── studies/<id>.json  what to do with a case
├── shared/            the ROM and defaults.json
├── slurm/             launchers
├── tests/             the five checks
├── docs/              CODES.md, the design memo, the examples sheet
└── results/<case>/<study>/   everything a run produces
```

---

## Quick start

```matlab
run('/path/to/qs_framework/qs_startup.m')     % once per session, from anywhere

% check the plumbing before anything expensive
run_all_tests([1 2 5])

% what does label 3 mean? what is Amp_wh?
qs_codes

% a sweep
qs_init_study('heated_c1_NoSBLI_Periodic','sweep_coarse_v1')   % once
qs_run_sweep ('heated_c1_NoSBLI_Periodic','sweep_coarse_v1',1,1)  % or via SLURM
qs_merge     ('heated_c1_NoSBLI_Periodic','sweep_coarse_v1')
qs_plot_map  ('heated_c1_NoSBLI_Periodic','sweep_coarse_v1')

% refine where the map is undecided, then resubmit
qs_refine_points('heated_c1_NoSBLI_Periodic','sweep_coarse_v1', true)   % dry run first
qs_refine_points('heated_c1_NoSBLI_Periodic','sweep_coarse_v1')

% a long run, same code path
qs_init_study('heated_c1_NoSBLI_Periodic','long_nominal_20s_v1')
qs_run_long  ('heated_c1_NoSBLI_Periodic','long_nominal_20s_v1')
```

On the cluster:

```bash
bash slurm/submit_sweep.sh heated_c1_NoSBLI_Periodic sweep_coarse_v1 48
```

---

## Adding a heated case

1. `cases/heated_cN_.../aero/` — drop in the two Fluent exports
   (panel surface pressure; BL-edge pressure + Mach + sound speed).
   `.prof` as exported, or `.csv` with a header row.
2. `cases/heated_cN_.../ref/` — drop in the DIC results `.mat`.
3. Edit `case.json`: point `aero_build.surf_file` / `edge_file` at the two
   files. Everything else is already filled in from the challenge
   `caseN_conditions.json`.
4. Build the aero data and check it:
   ```matlab
   qs_build_aero('heated_c1_NoSBLI_Periodic')
   ```
   Read the `implied T0` line it prints. For a heated case it must come out
   near 388 K. If it says ~291 K, that is the old cold-flow RANS and the
   sweep will under-predict limit cycles by roughly 15% in damping.
5. `qs_init_study(...)` and go.

If the panel changed when it was installed (`f1_installed_Hz` is a
measurement, 263.67 ± 9.77 Hz), give that case its own ROM in
`cases/<id>/rom/` and point `rom_file` at it instead of `shared/`.

---

## The three ideas worth knowing

### 1. Point list, not a grid

The old sweep thought in terms of an `npc × ndT` array with linear indices.
Refinement and per-point run lengths break that. Here a study is a **list of
runs**, each with an ID derived from `(p_c, dT, t_end, ic_tag)`:

```
pc051927_dT0012p760_t003p00_flat
```

Consequences, all of which fall out for free:

- refinement **appends**; nothing is renumbered, so old results stay valid
- a slice skips IDs it already has, so **resume is just resubmitting**
- **per-point `t_end`** is one column, not a special mechanism
- merge is a **join on ID**; the rectangular grid is rebuilt only for plotting
- a 10 s re-run of a 3 s point is a **new row, not an overwrite** — both
  answers stay on the record, and `qs_rasterize` shows the longest one

### 2. Coarse → classify → refine

`sweep_coarse_v1` is deliberately coarse (0.25 kPa × 0.5 K, ~675 points,
~80 core-hours). Its job is to *locate* the boundary. `qs_refine_points` then

- bisects adjacent pairs of cells whose labels differ, down to
  `refine.min_dpc_kPa` / `min_ddT_K` — but **only where the boundary is
  locally coherent**, i.e. both cells agree with a majority of their own
  neighbours. Near the transition this map speckles (two attractors coexist
  and which one a flat start reaches varies erratically), and bisecting a
  speckle never converges — it just spends compute,
- re-queues at 10 s (then 20 s) any point flagged `nonstationary`,
  `near_threshold` or `transition`, and
- for speckled cells, re-runs them longer **and** from the mirrored initial
  condition (`refine.ic_probe`), which is the measurement that actually tests
  whether it is a basin effect.

Repeat merge → refine → resubmit until it reports nothing to add.

### 3. Mid-run transitions

Every run is cut into overlapping windows (0.5 s, 50% overlap) and each
window gets its own label. A point whose label changes and stays changed is a
**transition**, with a time attached. That matters because a calibrated
case-1 point was seen to snap to static at ~7 s — every 3 s label near a
boundary is provisional until the windows say it settled.

Storage follows from that: the **window signature is always kept** (a few
hundred bytes per point), the **full centre history only for flagged points**.
Keeping every history for a 2000-point sweep would be ~1.6 GB.

The primary label is taken from the **final window**, not from statistics over
the whole post-transient record: with a transition in the record, those
statistics blend two regimes and describe neither.

Transition times are biased late by up to about half a window plus the time
the amplitude took to cross the threshold. `trans_t_lo` / `trans_t_hi` bracket
it. Don't quote one to better than the window length.

---

## Configuration

Three layers, later wins:

```
shared/defaults.json  <-  cases/<id>/case.json  <-  studies/<id>.json  <-  call-site overrides
```

`case.json` is the physics of a CFD case. `studies/<id>.json` is what to do
with it. **Those are the only two files you edit.** The resolved settings and
the git commit are written into `results/.../manifest.json`, so every result
carries the configuration that produced it.

Keys beginning with `_` are comments and are ignored.

---

## Preflight

Runs automatically before a long run; set `preflight.run_on_sweep = true` to
run it on sweeps too. Five checks — the first three each cost a full sweep on
this project before they existed, and the fourth was added after a study was
pointed at a case whose operating point sat 26σ outside its own swept box:

| check | what it catches |
|---|---|
| implied `T0` from `a_edge`, `M_edge` | a cold-flow RANS behind a heated case: `a_edge` low by √(388/291) makes the piston-theory damping term `Zdot/a_l` ~15% too large, which suppresses limit cycles and shrinks the flutter region. The map still looks plausible. |
| `mean(p_l) − p_c` vs the DIC mean deflection | the panel being pushed the wrong way, which makes every mean-field RMSE meaningless. Set `case.net_load_sign_convention` to `+1`/`-1` to enable it. |
| `dT_cr` from `Ke v = dT·Kt v` | sweeping mostly below the buckling threshold |
| the study box contains the case nominal | a study written for one case pointed at another whose p_c sits somewhere else entirely. The nominal must be **inside** the box, with `preflight.pc_margin_sigma` (3) of p_c margin. `dT_margin_sigma` defaults to 0 = off, because σ(ΔT) ≈ 3.4 K is large while the LCO onset is narrow and mid-box — a σ-margin there has no physical content |
| minimum BL-edge Mach clear of 1 | `c1 = M/√(M²−1)` is singular at M = 1; subsonic edge cells make the coefficients complex and nothing downstream would tell you |

`preflight.strict = true` turns a finding into an error.

---

## Tests

```matlab
run_all_tests([1 2 5])   % fast: plumbing only, no integration
run_all_tests            % adds the regression and round-trip tests
```

`qs_startup` puts `tests/` on the path, so these work from any directory.

| test | what it proves | cost |
|---|---|---|
| t01 | config resolution, grid construction, IDs, slicing, append | < 1 s |
| t02 | model build matches the original driver; preflight fires | ~10 s |
| t03 | **metrics reproduce the original case-1 map to 1e-6**, with absolute floors so near-zero quantities are not compared as ratios | ~3 min/point |
| t03b | `t03b_attribute(pc,dT)` — why one point differs: settings, or trajectory sensitivity | ~6 min |
| t04 | init → slices → resume → merge → rasterize → refine | ~1 min |
| t05 | windows / classification / transition on known signals | ~5 s |

A report is written to `tests/last_test_report.txt` with the failing values,
not just the word FAIL. Send that file back if something breaks.

**t03 is the one that matters.** The six core files (`rc19_setup_rom`,
`rc19_solve_structural`, `rc19_nonlinear_force`, `rc19_reconstruct`,
`cm_setup_aero`, `pm_quasisteady_ept`) moved into `src/` byte-for-byte, so if
the configuration plumbing did not change a number anywhere, every metric must
come back identical. Anything else means a value that used to be hardcoded is
now being read from the wrong place.

t02 asserts that preflight **fails** on `unheated_c1`: that case's aero file
is the known cold-flow RANS, so a pass there would mean the gate is broken.

---

## Cluster notes

`submit_sweep.sh` uses N independent `sbatch` calls with `SLICE`/`NSLICE`
exported, **not** a job array, and passes `--exclude=ec77` on the command line.
Node `ec77` kills every job it receives at launch (1 s elapsed, `CANCELLED by
0`, before `module load` runs), and the no-array fan-out is the version that
has been verified to dispatch every task. Node names are not baked into the
`.slurm` files because they go stale.

After any fan-out, check where jobs landed, not just their state:

```bash
sacct -j <ids> -X --format=JobID%12,State%20,Elapsed,NodeList%20
sinfo -R
```

Resubmitting a lost slice is safe: finished points are skipped.

---

## Syncing

`results/` is excluded from `.gitignore` and from the up-sync. Push code and
cases up, pull results down:

```bash
rsync -av --exclude results/ qs_framework/ cluster:/path/qs_framework/
rsync -av cluster:/path/qs_framework/results/ qs_framework/results/
```

---

## Codes and conventions

Results tables are full of numbers that stand for things — `label = 3`,
`level = 1`, `signflip = -1`. What each one means is written down in one place,
**`docs/CODES.md`**, and printed by:

```matlab
qs_codes                % everything
qs_codes('labels')      % labels flags metrics points config files units
```

`qs_codes` builds the label section by calling `qs_label_name`, so it cannot
drift from what the classifier assigns. **Adding a coded value means editing
`qs_codes.m` and `docs/CODES.md` in the same commit.**

The short version:

| label | | flags (independent of the label) |
|---|---|---|
| 0 `unknown` | | `nonstationary`, `drift_tested`, `near_threshold`, |
| 1 `static` | | `divergent`, `short_record`, `trans_flag`, `ok` |
| 2 `LCO` | | |
| 3 `broadband` | | |
| 4 `divergent` | | |

There is no label for "nonstationary" — a run can be `static` **and**
`nonstationary`, which is the whole point of keeping them apart: one says what
the panel did, the other says whether to believe it.

`Amp_wh` is `(max − min)/(2h)` of the **centre node** over the whole
post-transient record — a half-range, so for a sinusoid it is the single-sided
amplitude in panel thicknesses. `amp_last` is the mean window amplitude over the
last half, and is the one the label actually describes; the map plots that.
`MeanPeak` / `StdPeak` are time statistics within one run, peaked over the
panel. **Nothing is ever averaged across sweep points.**

---

## Known state

- `cases/unheated_c1_NoSBLI_Periodic` is the **regression baseline** and uses
  the cold-flow aero file. Do not "fix" its numbers; make a new case.
- The heated case folders are **stubs** waiting for the profile exports and
  the DIC files.
- `Minf` is 1.92 throughout. AFRL's nominal is 2.0, and the choice moves
  `p_inf` by ~5.9 kPa. That is the open question behind the static-load sign
  issue and is not resolved here.
