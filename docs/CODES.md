# Codes and conventions

Every coded value the framework writes into a results table, and what it means.

From MATLAB the same legend is available without opening this file:

```matlab
qs_codes                % everything
qs_codes('labels')      % one section: labels flags metrics points config files units
```

`qs_codes` derives the label names by calling `qs_label_name`, so that section
cannot drift from what the classifier actually assigns. **If you add a coded
value, add it to `qs_codes.m` and to this file in the same commit.**

---

## Regime labels

Fields carrying a label code: `label`, `label_start`, `trans_from`, `trans_to`,
and every entry of the per-window `W.label` array.

| code | name | test |
|---|---|---|
| 0 | `unknown` | spectrum could not be formed (record too short) |
| 1 | `static` | window amplitude < `classify.amp_tol_wh` (default 0.10 w/h) |
| 2 | `LCO` | above tolerance **and** `pfrac ≥ classify.peak_power_frac` (0.35) — one dominant tone |
| 3 | `broadband` | above tolerance, no dominant tone |
| 4 | `divergent` | whole record: log-linear fit of window amplitude grows by more than `classify.growth_tol` (3×) and ends above tolerance |

**These five are the complete set.** There is no label for "nonstationary" —
that is a flag. A run can be `static` and `nonstationary` at the same time, and
keeping them separate is the point: one says what the panel was doing, the other
says whether to believe it.

A label is assigned **per analysis window**. A run's label is the label of its
**final window**, because when a record contains a mid-run transition the
whole-record statistics blend two regimes and describe neither.

Plot colours, in `qs_plot_map`: 1 pale grey, 2 blue, 3 orange, 4 red.

---

## Flags

Logical, one per run, all independent of the label.

| flag | meaning |
|---|---|
| `nonstationary` | amplitude drifted by more than `classify.stationarity_tol` (0.15) between the first and last half of the record, **or** the window labels were not all equal |
| `drift_tested` | whether that drift test was applied at all. See the note below — a clean `nonstationary = 0` with `drift_tested = 0` means "never asked", not "settled" |
| `near_threshold` | `\|amp_last − amp_tol_wh\|` lies within `classify.near_threshold_band` (±50%) of the tolerance: the static/LCO call could flip given more time |
| `divergent` | the same test that assigns label 4 |
| `short_record` | fewer windows than the transition test needs |
| `trans_flag` | a regime change was found that persisted to the end of the run |
| `ok` | the integration completed. False means read `errmsg` |

**Why `drift_tested` exists.** The drift test is a *ratio*. A panel sitting
still has a residual of about 1e-8 panel thicknesses, and the ratio between the
first and last half of that is roundoff — it exceeds 15% every single time.
Without a floor the flag fires on every static cell in the map, and since
`qs_refine_points` promotes flagged points to 10 s and then 20 s runs, that
spends hundreds of core-hours re-running dead panels. So the test is skipped
below `classify.stationarity_amp_floor` (0.01 w/h — a tenth of the static
tolerance, low enough that a borderline 0.07 w/h point is still judged, high
enough that numerical residue is not), and `drift_tested` records whether it ran.

---

## Metrics

One per run. **Nothing is averaged across the sweep** — every point is computed
in complete isolation. Within a run, `mean` and `std` are taken over **time**,
and the "peak" is then taken over **space**, across the panel. The amplitude
quantities are the **centre node only**.

| metric | definition | units |
|---|---|---|
| `Amp_wh` | `(max − min) / (2h)` of the centre node over the **whole** post-transient record. Half-range, so for a sinusoid it is the single-sided amplitude in panel thicknesses | w/h |
| `amp_last` | mean **window** amplitude over the last half of the record. This is the amplitude the label describes — prefer it on a map | w/h |
| `amp_first` | the same over the first half | w/h |
| `amp_ref` | `max(amp_first, amp_last)`; what gates the drift test | w/h |
| `growth_ratio` | `amp_last / amp_first` | – |
| `Freq` | dominant frequency of the centre-node spectrum | Hz |
| `pfrac` | power within ±5% of the peak frequency, over the total | – |
| `MeanPeak` | max \|time-mean deflection\| anywhere on the panel | mm |
| `StdPeak` | max time-std of deflection anywhere on the panel | mm |
| `MeanRMSE` | ‖model mean field − DIC mean field‖ / ‖DIC mean field‖, at the DIC points | – |
| `StdRMSE` | the same for the std field | – |
| `signflip` | +1 or −1: which orientation was used to match the DIC sign | – |
| `trans_time` | transition time estimate | s |
| `trans_t_lo`, `trans_t_hi` | the bracketing window's start and end | s |
| `wall_s` | wall-clock time for this point | s |

**`Amp_wh` and `amp_last` are not the same thing.** `Amp_wh` covers the whole
post-transient record; the label comes from the final window. On a point that
oscillated and then died, `Amp_wh` stays large while the label reads `static`.
That is why `qs_plot_map` draws `amp_last` next to the label panel — so the two
panels cannot contradict each other.

**Transition times are biased late** by up to half a window plus however long
the amplitude took to cross the threshold. Do not quote one to finer than the
window length; `trans_t_lo`/`trans_t_hi` bracket it.

**The classification uses none of `MeanPeak`, `StdPeak`, `MeanRMSE`,
`StdRMSE`.** Those are reported for comparison against the measurement. Labels
come only from the centre-node window analysis.

### Window table

Kept for **every** point, not just interesting ones: `tmid`, `amp_wh`, `wmean`,
`fdom`, `pfrac`, `label`. A few hundred bytes per run, which is what makes
transitions visible without storing full histories.

---

## Point list columns

| column | meaning |
|---|---|
| `id` | `pc<Pa>_dT<K>_t<s>_<ic>`, e.g. `pc051927_dT0012p760_t003p00_flat` |
| `pc_Pa`, `dT_K` | the operating point |
| `t_end`, `t_trans` | run length and discarded transient, seconds — **per point** |
| `level` | 0 = original coarse grid; 1, 2, … = refinement pass number |
| `origin` | `coarse` \| `refine` \| `manual` \| `long` |
| `ic_tag` | `flat` \| `flat_neg` \| `rest` \| `from:<point id>` |
| `save_hist` | force keeping the full centre history for this point |

Initial conditions, in `qs_initial_state`:

| tag | meaning |
|---|---|
| `flat` | `+seed` on every modal coordinate. What every original sweep used |
| `flat_neg` | `−seed`. The mirrored start; the cheapest possible basin probe |
| `rest` | exactly zero. Linear checks only — never for a buckled panel |
| `from:<id>` | continue from another point's saved final state |

`seed` is `ic.seed`, default 1e-4.

---

## Configuration values that are names, not numbers

| setting | allowed values |
|---|---|
| `study.mode` | `sweep` \| `long` |
| `solver.integrator` | `ode15s` \| `ode45` \| `ode23s` \| `ode113` |
| `ic.mode` | `flat` \| `flat_neg` \| `rest` |
| `aero.c3_mode` | `vandyke2` (c3 = 0, the published Eq. 8) \| `classical_cubic` ((γ+1)/12) |
| `store.save_history_if` | any of `transition` \| `nonstationary` \| `near_threshold` \| `divergent` \| `all` \| `none` |
| `store.history_precision` | `single` (default) \| `double` |
| `duration_policy.promote_if` | any of `boundary` \| `nonstationary` \| `near_threshold` \| `transition` \| `divergent` |
| `fields.mode` | `per_second` (default) \| `per_step` \| `single` |
| `fields.precision` | `single` \| `double` |
| `case.net_load_sign_convention` | `+1` or `−1`. Absent means the load-sign gate stays dormant |

`store.history_precision = single` puts a floor of about 1e-7 w/h on any
run-to-run comparison. Diagnostics that compare two runs need `double`;
`t03b_attribute` sets it automatically.

### Preflight gates

`qs_preflight` returns `F.checks` with one entry per gate, each with a `.pass`:

| gate | what it catches |
|---|---|
| `T0` | implied stagnation temperature from `a_edge`, `M_edge` — a cold-flow RANS behind a heated case |
| `load_sign` | `mean(p_l) − p_c` against the measured mean-deflection direction. Dormant unless `net_load_sign_convention` is set |
| `buckling` | `dT_cr` from `Ke v = dT·Kt v`, against the study's dT range |
| `box` | the case nominal must lie **inside** the swept box, with `preflight.pc_margin_sigma` (3) of p_c margin. `dT_margin_sigma` defaults to 0 = off |
| `mach` | minimum BL-edge Mach clear of 1, where `c1 = M/√(M²−1)` is singular |

`preflight.strict = true` turns a finding into an error instead of a warning.

---

## Files a run produces

| pattern | what it is |
|---|---|
| `pointlist.mat` / `.csv` | the study's list of runs |
| `manifest.json` | resolved settings plus the git commit that produced them |
| `points/partial_<k>_of_<n>.mat` | one worker's results, checkpointed |
| `merged/points.mat` / `.csv` | all workers joined by run id |
| `histories/<id>_hist.mat` | full centre history, for flagged points |
| `histories/<id>_state.mat` | final `[q; qdot]`, for continuing a run |
| `figures/fig_regime_map.*` | the three-panel map, plus its data |
| `<id>/postproc_data.mat` | long runs: the full modal history |
| `<id>/fields/fields_sec_NNN.mat` | long runs: one file per simulated second |

---

## Units — where kPa becomes Pa

| quantity | in JSON | inside the code |
|---|---|---|
| cavity pressure | kPa (`grid.pc_kPa`) | Pa (`pc_Pa`, `pc_nom_Pa`) |
| temperature rise | K | K |
| amplitude | – | w/h, panel thicknesses |
| deflection peaks | – | mm |
| pressures in results | – | Pa |
| time | s | s |

Rule of thumb: anything whose name ends `_Pa`, `_K` or `_m` is SI. A JSON grid
block is in kPa because that is what you read off a gauge.
