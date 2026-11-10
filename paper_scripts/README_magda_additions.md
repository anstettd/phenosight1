# Additions from Magda — November 2026

Scripts `4`–`7` and the figures under `paper_graphs/magda/`. All paths are relative
to the repository root, so run them from the project root (or with the `.Rproj` open).

## Change to `0_data_setup.R`

Dan — the growth-rate block is now filled in, as you asked in the comment there.

- The TOE cut-off is a single variable, `TOE_MIN <- 140`, at the top of the script,
  and it is now used for **both** the PS2 filter and the growth-rate fit.
- GR is fitted **per plant** and only then averaged to line means, rather than
  fitting one slope to the line-mean series. It keeps the within-line variation and
  lets the `R2 > 0.7` filter drop individual plants whose growth was not linear;
  fitting line means first would hide both. So it runs off `PS2_early2` /
  `PS2_late2` (plant level), not `PS2_early3` / `PS2_late3`.
- It writes `data/GR_early3.csv`, `data/GR_late3.csv` and a combined `data/GR_all.csv`.

Note this changes the numbers slightly from `Experiment01_*Growth_Mimulus.csv`:
those used `TOE > 120`, and only for Early — the Late frame starts at 411 h, so the
cut-off never applied there. Both windows now use `TOE_MIN`, which drops the 132 h
point from Early.

## Why 144 h

Projected shoot area bottoms out at exactly 144 h in all four site × year groups.
Before that the splines are erratic on very small plants, and in an earlier pass
most of the apparent 100%-contrast signal came from those first few points. With
`TOE_MIN <- 140` the first grid point kept is 144 h, so every series starts at its
own minimum.

## The scripts

| Script | What it makes |
|---|---|
| `4_treering_dose_response.R` | Tree-ring figure: radius = TOE, angle = trait, one 90° quadrant per site × year group. Fill is −log10(p) of the dose slope, red stars at FDR < 0.05. |
| `5_treering_dose_response_byLine.R` | Same, but each quadrant is split into its maternal lines. Oregon 2014 has two lines, so its third wedge is blank. |
| `6_timeseries_byGroup.R` | Per trait: 2 × 2 panels (rows = site, columns = year), all five watering regimes as mean ± SE, with a −log10(p) strip under each panel. |
| `7_timeseries_byLine.R` | Same, one panel per maternal line (rows = group, columns = lines). |

Statistics, in `6` and `7`: one-way ANOVA across the five regimes at each timepoint,
then Tukey HSD, keeping the contrasts against **100% SWHC**. The
plotted p is Tukey's p-adj, so the four comparisons at that timepoint are already
accounted for; the omnibus ANOVA p is kept in the output table as `p_anova`. Each
line in the strip is labelled by the regime it represents.

In `4` and `5`: `trait ~ StressDose` fitted at every timepoint, where
`StressDose = 100 − % SWHC` (so larger = drier), with BH correction
across timepoints within each trait × group.

## Two things to know about the data

**The two runs are different plants.** `Experiment01_EarlyPS2_Mimulus.csv` (Exp 20243,
0–408 h) and `Experiment01_LatePS2_Mimulus.csv` (Exp 20244, 411–747 h) share no
`Plant.ID` — the tray IDs differ (`Tray01` vs `Tray01A`). Every figure here breaks
its lines at 409.5 h rather than joining across it, and marks the break. The jump
there is a cohort change, not growth.

**The p-values run optimistic.** These tables are the spline outputs on the 12 h grid
and carry no `Data.Support` column, so the fits include smoothed values. Smoothing
shrinks the residual. Fine for ranking traits; not for quoting.

## Provenance of the PDFs

The PDFs committed under `paper_graphs/magda/` were produced from the same analysis
but with a matplotlib implementation, so re-running scripts `4`–`7` will reproduce
the analysis and the layout in ggplot without being pixel-identical.

## Results, in one paragraph

Against 100%, only **20%** separates consistently, in every group, and for shoot
area and Fv/Fm it is the only regime that does. **AriIdx and NDVI** behave
differently: 40% and 60% separate about as strongly as 20%, so they track a graded
response across the whole gradient rather than a drought threshold. **80%** is
almost never different from 100% except in NDVI (Oregon 2014, California 2010).
**qL responds nowhere** and qP almost nowhere. Fv/Fm is *higher* under 20% watering
in all four groups.

At maternal-line level there are only 3–5 plants per regime per line, so those
tests are much weaker — read the curves, not the strip. Two things show up only at
that resolution: line **H** is flat and much smaller than G and I, and in California
2010 line **A** is flat across the whole gradient while B and C are not, which is
what made that group look quiet when the lines were pooled.
