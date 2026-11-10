# Additions from Magda — November 2026

Two sets of things: the **paper figures** (`1.1`, `2.2`, `4.1`, with a legend file
each) and the **earlier exploratory scripts** (`4`–`7`, with the PDFs under
`paper_graphs/magda/`). All paths are relative to the repository root, so run them
from the project root (or with the `.Rproj` open).

Read the paper scripts first — `4`–`7` are kept because they cover traits and
maternal lines that the paper figures do not, but where the two overlap the paper
versions supersede them (see "How the paper scripts differ" below).

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
- It also renames `Roi.Mask` to `PSA`, so that is the column name downstream.

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

## The paper figures

Everything in these three scripts is fitted on **maternal-line means**, never on
individual plants — 15 points per group per timepoint (5 regimes × 3 lines), 10 for
Oregon 2014, which has only 2 lines. Each writes its statistics to `data/` as CSV
and its figure to `paper_graphs/`, and each has a companion `*_legend.txt` in
`paper_graphs/` with the draft legend, the full n accounting, what every graphical
element means, and the caveats worth keeping in the text.

| Script | Figure | Legend file |
|---|---|---|
| `1.1_Fig1_PSA_GR_MMJ.R` | `Figure 1.pdf` — (A) PSA over time per watering regime with a −log10(P) row under each row of curves, (B) growth rate for the early and late intervals with ANOVA P and Tukey letters, (C) GR regressed on % SWHC per group and interval. | `Fig1_PSA_GR_SWHC_legend.txt` |
| `2.2_Fig2_FvFm_ChlIdx_MMJ.R` | `Figure 2.pdf` — the same construction as Fig 1A for (A) Fv/Fm and (B) Chlorophyll Index. | `Fig2_FvFm_ChlIdx_legend.txt` |
| `4.1_Fig4_Regression_TreeRing_MMJ.R` | `Figure 4.pdf` — (A) NDVI regressed on % SWHC at 144 h for all four groups, (B) the tree-ring: that regression repeated at every timepoint and for every trait. | `Fig4_Regression_TreeRing_legend.txt` |

Figure 3 is the Fv/Fm vs ChlIdx lag analysis and is yours — nothing here touches it.

Statistics in Fig 1A and Fig 2: one-way ANOVA across the five regimes at each
timepoint within each group, then Tukey HSD, plotting only the four contrasts
against **100% SWHC**. The plotted value is −log10 of Tukey's p-adj, so the four
comparisons at that timepoint are already accounted for; the omnibus ANOVA p is
kept in the output table as `p_anova`.

Statistics in Fig 4: `lm(trait ~ SWHC)` at every timepoint × trait × group, with BH
correction across the 52 timepoints within each trait × group.

The run break at 409.5 h is drawn in every panel — a dashed vertical line in
Figures 1 and 2, a dashed ring in the tree-ring. It is where the plants outgrew the
imaging trays and the setup changed, so any step across it is technical.

### How the paper scripts differ from `4`–`7`

Three things changed, and all three matter if you compare the two sets of figures:

1. **Sign convention.** `4` and `5` fit `trait ~ StressDose`, where
   `StressDose = 100 − % SWHC`, so larger = drier. `4.1` fits `trait ~ SWHC`
   directly. **The slopes are therefore opposite in sign between the two.** In
   `4.1`, a negative slope means the trait falls as the pots get wetter, i.e. it is
   higher under drought.
2. **Fill.** `4` and `5` fill the ring with unsigned −log10(p), which hides
   direction. `4.1` fills with `sign(slope) × −log10(P)` on a diverging scale
   centred on zero: green = the trait rises with wetter soil, purple = it falls,
   near-white = no response either way. The raw slope cannot go on one colour scale
   across the 14 traits because the units differ (PSA is an area in a.u., Fv/Fm a
   ratio between 0 and 1), so the sign is taken from the slope and the magnitude
   from a quantity that is comparable. `FILL_MODE` at the top of `4.1` switches to
   `"signed_r"` for `sign(slope) × sqrt(R2)` (effect size instead of evidence) or
   `"logp"` for the old unsigned version.
3. **Tile height.** `geom_tile` takes its height from the *smallest* gap between
   TOE values. The grid is 12 h everywhere except across the run break
   (408 → 411 h = 3 h), so left to itself it draws every tile a quarter of its
   proper size and the ring comes out as slivers separated by white. `4.1` sets the
   height to the modal spacing. `4` and `5` still have this.

## The earlier scripts

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

These four run on **all plants**, not line means, which is why their strips are
much busier than Figure 1A. They are the right place to look for a trait that is
not in Figures 1, 2 or 4, or for what an individual maternal line is doing.

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

`Figure 1.pdf`, `Figure 2.pdf` and `Figure 4.pdf` are the ggplot output of the
scripts above and will reproduce exactly.

The PDFs committed under `paper_graphs/magda/` were produced from the same analysis
as scripts `4`–`7` but with a matplotlib implementation, so re-running those scripts
will reproduce the analysis and the layout in ggplot without being pixel-identical.

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

At line means, the four groups differ enormously in how much of the trait panel
responds at all: in Figure 4B, 902 of 2912 trait × group × timepoint regressions
pass FDR < 0.05, but California 2014 accounts for 465 of them (13 of 14 traits)
against 88 for Oregon 2010 (4 traits). Part of the Oregon 2014 shortfall is power —
2 lines, so 10 points per regression rather than 15.
