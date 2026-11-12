# Additions from Magda — November 2026

Dan's figures are the **main** figures for 1, 2 and 3. The main figures here are
the tree-ring **Figure 4** (`4.1`) and the two contrast half-rings **Figures 5
and 6** (`4.2`), all in `paper_graphs/` with the other main figures. The rest is
**supplemental**: `1.1` and `2.2` under `paper_graphs/Supplemental/magda/`, plus
the **earlier exploratory scripts** (`4`–`7`, with the PDFs under
`paper_graphs/magda/`). All paths are relative to the repository root, so run
them from the project root (or with the `.Rproj` open).

Read `1.1`, `2.2`, `4.1` and `4.2` first — `4`–`7` are kept because they cover
traits and maternal lines those four do not, but where the two overlap the newer
versions supersede them (see "How the newer scripts differ" below).

**Every supplemental figure is written twice**, at two replication levels — see
"Two replication levels" below. Run order is `0_data_setup.R` first, then any of
`1.1`, `2.2`, `4.1` in any order.

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
- Two further tables are written for the plant-level figures, and nothing else
  uses them: `data/PS2_all_plants.csv` (the rows that feed `PS2_early3` /
  `PS2_late3`, not averaged) and `data/GR_all_plants.csv` (the per-plant slopes
  before they are averaged into line means). Both additions are additive — your
  existing outputs are unchanged.

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

## The figures

Each script writes its statistics to `data/` as CSV. Each `*_legend.txt` holds the
draft legend, the full n accounting, what every graphical element means, and the
caveats worth keeping in the text.

| Script | Figures | Where | Legend file |
|---|---|---|---|
| `4.1_Fig4_Regression_TreeRing_MMJ.R` | `Figure4_TreeRing_{lineMeans,plants}.pdf` — the tree-ring: `trait ~ SWHC` at every timepoint for all 14 traits, with the single-timepoint NDVI regression inset in each corner. | `paper_graphs/` (**main**) | `Figure4_TreeRing_legend.txt` |
| `1.1_Fig1_PSA_GR_MMJ.R` | `FigS1_PSA_GR_{lineMeans,plants}.pdf` — (A) PSA over time per watering regime with a −log10(P) row under each row of curves, (B) growth rate for the early and late intervals with ANOVA P and Tukey letters, (C) GR regressed on % SWHC per group and interval. | `Supplemental/magda/` | `FigS1_PSA_GR_legend.txt` |
| `2.2_Fig2_FvFm_ChlIdx_MMJ.R` | `FigS2_FvFm_ChlIdx_{lineMeans,plants}.pdf` — the same construction as S1A for (A) Fv/Fm and (B) Chlorophyll Index. | `Supplemental/magda/` | `FigS2_FvFm_ChlIdx_legend.txt` |
| `4.2_Fig5_6_Contrast_TreeRing_MMJ.R` | `Figure5_*` and `Figure6_*` — eight figures in one run; see "The contrast figures" below. | `paper_graphs/` (**main**) | four files, one per figure × quantity |

## The contrast figures (5 and 6)

`4.2` is a derivative of Figure 4: same geometry and same reading, but each cell
carries a **contrast between two collections** rather than a trait. It writes
eight figures in one run — two contrasts × two quantities × two replication
levels — and needs no editing to do so.

**The contrast** decides what is compared and what the quadrants are. Both are
HALF rings, one quadrant per level of the other factor:

- **Figure 5**, by year: 2014 against 2010, within each site. Quadrants are
  sites, California upper right and Oregon upper left, as in Figure 4.
- **Figure 6**, by site: California against Oregon, within each collection year.
  Quadrants are years, 2014 upper right and 2010 upper left, so the older
  collection reads on the left.

**The quantity** decides what fills a cell:

- **`DeltaPlasticity`** — `(slope_comparison − slope_reference) / residual SD`.
  `trait ~ SWHC` is fitted separately within each level — these are Figure 4's
  regressions — and the slopes are subtracted, then divided by the residual SD of
  the model holding both, which makes the number an effect size and so comparable
  across traits with different units. The asterisks are the `C × SWHC`
  interaction, so fill and asterisks test the same hypothesis.
- **`EvolIndex`** (Figure 5) / **`LevelRatio`** (Figure 6) — an index
  `mean(comparison) / mean(reference)` is formed per watering regime, regressed on
  % SWHC, and the cell carries that slope. The asterisks are the main effect of
  the contrast, which is *not* a test of what the colour shows; the interaction P
  is in the CSV as `P_con_swhc` if you want it. Figure 6's is not called an
  evolutionary index, because a between-site difference is divergence rather than
  evolutionary change.

Why both: the index is mostly a **level** shift and Δ plasticity a **plasticity**
shift, and a population can evolve either without the other. Read radially too —
radius is experimental time, so a wedge that is near-white at the centre and
deepens outward is a difference that *builds* as the experiment runs, rather than
one present from the first image.

**What they show, briefly.** At maternal-line level **no cell** passes FDR < 0.05
for the interaction, in either figure (0 of 1456) — the Δ plasticity panels carry
no asterisks at all — while the level contrasts pass in 836 of 1456 cells for
year and 1080 for site. Taken together: the collections differ in trait level, and
there is no resolvable evidence that they differ in how they respond to water.
Two cautions, both spelled out in the legends. The interaction test is far weaker
than the main-effect test on the same rows, so a blank panel is not evidence of no
change in plasticity. And the main-effect counts are optimistic, since those terms
are fitted across 25–30 rows drawn from only five maternal lines — 57% and 74% of
the panel "significant" is a reason for caution, not celebration. The one trait
that breaks through is NDVI at plant level, consistently negative but absent at
line level, so a candidate rather than a result.

Every cell's numbers, both quantities and both P values, are in
`data/Figure5_contrast_regressions_<unit>.csv` and the Figure 6 equivalent.

## Two replication levels

Every figure here is written twice, and the two files differ in one thing only:
what counts as a replicate.

- `_lineMeans` — the **maternal line** is the replicate. 3 points per regime per
  group (2 for Oregon 2014). Runs off `PS2_all.csv` / `GR_all.csv`.
- `_plants` — the **individual plant** is the replicate, nothing pre-averaged.
  9–15 points per regime per group. Runs off `PS2_all_plants.csv` /
  `GR_all_plants.csv`.

Same traits, same models, same panels, same code — only n changes. The title at
the top of each page states which one it is.

The line-level version is the conservative one and is what should carry any claim
in the text: plants within a maternal line are siblings grown side by side, so
they are not independent replicates of the line's response, and treating them as
such inflates the degrees of freedom without adding information about the line.
The plant-level version is here to show how much of the apparent signal is real
structure and how much is the line-level design not having the power to resolve
it. A contrast significant at plant level and absent at line level is a
candidate, not a result.

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

### How the newer scripts differ from `4`–`7`

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
4. **Layout.** `4.1` is one panel, not two: the single-timepoint regression that
   used to be a full-width panel A is now four small insets, each in the corner its
   own quadrant points into. The insets are squares of `INSET_SIZE` × the ring
   panel; 0.19 is the largest that geometrically cannot touch the ring, and the
   derivation is in a comment at the top of the script. The quadrant names are set
   along the arc, centred on each quadrant's own block of traits, and the cell each
   inset comes from is outlined in red with a red radial pointer running out to it
   from the centre.

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

The two `Figure4_TreeRing_*.pdf` and the eight `Figure5_*` / `Figure6_*` files in
`paper_graphs/`, and the four `FigS*.pdf` under
`paper_graphs/Supplemental/magda/`, are the ggplot output of the scripts above and
will reproduce exactly.

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
responds at all: in Figure 4, 902 of 2912 trait × group × timepoint regressions
pass FDR < 0.05, but California 2014 accounts for 465 of them (13 of 14 traits)
against 88 for Oregon 2010 (4 traits). Part of the Oregon 2014 shortfall is power —
2 lines, so 10 points per regression rather than 15.
