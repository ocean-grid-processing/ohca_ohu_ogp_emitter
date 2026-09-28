# ohca_ohu_ogp_emitter

`ohca_ohu_ogp_emitter` packages one `ogp_derive` blob into the annual **OHCA** (ocean heat content
anomaly) and **OHU** (ocean heat uptake) deliverable — one NetCDF per synthetic level,
`ohca_ohu_<tag>_<lo>_<hi>_dbar_<data>_tw<baseline>_<product_name>_<author>.nc`.

```
localgp_ogp_ingest ─▶ publish ─▶ ogp_derive (--quantities ohca,ohu,ohca_trend,ohu_trend) ─▶ ohca_ohu_ogp_emitter ─▶ per-level .nc
```

The analysis is all upstream. `ogp_derive` does the `n_fac` cross-layer combine, the annual means, the
OHCA baseline window, the OLS trends, and the ensemble → SD collapse. Its blob hands over `ohca`/`ohu`
(with their `_sd` and trends) as **basin-integrated extensive** quantities (TJ, and TJ per month),
plus `area_m2`, the `level`, and the `time_window` it was built with. This emitter is only the
packaging: divide by the area, carry the units to the target's per-area densities, relabel, and write.

This emitter is OHC-specific by design: it reads the `ohca`/`ohu` recipes and nothing else, so it has
no `[quantity]`-table generality to speak of. A different quantity gets a different emitter
(`map_ogp_emitter`, `mld_ogp_emitter`); this one exists to match one published target.

> **Units:** output matches the target (Zenodo 14720478 v4.0.0) — **`ohca` in J/m²**, **`ohu` in
> W/m²** (per-area densities; the `trend` attrs are W/m² on `ohca` and W/m²/s on `ohu`). Conversions:
> `ohca = ohca[TJ]/area × 1e12`; `ohu = ohu[TJ/mo]/area × 1e12 / sec_per_month`, using a **round
> 30-day month** (= 360-day year) to match the target (their OHU is 1.0146× a `365.25/12` month);
> trends divide by the seconds in one step of their `per` cadence (a 365-day year) to reach per-second.

## What it computes

Per blob (one synthetic level):

```
ohca(t) = blob.ohca / area_m2 × 1e12                    # J/m²
ohu(t)  = blob.ohu  / area_m2 × 1e12 / sec_per_month    # W/m²
*_std   = the matching blob.*_sd, same conversion       # present when the ensemble was on
```

The values ride a `time_ohca` axis — days since 2004-06-01, each year anchored at its 1-June — built
from the blob's `year` coord. The `low`/`high` in the filename come from the blob's `level` attr; the
`<data>` span from the blob's own `year` axis; the OHCA baseline label (`tw<baseline>`, and the
`ohca` long_name) from its `time_window` — a windowless derive run labels the baseline with the data
span itself.

**Trends.** `ohca_trend`/`ohu_trend` (and their `_sd`) are hung on the `ohca`/`ohu` variables as
`trend` / `trend_units` / `trend_std` attrs, converted per-second using each trend's `per` attr
(`year` → 365 d; `month` → 30 d): `trend` on `ohca` is W/m², on `ohu` is W/m²/s.

**OHU first year.** OHU is the annual mean of a month-to-month difference, whose leading step has no
prior month. `ogp_derive` voids that whole first year (NaN, not an 11-month partial mean) and keeps it
out of the trend fit; this emitter carries the NaN through unchanged, and it lands on disk as the
target's `-999` fill. Nothing is blanked here.

## Building the input

`ohca_ohu_ogp_emitter` consumes one `ogp_derive` blob per synthetic level, built along the LocalGP OHC
happy path in [`ogp_derive/examples/derive_ohc.slurm`](../ogp_derive/examples/derive_ohc.slurm):

```bash
python ../ogp_derive/run.py OHC_<tag>*.nc \
    --levels levels/localgp.toml --level 0_2000 --time-window 2005:2024 \
    --quantities ohca,ohu,ohca_trend,ohu_trend --mask contiguous_from_top \
    --bathy etopo60.cdf --tag <tag> --code-version URL \
    --product-name LocalGP --author Giglio_etal2026 --citation "…" --out <dir>
```

The `--time-window` is the OHCA baseline and the trend-fit years; it becomes this deliverable's
`tw<baseline>` token and `time_window` attr, so one derive run per baseline gives one deliverable per
baseline. Leave the ensemble on (the default) for the `_std` companions.

Each blob **must** carry:

- data vars **`ohca`** and **`ohu`** (annual, extensive; the emitter exits if either is absent);
  **`ohca_trend`** / **`ohu_trend`** (each with a `per` attr) for the trend attributes; and the `_sd`
  companions when the derive run kept the ensemble. Trend and `_sd` outputs that aren't present are
  simply skipped.
- attrs **`area_m2`**, **`level`**, and **`time_window`**.

Nothing else in the blob is read: the `quantity` table, `field_units` and `reduction` stamps ride
along inside the forwarded provenance but this emitter does not consult them.

## Usage

### Test
```bash
docker image build -t ohca_ohu_ogp_emitter:test .
docker container run -v $(pwd):/app ohca_ohu_ogp_emitter:test pytest
```

### Run

One blob in, one deliverable out, per level. The happy path is [`emit.slurm`](emit.slurm), which
takes `<baseline_window> <level>` and globs the matching derive blob out of the results directory;
[`run.sh`](run.sh) loops it over the windows and levels of a release:

```bash
sbatch emit.slurm 2005_2024 0_2000
```

#### emit.py options

All configuration is on the command line — no env, no config file. Every resolved option lands in
`config_record` (under this stage's `run_config`), except `--citation`, which has its own attr.

| option | required | default | what it does |
|---|:--:|---|---|
| `derive_*.nc` (positional, 1+) | **yes** | | `ogp_derive` blobs, one per synthetic level (`derive_<tag>_<data>_tw<baseline>_<level>.nc`). Each must carry `ohca` and `ohu` |
| `--tag` | **yes** | | run token in the filename and the `provenance_tag` attr. Used verbatim; should match the tag the blob was derived under |
| `--code-version` | **yes** | | URL to the exact `ohca_ohu_ogp_emitter` code (commit/release); recorded as this stage's `code_version` inside `config_record` |
| `--product-name` | **yes** | | product_name string; first of the filename's trailing pair (whitespace-stripped, case preserved), a standalone top-level `product_name` attr, and recorded in `config_record` |
| `--author` | **yes** | | author string; last of the filename's trailing pair (e.g. `Giglio_etal2026`) and recorded in `config_record` |
| `--citation` | **yes** | | citation sentence; written to the standalone top-level `citation` attr (kept out of `config_record` so it isn't duplicated) |
| `--provenance-link` | | *(none)* | URL/path to the provenance record; written to the `provenance_link` attr |
| `--out` | | `.` | output directory (created if absent) |

## Output and provenance

Global attrs on each file: `level`, `time_window`, `provenance_tag`, `provenance_link` (when given),
`citation`, `product_name`, and one `config_record`.

**Provenance chain.** Each deliverable is built from one derive blob, so this step is a 1-in-1-out
courier: it rolls that blob's whole provenance chain forward (every `*_run_config` / `*_run_facts` /
`*_code_version` — the grouped `localgp_ingest_*` / `localgp_publish_*` and the `ohc_derive_*` blocks)
and folds in its own block, emitting the lot as **one** `config_record` attribute keyed by stage:

```
config_record = {
  "localgp_ingest":       {"run_config": {…}, "run_facts": {…}, "code_version": "…"},
  "localgp_publish":      {…},
  "ohc_derive":           {…},
  "ohc_ohca_ohu_emitter": {"run_config": {resolved args}, "run_facts": {level, time_window, area_m2, quantities_present, ensemble, source_blob}, "code_version": "…"}
}
```

The stage keys are the pipeline's provenance contract and predate the repo renames: `ohc_derive` and
`ohc_ohca_ohu_emitter` are written by `ogp_derive` and this emitter respectively, and a submission that
arrived via `--contract ME4OH` carries no `localgp_*` blocks at all — only what the derive stage
recorded (including its `inferred_config`).

*Why one attribute:* a dozen separate global attributes tips HDF5 into **dense (fractal-heap) attribute
storage**, whose exact layout some netcdf builds mis-read; a single attribute keeps the file at ≤ 8
global attributes, i.e. **compact** storage, which every reader handles. The global `provenance_tag` /
`provenance_link` stay separate (they're the run's discoverable identity).

The per-constituent `localgp_*` blocks are **DRY'd**: each `{15_20:{…}, 15_300:{…}, …}` fan-out is
factored into `{"shared": {common config}, "per_constituent": {only what differs}}`, and a fan-out
whose entries fully agree (e.g. a shared `code_version`) collapses to a bare value. Lossless and
reversible (a constituent's block is `shared` merged with its `per_constituent` entry), driven by the
`constituents` roster in `ohc_derive.run_facts` — so only genuine fan-outs are touched and a value like
`n_fac`, nested inside a non-fanned block, is never mistaken for one.

## Notes

- The **30-day-month** OHU conversion is the magic number inferred from earlier implementations of this pipeline (their OHU is a constant
  1.0146× a `365.25/12` month); `ohca` is unaffected.
- **Trends are attributes**, not variables, mirroring the target's per-variable `trend` attr.
- **Fill value `-999`** (the target's) on every data variable on write, so the voided OHU first year
  lands as `-999` on disk (and decodes back to NaN on read).
- **Annual means are calendar-year** upstream (`ogp_derive` groups by `time.year`), and are stamped
  here at 1-June, which is the target's own convention: on the shared `time_ohca` axis our `ohca` and
  `ohu` reproduce the target level files index-for-index, and a ±1-year shift breaks the match.
