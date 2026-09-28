#!/usr/bin/env python3
"""OHCA/OHU packaging: one ogp_derive blob -> the target per-area deliverable.

The factory has already done the analysis — the n_fac cross-layer combine, the annual means, the OHCA
baseline window, and the OLS trends. Its blob carries `ohca`/`ohu` (with their `_sd` and trends) as
basin-integrated extensive quantities (TJ, and TJ per month), plus `area_m2`, the `level`, and the
`time_window` it was built with. This step is only the packaging: divide by the area, carry the units
to the target's per-area densities, relabel, and write one file per level.

Target (Zenodo 14720478 v4.0.0): `ohca` in J/m2 and `ohu` in W/m2 on a `time_ohca` axis (days since
2004-06-01, each year anchored at its 1-June), with the linear trends as attrs on each variable.

OHU's first year arrives NaN from the factory (its leading tendency step is undefined, so that year is
voided upstream and excluded from the trend); it carries through here and lands as -999 on write.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import xarray as xr

EPOCH = "2004-06-01"        # time_ohca reference; each year anchored at its 1-June
TERA = 1e12                 # TJ -> J
# OHU per-month -> per-second: a round 30-day month (a 360-day year), matching the target's convention.
SEC_PER_MONTH = 30.0 * 86400.0
# Trend per year-step -> per-second: a 365-day year, matching the target's trend axis.
SEC_PER_YEAR = 365.0 * 86400.0
SEC_PER = {"year": SEC_PER_YEAR, "month": SEC_PER_MONTH}

# This step's identity, used to key its block inside the consolidated `config_record`. Each output file
# is built from one derive blob (one level), so this step is a 1-in-1-out courier: it rolls the blob's
# whole provenance chain forward and folds its own block in, emitting the lot as one `config_record`
# attribute (one attribute keeps the file in HDF5 compact storage — see `stamp_config_record`).
STAGE = "ohc_ohca_ohu_emitter"
_PROV_SUFFIXES = ("_run_config", "_run_facts", "_code_version")


def _compact(obj):
    """One-line JSON — reads as a single clean line in `ncdump -h`."""
    return json.dumps(obj, separators=(",", ":"), default=str)


def _shared_and_per(group_map):
    """{group: block} -> (shared, per): keys present in every group with an equal value go to `shared`;
    everything else stays per group. Lossless — block[g] == {**shared, **per[g]}."""
    groups = list(group_map)
    common = set(group_map[groups[0]])
    for g in groups[1:]:
        common &= set(group_map[g])
    shared = {}
    for k in sorted(common):
        vals = [group_map[g][k] for g in groups]
        if all(v == vals[0] for v in vals):
            shared[k] = vals[0]
    per = {g: {k: v for k, v in group_map[g].items() if k not in shared} for g in groups}
    return shared, per


def _compact_block(block, axis):
    """Factor one fan-out `{group: value}`: object values -> shared + per_<axis> (a fully-shared block
    collapses to the bare shared object); scalar values -> the bare value if all agree, else per_<axis>."""
    values = list(block.values())
    if all(isinstance(v, dict) for v in values):
        shared, per = _shared_and_per(block)
        if not any(per.values()):
            return shared
        return {"shared": shared, "per_" + axis: per}
    if all(v == values[0] for v in values):
        return values[0]
    return {"per_" + axis: dict(block)}


def _maybe_json(v):
    """Parse a forwarded block back to JSON so it nests as a real object; leave non-JSON (a bare
    code_version URL) as-is."""
    try:
        return json.loads(v)
    except (TypeError, ValueError):
        return v


def _dry_constituent_fanouts(record):
    """Within the assembled record, factor any per-constituent fan-out into shared + per_constituent
    (lossless). Driven by the constituents roster in `ohc_derive.run_facts`, so only genuine fan-outs
    are touched — a value like `n_fac`, nested inside a non-fanned block, is never a candidate."""
    facts = record.get("ohc_derive", {}).get("run_facts")
    roster = (set(facts["constituents"]) if isinstance(facts, dict)
              and isinstance(facts.get("constituents"), list) else None)
    if not roster:
        return record
    for parts in record.values():
        if not isinstance(parts, dict):
            continue
        for name, val in list(parts.items()):
            if isinstance(val, dict) and len(val) > 1 and set(val) <= roster:
                parts[name] = _compact_block(val, "constituent")
    return record


def stamp_config_record(out, blob, cfg, source_path):
    """Assemble the whole provenance chain into ONE `config_record` attribute, keyed by stage and DRY'd
    per constituent. A single attribute keeps the file at <=8 global attributes, i.e. HDF5 *compact*
    attribute storage — which every reader handles. Emitting a dozen separate `*_run_config` etc. tips
    HDF5 into dense (fractal-heap) storage, whose exact layout some netcdf builds mis-read."""
    record = {}
    # the forwarded chain: group the blob's *_run_config/_run_facts/_code_version by stage
    for k, v in blob.attrs.items():
        for suffix in _PROV_SUFFIXES:
            if k.endswith(suffix):
                record.setdefault(k[:-len(suffix)], {})[suffix[1:]] = _maybe_json(v)
                break
    # this step's own block. citation has its own top-level attr, so keep it out of the brick (not
    # duplicated); product_name/author stay in run_config for the record.
    record[STAGE] = {
        "run_config": {k: v for k, v in vars(cfg).items() if k != "citation"},
        "run_facts": {
            "level": blob.attrs.get("level"),
            "time_window": blob.attrs.get("time_window", "all"),
            "area_m2": float(blob.attrs["area_m2"]),
            "quantities_present": [q for q in ("ohca", "ohu", "ohca_trend", "ohu_trend") if q in blob],
            "ensemble": any(q + "_sd" in blob for q in ("ohca", "ohu")),
            "source_blob": os.path.abspath(source_path),
        },
        "code_version": cfg.code_version,
    }
    out.attrs["config_record"] = _compact(_dry_constituent_fanouts(record))


def _trend_seconds(trend_var):
    """Seconds in one step of a trend's cadence, read from its `per` attr (default year)."""
    return SEC_PER[trend_var.attrs.get("per", "year")]


def to_jm2(tj, area):
    """Basin-integrated TJ -> per-area J/m2."""
    return tj / area * TERA


def to_wm2(tj_per_month, area):
    """Basin-integrated TJ per month -> per-area W/m2."""
    return tj_per_month / area * TERA / SEC_PER_MONTH


def _time_ohca(years):
    """Integer years -> days since EPOCH, each year anchored at its 1-June."""
    ref = pd.Timestamp(EPOCH)
    days = np.array([(pd.Timestamp("%d-06-01" % y) - ref).days for y in years], dtype="float64")
    return xr.DataArray(days, dims=("time_ohca",),
                        attrs={"units": "days since %s 00:00:00" % EPOCH,
                               "calendar": "proleptic_gregorian", "long_name": "time"})


def build_dataset(blob, tag, provenance_link, citation="", product_name=""):
    """A derive blob -> the OHCA/OHU deliverable Dataset over `time_ohca`."""
    area = float(blob.attrs["area_m2"])
    window = blob.attrs.get("time_window", "all")
    baseline = "all-time mean" if window == "all" else "%s mean" % window

    time = _time_ohca(blob["ohca"]["year"].values.astype("int64"))
    dv = {
        "ohca": xr.DataArray(to_jm2(blob["ohca"].values, area), dims=("time_ohca",),
                             attrs={"units": "J/m2", "area_m2": area,
                                    "long_name": "annual OHC anomaly (%s removed)" % baseline}),
        "ohu": xr.DataArray(to_wm2(blob["ohu"].values, area), dims=("time_ohca",),
                            attrs={"units": "W/m2", "area_m2": area,
                                   "long_name": "annual mean ocean heat uptake"}),
    }
    if "ohca_sd" in blob:
        dv["ohca_std"] = xr.DataArray(to_jm2(blob["ohca_sd"].values, area), dims=("time_ohca",),
                                      attrs={"units": "J/m2"})
        dv["ohu_std"] = xr.DataArray(to_wm2(blob["ohu_sd"].values, area), dims=("time_ohca",),
                                     attrs={"units": "W/m2"})

    # Linear trends as attrs, per second: the trend's `per` attr picks the seconds in one step.
    if "ohca_trend" in blob:
        sec = _trend_seconds(blob["ohca_trend"])
        dv["ohca"].attrs.update({"trend": to_jm2(float(blob["ohca_trend"]), area) / sec, "trend_units": "W/m2"})
        if "ohca_trend_sd" in blob:
            dv["ohca"].attrs["trend_std"] = to_jm2(float(blob["ohca_trend_sd"]), area) / sec
    if "ohu_trend" in blob:
        sec = _trend_seconds(blob["ohu_trend"])
        dv["ohu"].attrs.update({"trend": to_wm2(float(blob["ohu_trend"]), area) / sec, "trend_units": "W/m2/s"})
        if "ohu_trend_sd" in blob:
            dv["ohu"].attrs["trend_std"] = to_wm2(float(blob["ohu_trend_sd"]), area) / sec

    out = xr.Dataset(dv, coords={"time_ohca": time})
    out.attrs["level"] = blob.attrs["level"]
    out.attrs["time_window"] = window
    out.attrs["provenance_tag"] = tag
    if provenance_link is not None:
        out.attrs["provenance_link"] = provenance_link
    out.attrs["citation"] = citation
    if product_name:
        out.attrs["product_name"] = product_name   # top-level discoverable key (also in config_record)
    return out


def _data_span(blob):
    """`YYYY_YYYY` for the years the blob's own axis spans (`year` annual, `time` monthly)."""
    if "year" in blob.coords:
        yrs = blob["year"].values.astype(int)
        return "%d_%d" % (int(yrs.min()), int(yrs.max()))
    if "time" in blob.coords:
        yrs = blob["time"].values.astype("datetime64[Y]").astype(int) + 1970
        return "%d_%d" % (int(yrs.min()), int(yrs.max()))
    return "all"


def _file_token(blob):
    """Combined filename token `<data>_tw<baseline>`: the blob's own data span, then its baseline
    window (defaulting to the whole data span when the derive run was windowless) — matching the derive
    filename, so runs differing only in baseline don't collide here either."""
    data = _data_span(blob)
    win = blob.attrs.get("time_window", "all")
    baseline = data if (not win or win == "all") else win.replace("-", "_")
    return "%s_tw%s" % (data, baseline)


def filename(level, tag, token, product_name, author):
    """Per-level name: ohca_ohu_<tag>_<lo>_<hi>_dbar_<data>_tw<baseline>_<product_name>_<author>.nc
    (tag leads after the step; product_name/author are the last thing before .nc)."""
    low, high = level.split("_")
    return "ohca_ohu_%s_%s_%s_dbar_%s_%s_%s.nc" % (tag, low, high, token, product_name, author)


def main():
    ap = argparse.ArgumentParser(description="OHCA/OHU packaging: ogp_derive blob -> target deliverable")
    ap.add_argument("blobs", nargs="+", help="ogp_derive output NetCDFs (derive_<tag>_<data>_tw<baseline>_<level>.nc)")
    ap.add_argument("--tag", required=True, help="provenance tag: filename token + provenance_tag attr")
    ap.add_argument("--provenance-link", default=None, help="URL/path to the provenance record")
    ap.add_argument("--code-version", required=True,
                    help="URL to the exact ohca_ohu_ogp_emitter code (commit/release); recorded as this "
                         "stage's code_version inside config_record")
    ap.add_argument("--product-name", required=True,
                    help="product_name string, the first of the filename's trailing pair and in config_record "
                         "(e.g. LocalGP)")
    ap.add_argument("--author", required=True,
                    help="author string, the last of the filename's trailing pair and in config_record "
                         "(e.g. Giglio_etal2026)")
    ap.add_argument("--citation", required=True,
                    help="citation sentence; written to the top-level `citation` attr")
    ap.add_argument("--out", default=".")
    cfg = ap.parse_args()
    cfg.product_name = "".join(cfg.product_name.split())                  # filename tokens: whitespace-stripped,
    cfg.author = "".join(cfg.author.split())                    # case preserved, no other munging
    os.makedirs(cfg.out, exist_ok=True)
    for path in cfg.blobs:
        blob = xr.open_dataset(path)
        if "ohca" not in blob or "ohu" not in blob:
            raise SystemExit("%s carries no ohca/ohu; run ogp_derive with --quantities ohca,ohu (+ trends)"
                             % path)
        dest = os.path.join(cfg.out, filename(blob.attrs["level"], cfg.tag, _file_token(blob),
                                              cfg.product_name, cfg.author))
        out = build_dataset(blob, cfg.tag, cfg.provenance_link, cfg.citation, cfg.product_name)
        stamp_config_record(out, blob, cfg, path)                  # whole chain -> one config_record attr
        enc = {v: {"_FillValue": -999.0} for v in out.data_vars}   # target fill (NaN -> -999)
        out.to_netcdf(dest, engine="netcdf4", encoding=enc)
        print("wrote", dest)


if __name__ == "__main__":
    main()
