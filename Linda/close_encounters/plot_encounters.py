# %% Imports and configuration
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CHANNELS = ["IR1", "IR2", "IR3", "IR4", "UV1", "UV2"]

# Same panel layout as quickview.py
PANEL_LAYOUT = {
    "IR1": (0, 0), "IR3": (0, 1), "UV1": (0, 2),
    "IR2": (1, 0), "IR4": (1, 1), "UV2": (1, 2),
}

CMAP = "magma"
TOL_S = 30.0  # accept nearest image within this many seconds of mats_obs_utc

CSV_PATH = Path("/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/close_encounters/out/criterion_D_tangent_in_cips_swath.csv")
OUT_DIR = Path("/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/close_encounters/out/plots")

# Local hourly L1b pickles produced by
# MATS-analysis/Linda/PlottingMonitoring/download_script.py (read_MATS_data, level='1b').
# Used instead of the S3 zarr store (mats-level-1b-limb-cropd), whose large ImageCalibrated
# chunks currently fail to download (server-side truncated responses).
LOCAL_DATA_DIR = Path("/Users/lindamegner/MATS/MATS-retrieval/data/daily_20230209-20230215")

# %% Load and filter events
df = pd.read_csv(CSV_PATH)
df["mats_obs_dt"] = pd.to_datetime(df["mats_obs_utc"], utc=True)

t0 = pd.Timestamp("2023-02-11", tz="UTC")
t1 = pd.Timestamp("2023-02-17 23:59:59", tz="UTC")
mask = (
    (df["mats_obs_dt"] >= t0)
    & (df["mats_obs_dt"] <= t1)
    & (df["mats_tp_lat_deg"] < -75)
)
filtered = df[mask].reset_index(drop=True)
print(f"Events matching criteria (Feb 11-17, TP lat < -75°): {len(filtered)}")
print(filtered[["mats_obs_utc", "mats_tp_lat_deg", "mats_tp_lon_deg",
                "distance_tp_to_aim_km"]].to_string(index=False))

# %% Index local hourly pickle files (df_<start>_<end>.pkl, one hour of all channels each)
_HOUR_FILE_RE = re.compile(r"df_(\d{8}_\d{6})_(\d{8}_\d{6})\.pkl$")


def _index_hour_files(data_dir: Path) -> list[tuple[pd.Timestamp, pd.Timestamp, Path]]:
    index = []
    for p in sorted(data_dir.glob("df_*.pkl")):
        m = _HOUR_FILE_RE.search(p.name)
        if not m:
            continue
        start = pd.to_datetime(m.group(1), format="%Y%m%d_%H%M%S", utc=True)
        end = pd.to_datetime(m.group(2), format="%Y%m%d_%H%M%S", utc=True)
        index.append((start, end, p))
    return index


HOUR_FILES = _index_hour_files(LOCAL_DATA_DIR)
if not HOUR_FILES:
    raise SystemExit(f"No hourly pickle files found in {LOCAL_DATA_DIR}")
print(f"Indexed {len(HOUR_FILES)} local hourly files spanning "
      f"{HOUR_FILES[0][0]} .. {HOUR_FILES[-1][1]}")

_hour_cache: dict[Path, pd.DataFrame] = {}


def _load_hour_file(path: Path) -> pd.DataFrame:
    """Load one hourly pickle, keeping only the columns plotting needs (cached)."""
    if path not in _hour_cache:
        print(f"    loading {path.name} ...", end=" ", flush=True)
        raw = pd.read_pickle(path)
        _hour_cache[path] = raw[["TMHeaderTime", "channel", "TPlat", "TPlon", "ImageCalibrated"]].copy()
        print("ok")
    return _hour_cache[path]


# %% Helper: find nearest image
def nearest_image(channel: str, target: pd.Timestamp, tol_s: float) -> pd.Series | None:
    """Return the local L1b row for `channel` nearest to `target`, or None if beyond tol_s."""
    tol = pd.Timedelta(seconds=tol_s)
    lo, hi = target - tol, target + tol
    candidate_files = [p for (start, end, p) in HOUR_FILES if start <= hi and end >= lo]

    best_row, best_d = None, tol
    for path in candidate_files:
        sub = _load_hour_file(path)
        sub = sub[sub["channel"] == channel]
        if sub.empty:
            continue
        diffs = (sub["TMHeaderTime"] - target).abs()
        idx = diffs.idxmin()
        d = diffs.loc[idx]
        if d <= best_d:
            best_d = d
            best_row = sub.loc[idx]
    return best_row


# %% Plot all events and save PNGs
OUT_DIR.mkdir(parents=True, exist_ok=True)

for i, row in filtered.iterrows():
    obs_time = pd.Timestamp(row["mats_obs_utc"])
    if obs_time.tzinfo is None:
        obs_time = obs_time.tz_localize("UTC")
    t_str = obs_time.strftime("%Y-%m-%dT%H:%M:%SZ")
    tp_lat = row["mats_tp_lat_deg"]
    tp_lon = row["mats_tp_lon_deg"]
    aim_lat = row["aim_center_lat_deg"]
    aim_lon = row["aim_center_lon_deg"]
    dt_min = row["time_diff_min"]

    print(f"[{i + 1}/{len(filtered)}] {t_str}  TP lat={tp_lat:.2f}°")

    fig, axes = plt.subplots(2, 3, figsize=(13, 6))
    fig.suptitle(
        f"MATS × CIPS encounter  |  MATS obs: {t_str}\n"
        f"MATS TP lat={tp_lat:.2f}°  lon={tp_lon:.2f}°  |  "
        f"AIM center lat={aim_lat:.2f}°  lon={aim_lon:.2f}°  |  "
        f"Δt={dt_min:+.1f} min",
        fontsize=10,
    )

    for ch, (ri, ci) in PANEL_LAYOUT.items():
        ax = axes[ri, ci]
        entry = nearest_image(ch, obs_time, TOL_S)

        if entry is None:
            ax.text(0.5, 0.5, f"{ch}\n(no data within\n{TOL_S:.0f} s)",
                    ha="center", va="center", transform=ax.transAxes, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        img = np.asarray(entry["ImageCalibrated"])
        mean = float(np.nanmean(img))
        std = float(np.nanstd(img))
        im = ax.imshow(img, origin="lower", aspect="auto", cmap=CMAP,
                       vmin=mean - 2 * std, vmax=mean + 2 * std)

        tp_lat_img = float(entry["TPlat"])
        tp_lon_img = float(entry["TPlon"])
        img_t = entry["TMHeaderTime"].strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
        ax.set_title(
            f"{ch}\nTMHeaderTime={img_t}\nTP=({tp_lat_img:.2f}°, {tp_lon_img:.2f}°)",
            fontsize=7,
        )
        ax.set_xticks([])
        ax.set_yticks([])

        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label("photon nm$^{-1}$ m$^{-2}$ sr$^{-1}$ s$^{-1}$", fontsize=7, rotation=270, labelpad=10)
        cb.ax.tick_params(labelsize=7)

    fig.tight_layout()
    fname = OUT_DIR / f"encounter_{t_str.replace(':', '').replace('-', '')}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {fname.name}")

print(f"\nDone. {len(filtered)} plot(s) saved to {OUT_DIR}/")

# %%
