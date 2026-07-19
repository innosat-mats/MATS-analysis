# %% Imports and configuration
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mats_l1b_tools.fetch_data import fetch_MATS_l1b_data

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

# %% Load and filter events
df = pd.read_csv(CSV_PATH)
df["mats_obs_dt"] = pd.to_datetime(df["mats_obs_utc"], utc=True)

t0 = pd.Timestamp("2023-02-11", tz="UTC")
t1 = pd.Timestamp("2023-02-13 23:59:59", tz="UTC")
mask = (
    (df["mats_obs_dt"] >= t0)
    & (df["mats_obs_dt"] <= t1)
    & (df["mats_tp_lat_deg"] < -75)
)
filtered = df[mask].reset_index(drop=True)
print(f"Events matching criteria (Feb 11-13, TP lat < -75°): {len(filtered)}")
print(filtered[["mats_obs_utc", "mats_tp_lat_deg", "mats_tp_lon_deg",
                "distance_tp_to_aim_km"]].to_string(index=False))

# %% Open all channels (run once — takes a moment)
datasets = {}
for ch in CHANNELS:
    print(f"  {ch} ...", end=" ", flush=True)
    datasets[ch] = fetch_MATS_l1b_data(ch)
    print("ok")


# %% Helper: find nearest image
def nearest_image(ds, target: np.datetime64, tol_s: float):
    """Return the image row whose time is nearest to target, or None if beyond tol_s."""
    times = ds["time"].values
    if len(times) == 0:
        return None
    tol = np.timedelta64(int(tol_s * 1e9), "ns")
    j = int(np.searchsorted(times, target))
    candidates = [j - 1, j]
    best, best_d = None, tol + np.timedelta64(1, "ns")
    for c in candidates:
        if 0 <= c < len(times):
            d = abs(times[c] - target)
            if d < best_d:
                best_d = d
                best = c
    if best is None:
        return None
    return ds.isel(time=best).load()


# %% Plot all events and save PNGs
OUT_DIR.mkdir(parents=True, exist_ok=True)

for i, row in filtered.iterrows():
    obs_time = np.datetime64(row["mats_obs_utc"].replace("Z", ""), "s").astype("datetime64[ns]")
    t_str = pd.Timestamp(obs_time).strftime("%Y-%m-%dT%H:%M:%SZ")
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
        entry = nearest_image(datasets[ch], obs_time, TOL_S)

        if entry is None:
            ax.text(0.5, 0.5, f"{ch}\n(no data within\n{TOL_S:.0f} s)",
                    ha="center", va="center", transform=ax.transAxes, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        img = np.asarray(entry["ImageCalibrated"].values)
        mean = float(np.nanmean(img))
        std = float(np.nanstd(img))
        im = ax.imshow(img, origin="lower", aspect="auto", cmap=CMAP,
                       vmin=mean - 2 * std, vmax=mean + 2 * std)

        tp_lat = float(entry["TPlat"].values) if "TPlat" in entry else float("nan")
        tp_lon = float(entry["TPlon"].values) if "TPlon" in entry else float("nan")
        img_t = pd.Timestamp(entry["time"].values).strftime("%H:%M:%S")
        ax.set_title(f"{ch}  {img_t}  TP=({tp_lat:.2f}°, {tp_lon:.2f}°)", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])

        units = str(datasets[ch]["ImageCalibrated"].attrs.get("units", "")).strip()
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        if units:
            cb.set_label(units, fontsize=7, rotation=270, labelpad=10)
        cb.ax.tick_params(labelsize=7)

    fig.tight_layout()
    fname = OUT_DIR / f"encounter_{t_str.replace(':', '').replace('-', '')}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {fname.name}")

print(f"\nDone. {len(filtered)} plot(s) saved to {OUT_DIR}/")
