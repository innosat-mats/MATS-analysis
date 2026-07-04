"""
Overview plots for the locally extracted MATS Level0 payload data
(produced by rac-extract-payload, read via read_MATS_payload_data_from_disk).

Generates:
  1. HTR heater temperatures over time
  2. PWR bus voltages and currents over time
  3. PM photometer counts over time
  4. CCD sensor temperature and exposure time over time
  5. Data coverage: packet counts per hour per subsystem
  6. Sample CCD images, one per CCDSEL
"""

import datetime as DT
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mats_utils.rawdata.read_data import read_MATS_payload_data_from_disk
from mats_l1_processing.read_parquet_functions import convert_image_data, remove_faulty_rows

EXTRACTED_PATH = str(Path(__file__).resolve().parents[2] / "extracted")
START = DT.datetime(2026, 7, 1)
STOP = DT.datetime(2026, 7, 4)
OUTDIR = Path(__file__).resolve().parent

SUBSYSTEMS = ["HTR", "PWR", "PM", "CPRU", "TCV", "STAT", "CCD"]


def load_all():
    data = {}
    for name in SUBSYSTEMS:
        try:
            data[name] = read_MATS_payload_data_from_disk(
                START, STOP, data_type=name, path=EXTRACTED_PATH
            )
            print(f"{name}: {len(data[name])} rows")
        except Exception as err:
            print(f"{name}: failed to load ({err})")
    return data


def plot_htr_temperatures(df):
    fig, ax = plt.subplots(figsize=(11, 5))
    heaters = ["HTR1", "HTR2", "HTR7", "HTR8"]
    channels = ["A", "B", "OD"]
    colors = plt.cm.tab10(np.linspace(0, 1, len(heaters)))
    for color, heater in zip(colors, heaters):
        for channel, style in zip(channels, ["-", "--", ":"]):
            col = f"{heater}{channel}"
            if col in df.columns:
                ax.plot(df["TMHeaderTime"], df[col], style, color=color,
                        lw=1, alpha=0.8, label=col)
    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("HTR heater temperatures")
    ax.legend(fontsize=7, ncol=4)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTDIR / "01_htr_temperatures.png", dpi=150)


def plot_pwr_bus(df):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for col in ["PWRP32V", "PWRP16V", "PWRM16V", "PWRP3V3"]:
        axes[0].plot(df["TMHeaderTime"], df[col], lw=1, label=col)
    axes[0].set_ylabel("Voltage (V)")
    axes[0].set_title("PWR bus voltages")
    axes[0].legend(fontsize=8)

    for col in ["PWRP32C", "PWRP16C", "PWRM16C", "PWRP3C3"]:
        axes[1].plot(df["TMHeaderTime"], df[col], lw=1, label=col)
    axes[1].set_ylabel("Current (A)")
    axes[1].set_xlabel("Time")
    axes[1].set_title("PWR bus currents")
    axes[1].legend(fontsize=8)

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTDIR / "02_pwr_bus.png", dpi=150)


def plot_pm_counts(df):
    fig, ax = plt.subplots(figsize=(11, 5))
    for col in ["PM1A", "PM1B", "PM1S", "PM2A", "PM2B", "PM2S"]:
        ax.plot(df["TMHeaderTime"], df[col], lw=1, alpha=0.8, label=col)
    ax.set_xlabel("Time")
    ax.set_ylabel("Counts")
    ax.set_title("PM photometer counts")
    ax.legend(fontsize=8, ncol=3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTDIR / "03_pm_counts.png", dpi=150)


def plot_ccd_housekeeping(df):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for ccdsel, sub in df.groupby("CCDSEL"):
        axes[0].plot(sub["TMHeaderTime"], sub["TEMP"], ".", ms=3, label=f"CCDSEL {ccdsel}")
        axes[1].plot(sub["TMHeaderTime"], sub["TEXPMS"], ".", ms=3, label=f"CCDSEL {ccdsel}")
    axes[0].set_ylabel("CCD TEMP (raw)")
    axes[0].set_title("CCD sensor temperature")
    axes[0].legend(fontsize=7, ncol=4)
    axes[1].set_ylabel("Exposure time (ms)")
    axes[1].set_xlabel("Time")
    axes[1].set_title("CCD exposure time")

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTDIR / "04_ccd_housekeeping.png", dpi=150)


def plot_data_coverage(data):
    fig, ax = plt.subplots(figsize=(11, 5))
    for name, df in data.items():
        if len(df) == 0:
            continue
        counts = df.set_index("TMHeaderTime").resample("1h").size()
        ax.plot(counts.index, counts.values, drawstyle="steps-mid", label=name)
    ax.set_xlabel("Time")
    ax.set_ylabel("Packets per hour")
    ax.set_title("Data coverage per subsystem")
    ax.legend(fontsize=8, ncol=4)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTDIR / "05_data_coverage.png", dpi=150)


def plot_sample_images(df):
    df = df.copy()
    convert_image_data(df)
    df = remove_faulty_rows(df)

    ccdsels = sorted(df["CCDSEL"].unique())
    ncols = 4
    nrows = int(np.ceil(len(ccdsels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, ccdsel in zip(axes, ccdsels):
        row = df[df["CCDSEL"] == ccdsel].iloc[0]
        im = ax.imshow(row["IMAGE"], cmap="viridis", aspect="auto")
        ax.set_title(f"CCDSEL {ccdsel}\n{row['TMHeaderTime']}", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046)
    for ax in axes[len(ccdsels):]:
        ax.axis("off")

    fig.suptitle("Sample CCD image per CCDSEL")
    fig.tight_layout()
    fig.savefig(OUTDIR / "06_ccd_sample_images.png", dpi=150)


def main():
    data = load_all()

    if "HTR" in data and len(data["HTR"]):
        plot_htr_temperatures(data["HTR"])
    if "PWR" in data and len(data["PWR"]):
        plot_pwr_bus(data["PWR"])
    if "PM" in data and len(data["PM"]):
        plot_pm_counts(data["PM"])
    if "CCD" in data and len(data["CCD"]):
        plot_ccd_housekeeping(data["CCD"])
        plot_sample_images(data["CCD"])
    plot_data_coverage(data)

    print(f"Saved plots to {OUTDIR}")
    plt.show()


if __name__ == "__main__":
    main()
