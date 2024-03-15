import numpy as np
import pandas as pd
import os
import pickle
import datetime as DT
import argparse
from mats_utils.rawdata.read_data import read_MATS_data
# from mats_l2_processing.forward_model import prepare_profile
from mats_utils.geolocation.coordinates import col_heights
from matplotlib import pyplot as plt
from scipy.ndimage import binary_erosion


def get_args():
    parser = argparse.ArgumentParser(description="NLC identification and analysis tools",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data", type=str, default="verdec2d.pickle",
                        help="Pickle filename for L1b data. Fetch data from aws and pickle if not found.")
    parser.add_argument("--channel", type=str, default="UV2", help="Channel to retrieve.")
    parser.add_argument("--out", type=str, default="ssa_analysis", help="Name prefix for result files.")
    parser.add_argument("--start_time", type=int, nargs=5, default=[2023, 3, 31, 21, 0],
                        help="Start time for data set.")
    parser.add_argument("--stop_time", type=int, nargs=5, default=[2023, 3, 31, 22, 35],
                        help="Start time for data set.")
    parser.add_argument("--show_images", type=int, nargs="*", default=[], help="Image numbers to plot.")
    parser.add_argument("--rmean_hwidth", type=int, default=2, help="Half-width for running means.")
    parser.add_argument("--nlc_alts", type=float, nargs=2, default=[80, 90], help="NLC altitude range, km.")
    parser.add_argument("--bgr_alts", type=float, nargs=2, default=[58, 70], help="Background altitude range, km.")
    parser.add_argument("--tomo_padding", type=int, default=16, help="Number of images to use for tomography around each nlc.")

    return parser.parse_args()


def get_data(args):
    dftop = read_MATS_data(DT.datetime(*args.start_time), DT.datetime(*args.stop_time), level="1b", version="0.4")
    with open(args.data, 'wb') as handle:
        pickle.dump(dftop, handle)
    return dftop


def lc_str2float(local_time):
    res = []
    for lc in local_time:
        nums = lc.split(":")
        res.append(int(nums[0]) + int(nums[1]) / 60 + int(nums[2]) / 3600)
    return np.array(res)


def altmean(image, tanalts, alt_range, percentile=None):
    BLANK_ROWS = 10
    assert image.shape[0] - BLANK_ROWS == tanalts.shape[0], f"{image.shape} {tanalts.shape}"
    assert len(alt_range) == 2
    thr_row = [np.argmin(np.abs(tanalts * 0.001 - val)) for val in alt_range]
    mdata = image[thr_row[0]:thr_row[1], :]
    if percentile is None:
        return np.mean(mdata)
    else:
        perc = np.percentile(mdata, percentile)
        return np.mean(np.minimum(mdata, perc))


def rmean(data, hwidth, cont_data, cont_thr):
    assert len(data.shape) == 1
    assert len(cont_data.shape) == 1
    size = len(data)
    assert size == len(cont_data)

    diffs = cont_data[1:] - cont_data[:-1]
    didx = [i + 1 for i in range(size - 1) if np.abs(diffs[i]) > cont_thr]
    didx = [0] + didx + [size]
    res = np.ones_like(data) * np.mean(data)
    width = 2 * hwidth + 1
    kernel = np.ones(width) / width

    for i in range(len(didx) - 1):
        if didx[i + 1] - didx[i] < width:
            continue
        res[(didx[i] + hwidth):(didx[i + 1] - hwidth)] = np.convolve(data[didx[i]:didx[i + 1]], kernel, mode="valid")
    return res, didx[1:-1]


def dual_plot(d1, name1, d2, name2, fname, vlines=[], shaded=None):
    assert len(d1.shape) == 1
    assert len(d2.shape) == 1
    size = len(d1)
    assert size == len(d2)
    numbers = list(range(size))
    plt.clf()
    fig, ax1 = plt.subplots(figsize=(24, 12))
    ax1.set_xlabel('Image number')
    ax1.set_ylabel(name1, color='tab:blue')
    ax1.plot(numbers, d1, color='tab:blue', lw=1)
    ax2 = ax1.twinx()
    ax2.set_ylabel(name2, color='tab:red')
    ax2.plot(numbers, d2, color='tab:red', lw=1)
    d2_range = (np.min(d2), np.max(d2))
    for vline in vlines:
        ax2.axvline(x=vline, color='k', lw=1)
    if shaded is not None:
        assert len(shaded) == size
        d2_range = (np.min(d2), np.max(d2))
        shade = np.where(shaded, d2_range[1], d2_range[0])
        ax2.fill_between(list(range(len(d2))), shade, y2=d2_range[0], alpha=0.4)
    plt.savefig(fname, dpi=400)


def pvgrad(data, denoise_iter=0, roll=0):
    pad = np.zeros((1, data.shape[1]))
    diff = np.maximum(0, np.diff(data, axis=0, append=pad))
    if denoise_iter > 0:
        diff[binary_erosion(diff > 0, iterations=denoise_iter) == 0] = 0.0
    if roll > 0:
        diff = np.minimum(diff, np.roll(diff, roll))
    # denoised = np.minimum(diff, np.roll(diff, roll, axis=1)) if roll > 0 else diff
    return diff


# Get data
args = get_args()
dftop = pd.read_pickle(args.data) if os.path.isfile(args.data) else get_data(args)

# Technical data
df = dftop[dftop['channel'] == args.channel].dropna().reset_index(drop=True)
columns = np.arange(0, df["NCOL"][0], 1)
ref_col = columns[int(len(columns) / 2)]
rows = np.arange(0, df["NROW"][0] - 10, 1)
num_images = len(df)

# thr_row = np.zeros((4, num_images), dtype=int)
for i in range(num_images):
    tanalts = np.array(col_heights(df.iloc[i], ref_col, 10, spline=True)(rows))
    # _, tanalts = prepare_profile(df.iloc[i], ref_col, rows)
    # thr_row[:, i] = [np.argmin(np.abs(tanalts * 0.001 - val)) for val in args.bgr_alts + args.nlc_alts]

halt_range = (90, 100)
# Main stats
ssa = df["TPssa"].to_numpy()
sza = df["TPsza"].to_numpy()
tplat = df["TPlat"].to_numpy()
tplc = lc_str2float(df["TPlocaltime"])
# mean_signal = np.array([np.mean(df["ImageCalibrated"][i][thr_row[0, i]:thr_row[1, i], :]) for i in range(num_images)])
# hor_grad = np.array([np.mean(np.diff(df["ImageCalibrated"][i][thr_row[0, i]:thr_row[1, i], :], axis=1) ** 2) for i in range(num_images)])
# nlc_alt_signal = np.array([np.mean(df["ImageCalibrated"][i][thr_row[2, i]:thr_row[3, i], :]) for i in range(num_images)])

halt_signal = np.array([altmean(df["ImageCalibrated"][i], tanalts, halt_range, percentile=99) for i in range(num_images)])
corrected = np.array([df["ImageCalibrated"][i]  for i in range(num_images)])

bgr_signal = np.array([altmean(corrected[i], tanalts, args.bgr_alts, percentile=99) for i in range(num_images)])
nlc_signal = np.array([altmean(corrected[i], tanalts, args.nlc_alts, percentile=99) for i in range(num_images)])
vgrad_nlc = np.array([altmean(pvgrad(corrected[i], denoise_iter=1, roll=1), tanalts, args.nlc_alts, percentile=99) for i in range(num_images)])

padding_kernel = np.ones(2 * args.tomo_padding + 1)
for_tomo = np.convolve(vgrad_nlc, padding_kernel, mode='same') > 0
nlc_id = vgrad_nlc > 0


bgr_signal_rmean, cutoffs = rmean(bgr_signal, args.rmean_hwidth, ssa, 0.1)
nlc_signal_rmean, _ = rmean(nlc_signal, args.rmean_hwidth, ssa, 0.1)

# bgr_signal_rmean = np.convolve(bgr_signal, np.ones(5) / 5, mode='valid')
numbers = list(range(num_images))

# Plots
#dual_plot(bgr_signal, "Background signal", nlc_signal, "NLC signal", f"{args.out}_means.png")
#dual_plot(bgr_signal_rmean, "Background signal", nlc_signal_rmean, "NLC signal", f"{args.out}_rmeans.png")
dual_plot(bgr_signal, "Rayleigh radience, rel. units", vgrad_nlc, "NLC proxy, rel. units", f"{args.out}_vgrad.png",
          vlines=cutoffs, shaded=for_tomo)
dual_plot(bgr_signal, f"Mean radiance {args.bgr_alts[0]}-{args.bgr_alts[1]} km altitude, rel. units",
          halt_signal, f"Mean radience {halt_range[0]}-{halt_range[1]} km  altitude, rel. units",
          f"{args.out}_bgr_halt.png", vlines=cutoffs)
dual_plot(nlc_signal, f"Mean radiance {args.nlc_alts[0]}-{args.nlc_alts[1]} km altitude, rel. units",
          bgr_signal, f"Mean radience {args.bgr_alts[0]}-{args.bgr_alts[1]} km  altitude, rel. units",
          f"{args.out}_nlc_bgr.png", vlines=cutoffs, shaded=nlc_id)
# fig, ax1 = plt.subplots()
# ax1.set_xlabel('Image number')
# ax1.set_ylabel('Rayleigh scattering ssa factor', color='tab:blue')
# ax1.plot(numbers, 1 - np.cos(np.deg2rad(df["TPlat"])), color='tab:blue')
# ax1.set_ylim(0, 1000)
#ax2 = ax1.twinx()
# ax2.set_ylabel(f'Mean signal', color='tab:red')
# ax2.plot(numbers, bgr_signal, color='tab:red')
# plt.savefig(f"{args.out}_means.png", dpi=400)

plt.figure()
thr = 1e-5
plt.scatter(halt_signal[vgrad_nlc > thr], nlc_signal[vgrad_nlc > thr],
            color="r", s=2, marker='+', alpha=0.4, label="NLC")
plt.scatter(halt_signal[vgrad_nlc < thr], nlc_signal[vgrad_nlc < thr],
            color="k", s=2, marker='+', alpha=0.4, label="no NLC")
plt.xlabel(f"Mean radience {halt_range[0]}-{halt_range[1]} km  altitude, rel. units")
plt.ylabel(f"Mean radiance {args.nlc_alts[0]}-{args.nlc_alts[1]} km altitude, rel. units")
plt.legend()
plt.savefig(f"halt_vs_nlc_{args.out}.png", dpi=400)


plt.figure()
plt.scatter(ssa, bgr_signal)
plt.xlabel("Solar scattering angle, deg")
plt.ylabel("Background signal, deg")
plt.savefig(f"ssa_vs_background_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(sza, bgr_signal)
plt.xlabel("Solar zenith angle, deg")
plt.ylabel("Background signal, deg")
plt.savefig(f"sza_vs_background_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(tplc, bgr_signal)
plt.xlabel("Local time, h")
plt.ylabel("Background signal, deg")
plt.savefig(f"lc_vs_background_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(ssa, nlc_signal)
plt.xlabel("Solar scattering angle, deg")
plt.ylabel("Signal from NLC altitudes, deg")
plt.savefig(f"ssa_vs_nlc_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(tplat, bgr_signal, c=tplc, cmap="twilight_shifted", s=1)
cb = plt.colorbar()
cb.set_label("Local time, h")
plt.xlabel("Latitude, deg")
plt.ylabel("Background signal, deg")
plt.savefig(f"lat_vs_mean_{args.out}.png", dpi=400)

# plt.figure()
# plt.scatter(tplat, bgr_signal_rmean, c=tplc, cmap="twilight_shifted", s=1)
# cb = plt.colorbar()
# cb.set_label("Local time, h")
# plt.xlabel("Latitude, deg")
# plt.ylabel("Background signal, deg")
# plt.savefig(f"lat_vs_rmean_{args.out}.png", dpi=400)
# plt.figure()

plt.scatter(ssa, sza, c=bgr_signal, cmap="inferno", marker="+", s=3)
cb = plt.colorbar()
cb.set_label("Mean signal")
plt.xlabel("Solar scattering angle, deg")
plt.ylabel("Solar zenith angle, deg")
plt.savefig(f"ssa_vs_sza_vs_bgr_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(ssa, bgr_signal, c=sza, cmap="inferno", marker="+", s=3)
cb = plt.colorbar()
cb.set_label("Solar zenith angle, deg")
plt.xlabel("Solar scattering angle, deg")
plt.ylabel("Mean signal")
plt.savefig(f"ssa_vs_sza_vs_bgr2_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(sza, bgr_signal, c=ssa, cmap="inferno", marker="+", s=3)
cb = plt.colorbar()
cb.set_label("Solar scattering angle, deg")
plt.xlabel("Solar zenith angle, deg")
plt.ylabel("Mean signal")
plt.savefig(f"ssa_vs_sza_vs_bgr3_{args.out}.png", dpi=400)

plt.figure()
plt.scatter(ssa, bgr_signal, c=nlc_signal, cmap="inferno", marker="+", s=3)
cb = plt.colorbar()
cb.set_label("Signal from NLC altitudes, deg")
plt.xlabel("Solar scattering angle, deg")
plt.ylabel("Mean signal")
plt.savefig(f"ssa_vs_bgr_vs_nlc_{args.out}.png", dpi=400)

# Plot images
for i in args.show_images:
    plt.figure()
    plt.clf()
    plt.imshow(df["ImageCalibrated"][i], cmap="inferno")
    plt.colorbar()
    plt.savefig(f"raw_{args.out}_{i}.png", dpi=400)
