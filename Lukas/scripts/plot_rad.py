import numpy as np
import pandas as pd
import sys
from matplotlib import pyplot as plt
from mats_l2_processing.forward_model import prepare_profile
from scipy.interpolate import CubicSpline

# This script plots MATS images and optionally removes background by subtracting
# indicated images
# 
# Usage: python plot_rad.py <pickle file with MATS data> <channel>
#        [comma separated list of image numbers from which background will be constructed] 


def get_rads(rows, columns, data):
    rad = np.zeros((len(rows), len(columns)))
    tanalt = rad.copy()
    for i, col in enumerate(columns):
        profile, tanalts = prepare_profile(data, col, rows)
        rad[:, i] = profile
        tanalt[:, i] = tanalts
    return rad, tanalt


def interp_image(rads, tanalts, ref_tanalts):
    assert rads.shape == tanalts.shape
    assert len(ref_tanalts.shape) == 1
    assert ref_tanalts.shape[0] == tanalts.shape[0]
    res = np.zeros_like(rads)
    for i in range(rads.shape[1]):
        profile = CubicSpline(tanalts[:, i], rads[:, i], extrapolate=True)
        res[:, i] = profile(ref_tanalts)
    return res


kind = 'IR2' if len(sys.argv) < 3 else sys.argv[2]
ref = [] if len(sys.argv) < 4 else [int(x) for x in sys.argv[3].split(",")]
dftop = pd.read_pickle(sys.argv[1])
df = dftop[dftop['channel'] == kind].dropna().reset_index(drop=True)
columns = np.arange(0, df["NCOL"][0], 1)
rows = np.arange(0, df["NROW"][0] - 10, 1)


all_tanalts = np.zeros((len(df), len(rows), len(columns)))
all_rads = all_tanalts.copy()

for image in range(len(df)):
    all_rads[image, :, :], all_tanalts[image, :, :] = get_rads(rows, columns, df.iloc[image])

if len(ref) > 0:
    ref_tanalts = all_tanalts[ref[0], :, int(len(columns) / 2)]
    background = np.zeros_like(all_rads[0, :, :])
    for ir in ref:
        background += interp_image(all_rads[ir, :, :], all_tanalts[ir, :, :], ref_tanalts)
    background /= float(len(ref))

    for k in range(len(columns)):
        bcol = CubicSpline(ref_tanalts, background[:, k], extrapolate=False)
        for i in range(len(df)):
            all_rads[i, :, k] -= bcol(all_tanalts[i, :, k])

all_rads[all_rads < 0] = 0.0
vmax = 0.8 * np.nanmax(all_rads)

means = np.array([np.mean(all_rads[i, :, :][all_rads[i, :, :] > 1.0]) for i in range(all_rads.shape[0])])
nums = np.arange(means.shape[0])
plt.figure()
plt.plot(nums, means)
plt.xlabel("Image number")
plt.ylabel("Mean radiance")
plt.savefig("NLC_mean.png", dpi=400)

across = np.broadcast_to(columns[np.newaxis, :], all_tanalts.shape[1:])
for i in range(len(df)):
    # across = np.zeros_like(tanalt)
    # breakpoint()
    plt.figure()
    plt.pcolor(across, all_tanalts[i, :, :] * .001, all_rads[i, :, :], cmap="inferno", vmin=0, vmax=vmax)
    plt.xlabel("Image colums")
    plt.ylabel("Altitude, km")
    plt.colorbar()
    plt.savefig(f"nlc_noraileigh_{kind}_{i:04}.png", dpi=400)
    plt.close()
