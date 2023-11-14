import numpy as np
import sys
import pickle
from scipy.interpolate import CubicSpline
import pandas as pd

# Background subtraction script for tomography data

NCOLS = 43
REF_ALT = [9e4, 1e5]

def interp_image(rads, tanalts, ref_tanalts):
    assert rads.shape == tanalts.shape
    assert len(ref_tanalts.shape) == 1
    assert ref_tanalts.shape[0] == tanalts.shape[1]
    res = np.zeros_like(rads)
    for i in range(rads.shape[0]):
        profile = CubicSpline(tanalts[i, :], rads[i, :], extrapolate=True)
        res[i, :] = profile(ref_tanalts)
    return res


ref = [int(x) for x in sys.argv[3].split(",")]
jacobian = pd.read_pickle(sys.argv[1])
tan_alts = pd.read_pickle(f"tanalts_{sys.argv[1]}")
(y, _, _, _, _, _) = jacobian
y = y.reshape((-1, NCOLS, y.shape[1]))

tan_alts = tan_alts.reshape((-1, NCOLS, tan_alts.shape[1]))
ref_tanalts = tan_alts[ref[0], int(NCOLS / 2), :]
background = np.zeros_like(y[0, :, :])
for ir in ref:
    background += interp_image(y[ir, :, :], tan_alts[ir, :, :], ref_tanalts)
background /= float(len(ref))



for k in range(y.shape[1]):
    bcol = CubicSpline(ref_tanalts, background[k, :], extrapolate=False)
    for i in range(y.shape[0]):
        y[i, k, :] -= bcol(tan_alts[i, k, :])

y[np.isnan(y)] = 0.0
#y[y < 0.0] = 0.0

# "Calibration" by subtracting top altitude mean
#for k in range(y.shape[0]):
#    print(f"{np.mean(y[k, :, :]):.2E}",
#          f"{np.mean(y[k, :, :][np.logical_and(tan_alts[k, :, :] > REF_ALT[0], tan_alts[k, :, :] < REF_ALT[1])]):.2E}")
#    y[k, :, :] -= np.mean(y[k, :, :][np.logical_and(tan_alts[k, :, :] > REF_ALT[0], tan_alts[k, :, :] < REF_ALT[1])])
y = y.reshape((-1, y.shape[2]))
with open(f"y_{sys.argv[1]}_{sys.argv[2]}", "wb") as file:
    pickle.dump(y, file)
