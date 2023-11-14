import numpy as np
import pandas as pd
import pickle
# from matplotlib import pyplot as plt
import sys
from mats_l2_processing.grids import sph2cart, geoid_radius
from scipy.signal import savgol_filter
# from scipy import interpolate

# Tests for various tomography aptioris

ALONG_MARGIN_RATIO = 0.15
ALT_RANGE = (80000, 100000)
SIGMA = 0.075


def get_local_alts(radius_grid, acrosstrack_grid, alongtrack_grid, ecef_to_local):
    rr, acrr, alongg = np.meshgrid(radius_grid, acrosstrack_grid, alongtrack_grid, indexing="ij")
    lxx, lyy, lzz = sph2cart(rr, acrr, alongg)
    glgrid = ecef_to_local.inv().apply(np.dstack((lxx.flatten(), lyy.flatten(), lzz.flatten()))[0, :, :])
    altt = rr - geoid_radius(np.arcsin(glgrid[:, 2].reshape(rr.shape) / rr))
    return altt


def envelope(alts, min_alt, max_alt, rel_sigma):
    min_idx, max_idx = [np.argmin(np.abs(alts - v)) for v in [min_alt, max_alt]]
    edge_d = 0.5 - np.abs(np.linspace(-0.5, 0.5, max_idx - min_idx))
    res = np.zeros_like(alts)
    res[min_idx:max_idx] = 1.0 - np.exp(-(edge_d / rel_sigma) ** 2)
    # res /= np.mean(res)
    return res


def sg_apriori(res_pickle):
    x_hat, altitude_grid, alongtrack_grid, acrosstrack_grid, ecef_to_local = pd.read_pickle(res_pickle)

    alts = get_local_alts(altitude_grid, acrosstrack_grid, alongtrack_grid, ecef_to_local)
    margin = int(len(alongtrack_grid) * ALONG_MARGIN_RATIO)
    # alongtrack_grid = alongtrack_grid[margin:-margin]
    x = x_hat[:, :, margin:-margin]
    z = alts[:, :, margin:-margin]
    distr, zh = np.histogram(z, weights=x, bins=100, density=False)
    distr /= np.histogram(z, bins=100, density=False)[0]
    zh = 0.5 * (zh[1:] + zh[:-1])
    distr[distr < 0] = 0
    distr[zh < ALT_RANGE[0]] = 0
    distr[zh > ALT_RANGE[1]] = 0
    apr_profile = savgol_filter(distr, 25, 2, mode='nearest')
    apr_profile[apr_profile < 0] = 0
    apr_profile *= envelope(zh, ALT_RANGE[0], ALT_RANGE[1], SIGMA)
    # f = interpolate.interp1d(zh, apr_profile)
    return apr_profile, zh, distr


with open(f"{sys.argv[1]}_apr_sg.pkl", "wb") as file:
    pickle.dump(sg_apriori(f"{sys.argv[1]}.pkl"), file)
# envel = envelope(zh, ALT_RANGE[0], ALT_RANGE[1], SIGMA)
# plt.figure()
# plt.plot(distr, zh, color="blue")
# plt.plot(smooth, zh, color="red")
# plt.plot(envel * 5e13, zh, color="green")
# plt.plot(smooth * envel, zh, color="cyan")
# plt.show()
# breakpoint()
