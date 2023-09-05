#%%
import numpy as np
import numba

import fast_histogram

#%%
bins = (100, 100)
ranges = ((-1, 1), (-1, 1))
bins = np.asarray(bins).astype(np.int64)
ranges = np.asarray(ranges).astype(np.float64)

edges = (
    np.linspace(*ranges[0, :], bins[0] + 1),
    np.linspace(*ranges[1, :], bins[1] + 1),
)
#%%
np.random.seed(42)
vals = np.random.normal(size=[2, 1000000]).astype(np.float32)
vals1d = np.random.normal(size=[10000000]).astype(np.float32)
# %%
H = fast_histogram.histogramdd(vals.T, bins=100, range=((-1, 1), (-1, 1)))
# %%
