#%%
import numpy as np
import pandas as pd
from mats_l2_processing.forward_model import calc_jacobian
import pickle
from datetime import datetime
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l2_processing.inverse_model import do_inversion
from mats_l2_processing.grids import center_grid,localgrid_to_lat_lon_alt_3D, sph2cart, cart2sph
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import xarray as xr
folder = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/L2_test3/"
filenames = ["0","50","100","150","200","250"]



with open(folder + filenames[2] + ".pkl", "rb") as file:
    result = pickle.load(file)
[x_hat, y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local,non_uniform_ecef_grid_altitude,non_uniform_ecef_grid_lon,non_uniform_ecef_grid_lat,non_uniform_ecef_grid_r] = result

#%%

images = y.reshape(200,177,43)
plt.plot(images[0,:])
plt.show()
# %%
k_r = ks.reshape(200,177*43*178*44*101)

#%%
k_r = k_r.reshape(177,43,178,44,101)
k_r = k_r.reshape(200,177,43,178,44,101)
k_r = k_r.reshape(200,177,43,178,44,101)
k_r = k_r.reshape(200,177,43,178,44,101)
k_r = k_r.reshape(200,177,43,178,44,101)
# %%
