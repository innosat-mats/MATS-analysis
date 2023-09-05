#%%calculate jacobian for all measurements
import numpy as np
import pandas as pd
from mats_l2_processing.forward_model import calc_jacobian, cart2sph
import pickle
from datetime import datetime
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l2_processing.inverse_model import do_inversion
from mats_l2_processing.grids import center_grid,localgrid_to_lat_lon_alt_3D

# #%%
# starttime=datetime(2023,3,31,21,0)
# stoptime=datetime(2023,3,31,22,35)
# dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.4")

# df = dftop[dftop['channel'] == 'IR2'].dropna().reset_index(drop=True)

# offsets = np.arange(0,len(df)-100,50)
# #select part of orbit
# for offset in offsets:
#     num_profiles = 100 #use 50 profiles for inversion
#     df_batch = df.loc[offset:offset+num_profiles-1]
#     df_batch = df_batch.reset_index(drop=True)
#     columns = np.arange(0,df["NCOL"][0],2)
#     rows = np.arange(0,df["NROW"][0]-10,1)

#     y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local = calc_jacobian(df_batch,columns,rows)

#     y = y.reshape(-1)
#     y = np.matrix(y)

#     x_hat = do_inversion(ks,y)

#     filename = str(offset) + ".pkl"
#     with open(filename, "wb") as file:
#         pickle.dump((x_hat, y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local), file)


# %% reshape x_hat

filename = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/0.pkl"

with open(filename, "rb") as file:
    result = pickle.load(file)

[x_hat, y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local] = result

#%%


altitude_grid = center_grid(altitude_grid_edges)
alongtrack_grid = center_grid(alongtrack_grid_edges)
acrosstrack_grid = center_grid(acrosstrack_grid_edges)
x_hat_reshaped = np.array(x_hat).reshape(len(altitude_grid),len(acrosstrack_grid),len(alongtrack_grid))

[alt,lat,lon,r] = localgrid_to_lat_lon_alt_3D(altitude_grid,alongtrack_grid,acrosstrack_grid,ecef_to_local)
lat = np.rad2deg(lat)
lon = np.rad2deg(lon)

## Plot crossections

# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=0))
# plt.show()
# # %%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=1))
# plt.show()

# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=2))
# plt.show()
# # %%
# %%
import numpy as np
from scipy.interpolate import griddata

# Define your 3D data points and corresponding values
data_points = np.array([alt.flatten(),lon.flatten(),lat.flatten()]).T  # Shape: (N, 3)
values = x_hat       # Shape: (N,)

# Define the regular grid
alt_grid = np.arange(alt.min(), alt.max(), 1e3)
lon_grid = np.arange(lon.min(), lon.max(), 0.5)
lat_grid = np.arange(lat.min(), lat.max(), 0.5)
x_grid, y_grid, z_grid = np.meshgrid(alt_grid, lon_grid, lat_grid, indexing='ij')

# Interpolate using linear interpolation
interpolated_values = griddata(data_points, values, (x_grid, y_grid, z_grid), method='linear')
# %%
