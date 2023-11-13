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


radius_grid = center_grid(altitude_grid_edges)
alongtrack_grid = center_grid(alongtrack_grid_edges)
acrosstrack_grid = center_grid(acrosstrack_grid_edges)


# #%%
# plt.plot(ks.dot(x_hat)[::10],label = 'fitted')
# plt.plot(y[::10].T,'r--',label = 'measured')
# plt.legend()
# plt.title('y space metrics')
# plt.show()

# #%%
# plt.plot(ks.dot(x_hat)[::10],label = 'fitted')
# plt.plot(y[::10].T,'r',label = 'measured')
# plt.legend()
# plt.title('y space metrics zoom')
# plt.xlim([5e4,5.01e4])
# #plt.ylim([-1e15,2.5e16])
# plt.show()

# #%%
# plt.plot(ks.dot(x_hat)-y)
# plt.title('residuals')
# plt.show()

# #%%
# plt.plot(ks.dot(x_hat)-y)
# plt.title('residuals zoom')
# plt.xlim([5e5,5.05e5])
# plt.ylim([-1*1e15,1*1e15])
# plt.show()

#%%
x_hat_reshape1 = np.array(x_hat).reshape(len(radius_grid)
,len(acrosstrack_grid)
,len(alongtrack_grid))
# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=0))
# plt.clim([0,3e14])
# plt.show()
# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=1))
# plt.clim([0,2e14])
# plt.show()

# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=2))
# plt.clim([0,4e14])
# plt.show()
# %%
#crop edges

#x_hat_reshape1_cropped = x_hat_reshape1[10:-40,5:-5,20:-20]
x_hat_reshape1_cropped = x_hat_reshape1


# # %%
# plt.pcolor(x_hat_reshape1_cropped[:,21,:])
# plt.clim([0,3e12])
# plt.colorbar()
# plt.title('2D retrieval crossection')
# plt.show()


# # %%
# plt.pcolor(x_hat_reshape1_cropped[90,:,:])
# plt.clim([0,3e12])
# plt.colorbar()
# plt.title('2D retrieval crossection')
# plt.show()


#%%
alt_grid = np.arange(70e3, 120e3, 0.5e3)
along_grid = np.arange(-0.5, 0.5, 0.005)
across_grid = np.arange(-0.15, 0.15, 0.0025)

#%%
_, acrosstrack_grid_expanded,alongtrack_grid_expanded = np.meshgrid(radius_grid, acrosstrack_grid,alongtrack_grid, indexing='ij')

#%%
data_points = np.array([non_uniform_ecef_grid_altitude.flatten(),acrosstrack_grid_expanded.flatten(),alongtrack_grid_expanded.flatten()]).T  # Shape: (N, 3)
#data_points = np.array([alt.flatten(),lon.flatten(),lat.flatten()]).T  # Shape: (N, 3)
values = x_hat_reshape1.flatten()       # Shape: (N,)

#%%
x_grid, y_grid, z_grid = np.meshgrid(alt_grid, across_grid, along_grid, indexing='ij')

#%%
interpgrid = np.array([x_grid.flatten(),y_grid.flatten(),z_grid.flatten()]).T
# Interpolate using linear interpolation
#%%
interpolated_values = griddata(data_points, values, interpgrid, method='nearest',rescale=True)
#%%
interpolated_values = interpolated_values.reshape(len(alt_grid),len(across_grid),len(along_grid))


# %%
plt.scatter(data_points[:,2], data_points[:,0], 2, values)
plt.clim([0,5e13])
plt.show()
# %%
plt.pcolor(along_grid*6371,alt_grid*1e-3,interpolated_values[:,60,:])
plt.clim([0,5e13])
plt.xlim([-2000,2000])
plt.xlabel('km along track')
plt.ylabel('Altitude (km)')
plt.colorbar()
plt.title('Data interpolated onto regular grid (mid pixel)')
plt.show()
# %%
fig, axs = plt.subplots(3)

hdl = axs[0].pcolor(along_grid*6371,across_grid*6371,interpolated_values[50,:,:])
hdl.set_clim(0,5e13)
axs[0].set_ylim([-250,250])
axs[0].set_xlabel('km along track')
axs[0].set_xlim([-2000,2000])
axs[0].set_ylabel('km across track')
axs[0].set_aspect('equal', adjustable='box')
axs[0].set_title('Data interpolated onto regular grid at 95 km')

hdl = axs[1].pcolor(along_grid*6371,across_grid*6371,interpolated_values[40,:,:])
hdl.set_clim([0,5e13])
axs[1].set_ylim([-250,250])
axs[1].set_xlabel('km along track')
axs[1].set_xlim([-2000,2000])
axs[1].set_ylabel('km across track')
axs[1].set_aspect('equal', adjustable='box')
axs[1].set_title('Data interpolated onto regular grid at 90 km')

hdl = axs[2].pcolor(along_grid*6371,across_grid*6371,interpolated_values[30,:,:])
hdl.set_clim([0,5e13])
axs[2].set_ylim([-250,250])
axs[2].set_xlabel('km along track')
axs[2].set_xlim([-2000,2000])
axs[2].set_ylabel('km across track')
axs[2].set_aspect('equal', adjustable='box')
axs[2].set_title('Data interpolated onto regular grid at 85 km')
# %%
