#%%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
import xarray as xr
from numpy.linalg import inv
import oem as oem
import pickle
import scipy.sparse as sp
import time


def remove_empty_columns(matrix):

    _,I = np.nonzero(np.sum(matrix,0))
    _,I_0 = np.where(np.sum(matrix,0)==0)
    cleaned_matrix = matrix[:,I]

    return cleaned_matrix, I_0,I


def reinsert_zeros(vector, indexes):
    for idx in indexes:
        vector.insert(idx, 0)

filename = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/Donal/retrievals/jacobian_3.pkl"
with open(filename, "rb") as file:
    [profiles, edges, k] = pickle.load(file)


# %%
mr = []
error2_retrieval = []
error2_smoothing = []
ver = []
limb_fit = []
time_save = []
subsample = 1
k = k[::subsample,:]
k_reduced = k
k_reduced,empty_cols,filled_cols = remove_empty_columns(k)

y = profiles.reshape(-1)
y = y[::subsample]
xa=np.ones([k_reduced.shape[1]])
# Sa=sp.diags(np.ones([xa.shape[0]]),0).astype('float32') * (np.max(y)) * 10e-10
# Se=sp.diags(np.ones([k_reduced.shape[0]]),0).astype('float32') * (np.max(y))

Sa_inv=sp.diags(np.ones([xa.shape[0]]),0).astype('float32') * (1/np.max(y)) * 1e8
Se_inv=sp.diags(np.ones([k_reduced.shape[0]]),0).astype('float32') * (1/np.max(y))
xa=0*xa
xa = np.matrix(xa).T
y = np.matrix(y).T
# Se_inv_2 = sp.linalg.inv(Se)
# Sa_inv_2 = sp.linalg.inv(Sa)


#%%
start_time = time.time()
x_hat = oem.oem_basic_sparse_2(y, k_reduced, xa, Se_inv, Sa_inv, maxiter=1000)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")

# xa = reinsert_zeros(xa, empty_cols)
# x_hat = reinsert_zeros(x_hat, empty_cols)

np.save('xhat_3d.npy',x_hat)

#mr.append(A.sum(axis=1)) #sum over rows 
#error2_retrieval.append(np.diag(Sm.toarray()))
#error2_smoothing.append(np.diag(Ss.toarray()))
# %%
# import plotly.graph_objects as go
# X, Y, Z = np.meshgrid(edges[0][:-1], edges[1][:-1], edges[2][:-1],indexing='ij')

# fig = go.Figure(data=go.Volume(
#     x=X.flatten()-6421665,
#     y=Y.flatten(),
#     z=Z.flatten(),
#     value=x_hat.flatten(),
#     isomin=1e13,
#     opacity=0.2, # needs to be small to see through all surfaces
#     surface_count=20, # needs to be a large number for good volume rendering
#     ))
# fig.show()
# %%
#x_hat = np.load('xhat_3d.npy')
# %%
plt.plot(y)
plt.plot(k_reduced.dot(x_hat),':')
plt.show()
# %%
x_hat_old = x_hat
x_hat = np.zeros([k.shape[1],1])
x_hat[filled_cols] = x_hat_old

# %%
x_hat_reshape1 = np.array(x_hat).reshape(len(edges[0])-1
,len(edges[1])-1
,len(edges[2])-1)
plt.plot(x_hat_reshape1[:,3,:])
plt.show()

#%%

z_grid = ((edges[0][0:-1]+edges[0][1:])/2-6375000)*1e-3
across_grid = (edges[1][0:-1]+edges[1][1:])/2
along_grid = (edges[2][0:-1]+edges[2][1:])/2

# %%
plt.pcolor(along_grid,z_grid,np.sum(x_hat_reshape1[:,:,:],axis=1))
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([70,120])
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()

# %%
plt.pcolor(across_grid,along_grid,np.sum(x_hat_reshape1[:,:,:],axis=0).T)
plt.ylabel('Angle Along Orbit (radians)')
plt.xlabel('Across track (radians)')
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()

# %%
plt.pcolor(across_grid,z_grid,np.sum(x_hat_reshape1[:,:,:],axis=2))
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()
#%%
x_hat_reshape1_pos = x_hat_reshape1
x_hat_reshape1_pos[x_hat_reshape1_pos<0] = 0

#%%
import plotly.graph_objects as go

X, Y, Z = np.meshgrid(z_grid, across_grid, along_grid,indexing='ij')

# Create a custom colorscale with white for 0 values
colorscale = [
    [0, 'white'],
    [0.2, 'rgb(255, 0, 0)'],  # Example: Red at 0.2
    [0.4, 'rgb(0, 255, 0)'],  # Example: Green at 0.4
    [1, 'rgb(0, 0, 255)']     # Example: Blue at 1
]

fig = go.Figure(data=go.Volume(
    x = Z.reshape(-1),
    y = Y.reshape(-1),
    z = X.reshape(-1),
    value=x_hat_reshape1.reshape(-1),
    opacity=0.5, # needs to be small to see through all surfaces
    surface_count=10, # needs to be a large number for good volume rendering
    colorscale=colorscale,
    ))


# Setting axis limits
fig.update_layout(
    scene=dict(
        xaxis=dict(range=[-1, 0.5]),
        yaxis=dict(range=[-0.2, 0.2]),
        zaxis=dict(range=[50, 110])
    )
)
fig.show()

# # %%
# def volume_plot(volume_data):
#     # Create a figure and a 3D axis
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')

#     # Get the dimensions of the volume data
#     depth, height, width = volume_data.shape

#     # Create a 3D meshgrid for the volume plot
#     x, y, z = np.meshgrid(range(width), range(height), range(depth),indexing='ij')

#     # Flatten the arrays to create points for plotting
#     #x = x.flatten()
#     # Flatten the arrays to c()

#     # Plot the volume data using scatter3D
#     ax.scatter3D(x, y, z, c=volume_data, cmap='viridis', s=50)

#     # Set labels for the axes
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')

#     # Show the plot
#     plt.show()

# volume_plot(x_hat_reshape1)
# # %%

# %%
