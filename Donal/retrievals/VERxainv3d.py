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

def cart2sph(pos):
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / radius)

    return np.array([radius,longitude,latitude]).T

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
x_hat = np.load('xhat_3d.npy')
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

rs = ((edges[0][0:-1]+edges[0][1:])/2)
lons = (edges[1][0:-1]+edges[1][1:])/2
lats = (edges[2][0:-1]+edges[2][1:])/2

#%%
filename = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/Donal/retrievals/ecef_to_local.pkl"
with open(filename, "rb") as file:
    ecef_to_local = pickle.load(file)

#%%
ret_ecef=ecef_to_local.inv().apply(np.array([rs[0]*np.cos(lats),np.zeros(len(lats)),rs[0]*np.sin(lats)]).T)# %%
ret_lats=np.rad2deg(cart2sph(ret_ecef)[:,2])
#%%
def geoid_radius(latitude):
    '''
    Function from GEOS5 class.
    GEOID_RADIUS calculates the radius of the geoid at the given latitude
    [Re] = geoid_radius(latitude) calculates the radius of the geoid (km)
    at the given latitude (degrees).
    ----------------------------------------------------------------
            Craig Haley 11-06-04
    ---------------------------------------------------------------
    '''
    DEGREE = np.pi / 180.0
    EQRAD = 6378.14 * 1000
    FLAT = 1.0 / 298.257
    Rmax = EQRAD
    Rmin = Rmax * (1.0 - FLAT)
    Re = np.sqrt(1./(np.cos(latitude*DEGREE)**2/Rmax**2
                + np.sin(latitude*DEGREE)**2/Rmin**2)) 
    return Re

ret_rads = np.zeros([len(rs),len(ret_lats)])
for i  in range(len(ret_lats)):
    lat=ret_lats[i]
    ret_rads[:,i]=(rs-geoid_radius(lat))*1e-3



# %%
plt.pcolor(ret_lats[:-5],ret_rads[:,:-5],np.sum(x_hat_reshape1[:,:,:-5],axis=1))
plt.xlabel('Latitudes (degrees)')
plt.ylabel('Altitude (km)')
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()

# %%
plt.pcolor(lons,ret_lats[:-5],np.sum(x_hat_reshape1[:,:,:-5],axis=0).T)
plt.ylabel('Latitudes (radians)')
plt.xlabel('Across track (radians)')
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()

# %%
plt.pcolor(lons,ret_rads[:,int(ret_rads.shape[1]/2)],np.sum(x_hat_reshape1[:,:,:-5],axis=2))
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()
#%%
x_hat_reshape1_pos = x_hat_reshape1
x_hat_reshape1_pos[x_hat_reshape1_pos<0] = 0

#%%
import plotly.graph_objects as go

X, Y, Z = np.meshgrid(rs, lons, lats[:-5],indexing='ij')
X_new = np.tile(ret_rads[:, np.newaxis, :-5], (1, 19, 1))
x_hat_reshape1_new = x_hat_reshape1[:,:,:-5]

x_hat_reshape1_new_interp = np.zeros(x_hat_reshape1_new.shape)
for i in range(X_new.shape[1]):
    for j in range(X_new.shape[2]):
        x_hat_reshape1_new_interp[:,i,j] = np.interp((rs-geoid_radius(ret_lats[0]))*1e-3,X_new[:,i,j],x_hat_reshape1_new[:,i,j])

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
    value=x_hat_reshape1_new_interp.reshape(-1),
    opacity=0.5, # needs to be small to see through all surfaces
    surface_count=100, # needs to be a large number for good volume rendering
    colorscale=colorscale,
    ))


#Setting axis limits
# fig.update_layout(
#     scene=dict(
#         xaxis=dict(range=[-1, 0.5]),
#         yaxis=dict(range=[-0.2, 0.2]),
#         zaxis=dict(range=[50, 110])
#     )
# )
fig.show()


#%%
# Create a figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create a 3D scatter plot
x_values = Z.reshape(-1)-(rs-geoid_radius(ret_lats[0]))*1e-3
y_values = Y.reshape(-1)
z_values = X.reshape(-1)
values=x_hat_reshape1_new.reshape(-1)

# Define colormap and normalization
cmap = plt.get_cmap('viridis')

alpha_value = 0.02
scatter = ax.scatter(x_values, y_values, z_values, c=values, cmap='viridis', alpha=alpha_value)

# Customize plot (add labels, title, etc. if needed)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_title('3D Scatter Plot of 3D Data')

# Show the plot
plt.show()
# # %%
# import numpy as np
# from mayavi import mlab

# # Create a figure
# mlab.figure()

# # Create an isosurface with transparent regions for zero values
# source = mlab.pipeline.scalar_field(values)
# isosurface = mlab.pipeline.iso_surface(source, contours=[values.min(), 0, values.max()])
# isosurface.actor.property.opacity = 0.5  # Set transparency for the isosurface

# # Add a colorbar
# mlab.colorbar()

# # Display the plot
# mlab.show()


# %%
import numpy as np
import plotly.graph_objects as go

# Create a 3D scatter plot
x_values = Z.reshape(-1)
y_values = Y.reshape(-1)
z_values = X_new.reshape(-1)
values=x_hat_reshape1_new.reshape(-1)

# Create the plot
fig = go.Figure(data=go.Isosurface(
    x=x_values.flatten(),
    y=y_values.flatten(),
    z=z_values.flatten(),
    value=values.flatten(),
    colorscale='Viridis',
    isomin=1e12,
    isomax=2e12,
    opacity=0.5,  # Set transparency for the isosurface
))

# Show the plot
fig.show()
# %%
# Create the plot
fig = go.Figure(data=go.Scatter3d(
    x=x_values,
    y=y_values,
    z=z_values,
    mode='markers',
    marker=dict(
        size=3,
        color=values,
        opacity=0.05,
        colorscale='greens',
        showscale=True,
        cmin=0,
        cmax=1e12,
    ),
))

# Show the plot
fig.show()
# %%
