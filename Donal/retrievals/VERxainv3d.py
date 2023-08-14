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
subsample = 5
k = k[::subsample,:]
y = profiles.reshape(-1)
y = y[::subsample]
xa=np.ones([k.shape[1]])
Sa=sp.diags(np.ones([xa.shape[0]]),0).astype('float32') * np.max(y)
Se=sp.diags(np.ones([k.shape[0]]),0).astype('float32') * np.max(y) *1e10
xa=0*xa

Se_inv = sp.linalg.inv(Se)
Sa_inv = sp.linalg.inv(Sa)


#%%
x_hat = oem.oem_cgs(k, Se, Sa, y, xa,maxiter=1000)
np.save('xhat.npy',x_hat)
ver.append(x_hat)
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
x_hat = np.load('xhat.npy')
# %%
plt.plot(k*x_hat)
plt.plot(y)
plt.show()
# %%
x_hat_reshape1 = x_hat.reshape(len(edges[0])-1
,len(edges[1])-1
,len(edges[2])-1)
plt.plot(x_hat_reshape1[:,3,:])
plt.show()

# %%
plt.pcolor(edges[2][0:-1],(edges[0][0:-1]-6375000)*1e-3,x_hat_reshape1[:,3,:])
plt.xlim([-0.1,0.11])
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([70,120])
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()


# %%
plt.contourf(edges[2][0:-1],(edges[0][0:-1]-6375000)*1e-3,x_hat_reshape1[:,3,:],100)
plt.xlim([-0.1,0.11])
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([70,120])
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()

#%%
import plotly.graph_objects as go

X, Y, Z = np.meshgrid(edges[0][:-1], edges[1][:-1], edges[2][:-1])#,indexing='ij')

fig = go.Figure(data=go.Volume(
    x = X.reshape(-1)-6375000,
    y = Y.reshape(-1),
    z = Z.reshape(-1),
    value=x_hat_reshape1.reshape(-1),
    opacity=0.5, # needs to be small to see through all surfaces
    surface_count=100, # needs to be a large number for good volume rendering
    ))
fig.show()

# %%
def volume_plot(volume_data):
    # Create a figure and a 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Get the dimensions of the volume data
    depth, height, width = volume_data.shape

    # Create a 3D meshgrid for the volume plot
    x, y, z = np.meshgrid(range(width), range(height), range(depth))

    # Flatten the arrays to create points for plotting
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    volume_data = volume_data.flatten()

    # Plot the volume data using scatter3D
    ax.scatter3D(x, y, z, c=volume_data, cmap='viridis', s=50)

    # Set labels for the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Add a color bar to show the color scale
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('Intensity')

    # Show the plot
    plt.show()

volume_plot(x_hat_reshape1)
# %%
