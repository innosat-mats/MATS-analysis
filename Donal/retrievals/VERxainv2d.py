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
from oem_functions import linear_oem_sp
import pickle
import scipy.sparse as sp


ir2=xr.load_dataset('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/Donal/retrievals/IR2mars31vertest_2.nc')
filename = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/Donal/retrievals/jacobian_2.pkl"
with open(filename, "rb") as file:
    [edges, k] = pickle.load(file)


# %%
ch=ir2
mr = []
error2_retrieval = []
error2_smoothing = []
ver = []
limb_fit = []
time_save = []
xa=np.ones([k.shape[1]])
Sa=sp.csr_matrix(np.diag(xa) * ch.profile.max().values)
y=ch.profile.values.reshape(-1)
Se=sp.csr_matrix(np.diag(np.ones([k.shape[0]])) * ch.profile.max().values)*1e10
xa=0*xa
x_hat, G, A, Ss, Sm = linear_oem_sp(k, Se, Sa, y, xa)
ver.append(x_hat)
mr.append(A.sum(axis=1)) #sum over rows 
error2_retrieval.append(np.diag(Sm.toarray()))
error2_smoothing.append(np.diag(Ss.toarray()))
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
plt.plot(k*x_hat)
plt.plot(y)
plt.show()
# %%
x_hat_reshape1 = x_hat.reshape(50,-1)
plt.plot(x_hat_reshape1)
plt.show()

# %%
plt.pcolor(edges[2][0:-1],(edges[0][0:-1]-6375000)*1e-3,x_hat_reshape1)
plt.xlim([-0.1,0.11])
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([70,120])
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()


# %%
plt.contourf(edges[2][0:-1],(edges[0][0:-1]-6375000)*1e-3,x_hat_reshape1,100)
plt.xlim([-0.1,0.11])
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([70,120])
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()

# %%
