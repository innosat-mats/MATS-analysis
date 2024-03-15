# %%
import numpy as np
import pandas as pd
from numba import njit
from bisect import bisect_left
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos, findheight
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline, interp1d
import scipy.sparse as sp
from skyfield import api as sfapi
from skyfield.framelib import itrs
from os.path import expanduser
import xarray as xr
import pickle
import plotly.graph_objects as go

# from fast_histogram import histogramdd
import scipy.stats as stats
from jax import jit, grad, jacfwd, jacrev, value_and_grad
import jax.numpy as jnp


# %%
def geoid_radius(latitude):
    """
    Function from GEOS5 class.
    GEOID_RADIUS calculates the radius of the geoid at the given latitude
    [Re] = geoid_radius(latitude) calculates the radius of the geoid (km)
    at the given latitude (degrees).
    ----------------------------------------------------------------
            Craig Haley 11-06-04
    ---------------------------------------------------------------
    """
    DEGREE = np.pi / 180.0
    EQRAD = 6378.14 * 1000
    FLAT = 1.0 / 298.257
    Rmax = EQRAD
    Rmin = Rmax * (1.0 - FLAT)
    Re = np.sqrt(
        1.0
        / (
            np.cos(latitude * DEGREE) ** 2 / Rmax**2
            + np.sin(latitude * DEGREE) ** 2 / Rmin**2
        )
    )
    return Re


@njit(cache=True)
def cart2sph(pos=np.array([[]])):
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    # radius =np.linalg.norm(pos,axis=1)
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / radius)
    return radius, longitude, latitude


ch = xr.load_dataset(
    expanduser(
        "~donal/projekt/SIW/MATS-analysis/Donal/retrievals/IR1IR2test_400-520.nc"
    )
)
filename = expanduser(
    "~donal/projekt/SIW/MATS-analysis/Donal/retrievals/jacobianIR1IR2_400-520.pkl"
)
with open(filename, "rb") as file:
    [edges, k, ecef_to_local] = pickle.load(file)
msis = xr.load_dataset(
    expanduser(
        "~donal/projekt/SIW/MATS-analysis/Donal/retrievals/Datafiles/msis_cmam_climatology_200.nc"
    )
)
with open("runningfile_2", "rb") as file:
    [i, irow, ir1calcs, ir2calcs, ir1grads, ir2grads, profiles, ks] = pickle.load(file)
with open("apriori", "rb") as file:
    [VERarray, Tarray, o2array] = pickle.load(file)
# Tempory corrections due to bad apriori values     
VERarray *= 4*np.pi * 4*np.pi /1e6  #  divided by 4pi instead of mutiplying and was in m-3 
correction= 4*np.pi * 4*np.pi /1e6 * 1e4 # as above but also to change from cm-2 to m -2 
gradcorr= 1e4
ir1calcs = [intens * correction for intens in ir1calcs]  
ir2calcs = [intens * correction for intens in ir2calcs]  
# 

rs = edges[0]  # (edges[0][0:-1]+edges[0][1::])/2
lons = edges[1]  # (edges[1][0:-1]+edges[1][1::])/2
lats = edges[2]  # (edges[2][0:-1]+edges[2][1::])/2
ret_ecef = ecef_to_local.inv().apply(
    np.array([rs[0] * np.cos(lats), np.zeros(len(lats)), rs[0] * np.sin(lats)]).T
)
ret_lats = np.rad2deg(np.array(cart2sph(ret_ecef)).T[:, 2])
ret_lons = np.rad2deg(np.array(cart2sph(ret_ecef)).T[:, 1])
# %%
plt.figure()
plt.pcolor(ret_lats, rs, Tarray[:, 0, :])
plt.colorbar()
plt.title("Temperature -Apriori")
plt.figure()
plt.pcolor(ret_lats, rs, VERarray[:, 0, :])
plt.title("VER -Apriori")
plt.colorbar()
# %%
plt.figure()
plt.pcolor(np.array(ir1calcs).reshape(121, -1).T)
plt.title("IR1calc")
plt.colorbar()
# %%
ch["ir1profile"].plot(y="z")

# %%
plt.figure()
plt.pcolor(ret_lats, rs / 1000 - 6375, ir1grads[0][0][:, 0, :])
# %%
# lets test for IR1 retreival
ks = sp.lil_array((len(ir1calcs), len(ret_lats) * len(rs)))
for i in range(len(ir1calcs)):
    ks[i, :] =  gradcorr* ir1grads[i][0].flatten()  # correction for m-2
newir1 = ks @ VERarray.flatten()
plt.figure()
plt.pcolor(newir1.reshape(len(ret_lats), -1).T)
plt.title('Kintensity @ VERapriori')
plt.colorbar()
# %%
y = np.hstack([ir1calcs, ir2calcs])
y0 = np.hstack([ch.ir1profile.values.flatten(), ch.ir2profile.values.flatten()])
x0 = np.hstack([VERarray.flatten(), Tarray.flatten()])
x = x0.copy()  # First guess
ks = sp.lil_array((len(y), 2 * len(ret_lats) * len(rs)))
for i in range(len(ir1calcs)):
    ks[i, :] = np.hstack([ gradcorr * ir1grads[i][0].flatten(), gradcorr*ir1grads[i][1].flatten()])
    ks[i + len(ir1calcs), :] = np.hstack(
        [gradcorr*ir2grads[i][0].flatten(), gradcorr*ir2grads[i][1].flatten()]
    )
# %%
ks = ks.tocsc()
Sa_inv = sp.diags(
    np.hstack(
        [
            np.ones([int(x.shape[0] / 2)]) * (1 / np.max(y)),
            np.ones([int(x.shape[0] / 2)]) / 900,
        ]
    )
)
Se_inv = sp.diags(np.ones([ks.shape[0]]), 0).astype("float32") * (1 / np.max(y))
gamma = 0
# %%
# ktSei=ks.T @ Se_inv
# S=sp.linalg.inv((1+gamma)*Sa_inv + ktSei @ ks)
# xr=S @ (ktSei@(y-y0) - Sa_inv@(x-x0))
# %%
# hopefully faster
ktSei = ks.T @ Se_inv
S = (1 + gamma) * Sa_inv + ktSei @ ks
xr = sp.linalg.spsolve(S, (ktSei @ (y - y0) - Sa_inv @ (x - x0)))
xnew = x + xr
# %%
plt.figure()
plt.pcolor(xnew[:13431].reshape(111, 121), vmin=0, vmax=1e6)
plt.title('Retreived VER one iteration')
plt.colorbar()
plt.figure()
plt.pcolor(xnew[13431:].reshape(111, 121), vmin=120, vmax=400)
plt.title('Retreived T one iteration')
plt.colorbar()

# %%
