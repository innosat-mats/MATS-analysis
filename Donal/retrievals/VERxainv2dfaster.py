#%%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from scipy.interpolate import CubicSpline
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
import xarray as xr
from scipy.sparse.linalg import inv,spsolve
from oem_functions import linear_oem_sp
import pickle
import scipy.sparse as sp
from oem import oem_basic_sparse_2
from os.path import expanduser
import time
import scipy.stats as stats

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

def cart2sph(pos):
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / radius)

    return np.array([radius,longitude,latitude]).T


ch=xr.load_dataset(expanduser('~donal/projekt/SIW/MATS-analysis/Donal/retrievals/IR2mars31vertest_100-200.nc'))
filename = expanduser("~donal/projekt/SIW/MATS-analysis/Donal/retrievals/jacobian2_100-200.pkl")
with open(filename, "rb") as file:
    [edges, k, ecef_to_local] = pickle.load(file)


# %%
def remove_empty_columns(matrix):

    _,I = np.nonzero(np.sum(matrix,0))
    _,I_0 = np.where(np.sum(matrix,0)==0)
    cleaned_matrix = matrix[:,I]

    return cleaned_matrix, I_0,I


def reinsert_zeros(vector, indexes):
    for idx in indexes:
        vector.insert(idx, 0)
        
#ch=ir2
mr = []
error2_retrieval = []
error2_smoothing = []
ver = []
limb_fit = []
time_save = []
#k=np.matrix(k)
#k_reduced,empty_cols,filled_cols = remove_empty_columns(k)
k_reduced=k.tocsc()
y = ch.profile.values.reshape(-1)

rs=(edges[0][0:-1]+edges[0][1::])/2
lons=(edges[1][0:-1]+edges[1][1::])/2
lats=(edges[2][0:-1]+edges[2][1::])/2
ret_ecef=ecef_to_local.inv().apply(np.array([rs[0]*np.cos(lats),np.zeros(len(lats)),rs[0]*np.sin(lats)]).T)
ret_lats=np.rad2deg(cart2sph(ret_ecef)[:,2])
# %%
xa=np.zeros([len(rs),len(lats)]).astype('float32')
for i,lat in enumerate(ret_lats):
  xa[:,i]= stats.norm.pdf(rs,geoid_radius(lat)+88000,2500)
xa_norm=1e12/stats.norm.pdf(geoid_radius(lat)+88000,geoid_radius(lat)+88000,2500)  /2
#xa=np.ones([k_reduced.shape[1]])
xa=xa_norm*xa.reshape(-1)

Sa_inv=sp.diags(np.ones([xa.shape[0]]),0).astype('float32') * (1/np.max(y)) * 1e9
Se_inv=sp.diags(np.ones([k_reduced.shape[0]]),0).astype('float32') * (1/np.max(y))
# Sa=sp.diags(np.ones([xa.shape[0]]),0) * np.max(y) * 5e-10
# Se=sp.diags(np.ones([k.shape[0]]),0) * np.max(y) 
# print('start inverses')
# Se_inv = inv(Se)
# Sa_inv = inv(Sa)
# print ('Inverses complete')

#xa=0*xa


#xa = np.matrix(xa).Tk
#y = np.matrix(y).T


#%%
start_time = time.time()
x_hat = oem_basic_sparse_2(y, k_reduced, xa, Se_inv, Sa_inv, maxiter=1000)
end_time = time.time()
np.save('xhat2_snabb_100-200.npy',x_hat)
print(end_time-start_time,' sec')
#%%
#x_hat = np.load('xhat1_snabb_200-300.npy')
#%%
plt.figure()
plt.plot((k_reduced @ (x_hat)))
plt.plot(y)
plt.show()
#%%
plt.figure()
plt.plot((k_reduced @ (x_hat)))
plt.plot(y)
plt.plot(k_reduced @ xa)
plt.legend(['K x_hat','y','K xa'])
plt.xlim([15000,15500])
plt.show()

#x_hat_old = x_hat
#x_hat = np.zeros([k.shape[1],1])
#x_hat[filled_cols] = x_hat_old.T

#%%
x_hat_reshape1 = np.array(x_hat).reshape(len(edges[0])-1,-1 )
#x_hat_reshape1 = np.array(x_hat).reshape(-1,len(edges[0])-1 )

plt.plot(x_hat_reshape1[:,1:-2])
plt.show()
# %%
plt.pcolor((edges[2][0:-1]+edges[2][1::])/2,(edges[0][0:-1]-6360000)*1e-3,x_hat_reshape1)
#plt.xlim([-0.1,0.11])
plt.xlim([edges[2][5],edges[2][-7]])
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([50,120])
plt.clim()
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()


# %%
plt.contourf((edges[2][0:-1]+edges[2][1::])/2,(edges[0][0:-1]-6360000)*1e-3,x_hat_reshape1)
plt.xlim([edges[2][5],edges[2][-7]])
plt.colorbar()
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([50,120])
plt.title('2D tomographic retrieval nightglow IR1')
plt.show()

# %%


# %%
zmax=150.
zmin=50.
plot_data=np.zeros([len(ret_lats),int((zmax-zmin))])
#plot_data=np.zeros([len(ret_lats),len(rs)])

for i,lat in enumerate(ret_lats):

    xs=(rs-geoid_radius(lat))*1e-3
    print(i,lat,xs)
    print(x_hat_reshape1[:,i])
    print(np.interp(np.arange(zmin,zmax),xs,x_hat_reshape1[:,i]))
    plot_data[i,:] = np.interp(np.arange(zmin,zmax),xs,x_hat_reshape1[:,i])
    #plot_data[i,:] = (rs-geoid_radius(lat))*1e-3

plt.figure()
plt.contourf(lats[1:-1],np.arange(zmin,zmax),plot_data.T[:,1:-1])
#plt.clim([0,1e12])
plt.xlim([edges[2][1],edges[2][-2]])
plt.colorbar()
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Corrected Altitude (km)')
plt.ylim([50,120])

plt.title('2D tomographic retrieval nightglow IR2')
plt.show()
# %%
x_hat_ir1=np.load('xhat1_snabb_100-200.npy')/1e9*3.57
x_hat_ir2=np.load('xhat2_snabb_100-200.npy')/1e9*8.16
ratio=(x_hat_ir1/x_hat_ir2).reshape(len(edges[0])-1,-1 )

# %%
plt.figure()
plt.contourf((edges[2][0:-1]+edges[2][1::])/2,
             (edges[0][0:-1]-6360000)*1e-3,ratio,
             levels=np.linspace(0,2,11))
plt.clim([0,2])
plt.colorbar()
#plt.clim([0,1e12])
plt.xlim([edges[2][1],edges[2][-2]])

# %%
overlap=np.array([0.78584363, 0.7656091 , 0.74607243, 0.72726005, 0.70918753,
       0.69185982, 0.67527261, 0.659414  , 0.64426607, 0.62980651,
       0.61600987, 0.60284864, 0.59029418, 0.57831736, 0.56688911,
       0.55598082, 0.54556466, 0.53561371, 0.52610218, 0.51700545,
       0.50830014])
temps=np.array([100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220,
       230, 240, 250, 260, 270, 280, 290, 300])
tempspl=CubicSpline(overlap[-1::-1],temps[-1::-1])
# %%
plt.figure()


plt.pcolor((edges[2][0:-1]+edges[2][1::])/2,
             (edges[0][0:-1]-6360000)*1e-3,tempspl(0.6*ratio))#,levels=np.linspace(160,250,10) )
plt.clim([160,250])
plt.xlim([edges[2][1],edges[2][-2]])
plt.ylim([60,110])
plt.colorbar()
# %%
