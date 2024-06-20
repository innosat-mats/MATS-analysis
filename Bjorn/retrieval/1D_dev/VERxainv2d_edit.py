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
#from numpy.linalg import inv
#from oem_functions import linear_oem_sp
import pickle
import scipy.sparse as sp
from scipy.sparse.linalg import inv,spsolve
from oem import oem_basic_sparse_2
from os.path import expanduser

def cart2sph(pos):
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / radius)

    return np.array([radius,longitude,latitude]).T


ir2=xr.load_dataset(expanduser('~donal/projekt/SIW/MATS-analysis/Donal/retrievals/IR1mars31vertest_100-200.nc'))
filename = expanduser("~donal/projekt/SIW/MATS-analysis/Donal/retrievals/jacobian_100-200.pkl")
with open(filename, "rb") as file:
    [edges, k, ecef_to_local] = pickle.load(file)

#%% oem for sparse matrix
def linear_oem_sp(K, Se, Sa, y, xa):
    K = K.tocsc()
    K.eliminate_zeros()
    Sa = Sa.tocsc()
    Sa.eliminate_zeros()
    Se = Se.tocsc()
    Se.eliminate_zeros()
    print ('made the matrices')
    if len(y)<len(xa): # m form
        print ('doing the len y < len xa')
        G = Sa.dot(K.T).dot(inv(K.dot(Sa).dot(K.T) + Se))

    else:
        print('start inverses')
        Se_inv = inv(Se)
        Sa_inv = inv(Sa)
        print ('Inverses complete')
        G = inv(K.T @ Se_inv @ K + Sa_inv) @ K.T @ Se_inv
#        G = spsolve(K.T @ Se_inv @ K + Sa_inv, K.T @ Se_inv)
    x_hat = xa + G @ (y - K @ xa)
    A = G.dot(K)
    I = sp.identity(len(xa))
    Ss = (A - I).dot(Sa).dot((A - I).T) # smoothing error
    Sm = G.dot(Se).dot(G.T) #retrieval noise 
    
    return x_hat, G, A, Ss, Sm
# %%
ch=ir2
mr = []
error2_retrieval = []
error2_smoothing = []
ver = []
limb_fit = []
time_save = []
xa=np.ones([k.shape[1]]).T
y=ch.profile.values.reshape(-1).T
Sa=sp.diags(np.ones([xa.shape[0]]),0) * np.max(y) * 5e-10
Se=sp.diags(np.ones([k.shape[0]]),0) * np.max(y) 
xa=0*xa
#%%
x_hat, G, A, Ss, Sm = linear_oem_sp(k, Se, Sa, y, xa)
ver.append(x_hat)
np.save('xhat1_G_300-400.npy',x_hat)
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
x_hat = np.load('/Users/donal/projekt/xhat1_G_100-200.npy')
x_hat_reshape1 = x_hat.reshape(len(edges[0])-1,-1)

#%%
plt.plot((k@x_hat))
plt.plot(y)
plt.show()
# %%
plt.plot(x_hat_reshape1)
plt.show()

# %%
plt.pcolor((edges[2][0:-1]+edges[2][1::])/2,(edges[0][0:-1]-6375000)*1e-3,x_hat_reshape1)
#plt.xlim([-0.1,0.11])
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([70,120])
plt.clim()
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()


# %%
plt.contourf((edges[2][0:-1]+edges[2][1::])/2,(edges[0][0:-1]-6375000)*1e-3,x_hat_reshape1)
#plt.xlim([-1,1])
plt.xlim([edges[2][1],edges[2][-2]])
plt.colorbar()
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Altitude (km)')
plt.ylim([50,120])
plt.title('2D tomographic retrieval nightglow IR1')
plt.show()

# %%
rs=(edges[0][0:-1]+edges[0][1::])/2
lons=(edges[1][0:-1]+edges[1][1::])/2
lats=(edges[2][0:-1]+edges[2][1::])/2
ret_ecef=ecef_to_local.inv().apply(np.array([rs[0]*np.cos(lats),np.zeros(len(lats)),rs[0]*np.sin(lats)]).T)# %%

# %%
ret_lats=np.rad2deg(cart2sph(ret_ecef)[:,2])
# %%
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
# %%
zmax=150.
zmin=50.
plot_data=np.zeros([len(ret_lats),int((zmax-zmin))])
#plot_data=np.zeros([len(ret_lats),len(rs)])

#for i,lat in enumerate(ret_lats):
for i  in range(len(ret_lats)):
    lat=ret_lats[i]
    xs=(rs-geoid_radius(lat))*1e-3
    print(i,lat,xs)
    print(x_hat_reshape1[:,i])
    print(np.interp(np.arange(zmin,zmax),xs,x_hat_reshape1[:,i]))
    plot_data[i,:] = np.interp(np.arange(zmin,zmax),xs,x_hat_reshape1[:,i])
    #plot_data[i,:] = (rs-geoid_radius(lat))*1e-3

plt.figure()
plt.contourf(lats,np.arange(zmin,zmax),plot_data.T)
plt.xlim([edges[2][1],edges[2][-2]])
plt.colorbar()
plt.xlabel('Angle Along Orbit (radians)')
plt.ylabel('Corrected Altitude (km)')
plt.ylim([50,120])
plt.title('2D tomographic retrieval nightglow IR2')
plt.show()
# %%
