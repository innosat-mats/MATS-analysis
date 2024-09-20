#%%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
from mats_l2_processing.inverse_model import Sa_inv_tikhonov
#import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline,interp1d
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
import xarray as xr
from numpy.linalg import inv
from fast_histogram import histogramdd
from bisect import bisect_left
import matplotlib.pyplot as plt
import sys
import scipy.sparse as sparse
from scipy.stats import norm


def linear_oem(K, Se, Sa_inv, y, xa):
    # Adapted from Donal (uses Sa_inv instead)
    if len(y.shape) == 1:
        y = y.reshape(len(y),1)
    if len(xa.shape) == 1:
        xa = xa.reshape(len(xa),1)
        
    #if len(y)<len(xa): # m form
    #    G = Sa.dot(K.T).dot(np.linalg.inv(K.dot(Sa).dot(K.T) + Se))
        
    #else: # n form
    Se_inv = np.linalg.inv(Se)
    Sa_inv = Sa_inv.toarray()
    #Sa_inv = np.linalg.inv(Sa)
    G = np.linalg.inv(K.T.dot(Se_inv).dot(K) + Sa_inv).dot(K.T).dot(Se_inv)
#        G= np.linalg.solve(K.T.dot(Se_inv).dot(K) + Sa_inv, (K.T).dot(Se_inv))
    
    x_hat = xa + G.dot(y - K.dot(xa)) 
    A = G.dot(K)
    #Ss = (A - np.identity(len(xa))).dot(Sa).dot((A - np.identity(len(xa))).T) # smoothing error
    Sm = G.dot(Se).dot(G.T) #retrieval noise 
    
    return x_hat.squeeze()


def xainvert(ch,retrival_heights, weight_0, weight_1, weight_2, xa=None):
    # TBD comment
    # invert using oem

    ver = []

    for n in range(len(ch.time)): 
        profile=ch.profile[n].values
        Se=np.diag(profile)
        k=ch.k[n].values

        #### xa
        if xa is None:
            xa=np.ones_like(retrival_heights)
            xa=0*xa

        Sa_inv,_=Sa_inv_tikhonov(np.array([retrival_heights]), weight_0, [weight_1], volume_factors=False, store_terms=False)
        #Sa=np.diag(xa) * np.max(profile) / (1.3*10**11)

        x0 = linear_oem(k, Se, Sa_inv, profile, xa)
        ver.append(x0)
        #mr.append(A.sum(axis=1)) #sum over rows 
        #error2_retrieval.append(np.diag(Sm))
        #error2_smoothing.append(np.diag(Ss))

    return ver


#%%

retrival_heights= np.arange(60,100,1)
retrival_heights= np.arange(60,105,1)
s_140=1700e3
steps=100 #m steps
s_steps = np.arange(s_140,s_140 + 2e6,steps)
abs_bool=True

channels=['IR1','IR2']
channels=['IR1']

for channel in  channels:
    files='/media/waves/AVAGO/data/MATS/1D_inversions/v0.1/mar_w1/MAR_first_week/*.nc'
    #files=f'/media/waves/AVAGO/data/MATS/1D_inversions/v0.1/IR_comparison/{channel}_Mar1_l1bv0.5_sepabs.nc'

    ch=xr.open_mfdataset(files)

    #ch=ch.where(ch.)
    # xa = 0
    weight_0 = 0
    weight_1 = 0
    # normal xa (works well)
    weight_0 = 1e-3
    weight_1 = 6e-3

    # normal xa (try to get midlats)
    #weight_0 = 0
    #weight_1 = 9*1e-1

    #weight_0=0
    #weight_1=0

    xa=1e11*(2e3*1e6/1e11+norm.pdf(retrival_heights,83,4))

    xa=None
    ver = xainvert(ch,retrival_heights, weight_0, weight_1,xa)


    result_1d = xr.Dataset().update({
            'time': (['time'], ch.time.values),
            'z_r': (['z_r',], ch.z_r.values, {'units': 'km'}),
            'ver': (['time','z_r'], ver, {'long name': 'VER', 'units': 'photons cm-3 s-1'}),
            'latitude': (['time',], ch.TPlat.values),
            'longitude': (['time',], ch.TPlon.values),
            'channel': (['time',], ch.channel.values),
            })

    result_1d.ver.values = result_1d.ver.values * 4*np.pi / 10**6

    result_1d.to_netcdf(f"{channel}_MarW1_l1bv0.5_VER.nc")


# %%
