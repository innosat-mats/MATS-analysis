#%%
from mat4py import loadmat
import proplot as pplt
import pandas as pd
import numpy as np
import xarray as xr
from mats_utils.rawdata.read_data import read_MATS_data
from astropy.time import Time
from datetime import datetime, timezone
import datetime as DT
import calendar
from geopy.distance import geodesic
import cartopy.crs as ccrs
from time import strftime
from mats_utils.geolocation.coordinates import col_heights, satpos
import os
import pysolar
import matplotlib.pyplot as plt
import ref_index as RF

#%%
"""- lambda: OSIRIS' wavelength pixels

- L: example OSIRIS spectrum with A-Band dayglow

- L1-L4: the mean radiance that I calculate from 
that spectrum for the four MATS filter curves

- wIR1-wIR4: The "weights" that I apply when averaging 
the OSIRIS spectrum over the four MATS filter curves"""


data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/example_spectra/example_spectrum.mat')
lambdas=np.ndarray.flatten(np.array(data['lambda']))
Ls=np.ndarray.flatten(np.array(data['L']))
# %%
plt.plot(lambdas,Ls)
plt.xlim([310,875])
plt.ylim([0,8*10**14])
# %%

IR1=np.loadtxt('F-N-IR1-ABandCenter_transmission_air_6degr.dat',skiprows=1,unpack=True)
IR2=np.loadtxt('F-N-IR2-ABandTotal_air_6degr.dat',skiprows=1,unpack=True)
IR3=np.loadtxt('F-N-IR3-BgShort_transmission_air_6degr.dat',skiprows=1,unpack=True)
IR4=np.loadtxt('F-N-IR4-BgLong_transmission_air_6degr.dat',skiprows=1,unpack=True)


IR1[0,:]=RF.air2vac(IR1[0,:])
IR2[0,:]=RF.air2vac(IR2[0,:])
IR3[0,:]=RF.air2vac(IR3[0,:])
IR4[0,:]=RF.air2vac(IR4[0,:])

IR1[1,:]/=100
IR2[1,:]/=100
IR3[1,:]/=100
IR4[1,:]/=100

grid = np.arange(750, 780, 0.01)

filter1=np.interp(grid, IR1[0,:],IR1[1,:],left=0,right=0)
filter2=np.interp(grid, IR2[0,:],IR2[1,:],left=0,right=0)
filter3=np.interp(grid, IR3[0,:],IR3[1,:],left=0,right=0)
filter4=np.interp(grid, IR4[0,:],IR4[1,:],left=0,right=0)

Lsi=np.interp(grid,lambdas,Ls,left=0,right=0)

IR1_m = np.sum(filter1*Lsi)/filter1.sum()
IR2_m = np.sum(filter2*Lsi)/filter2.sum()
IR3_m = np.sum(filter3*Lsi)/filter3.sum()
IR4_m = np.sum(filter4*Lsi)/filter4.sum()

# %%

# plot filter1*Lsi
#plt.plot(grid,filter1*Lsi)
plt.plot(grid,filter2*Lsi)
#plt.plot(grid,filter3*Lsi)
#plt.plot(grid,filter4*Lsi)

# %%
