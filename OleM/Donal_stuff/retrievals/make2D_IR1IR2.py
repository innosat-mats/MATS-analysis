#%%
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
from scipy.interpolate import CubicSpline,interp1d
import scipy.sparse as sp
from skyfield import api as sfapi
from skyfield.framelib import itrs
from os.path import expanduser
import xarray as xr
import pickle
import plotly.graph_objects as go
#from fast_histogram import histogramdd
import scipy.stats as stats
from jax import jit,grad,jacfwd,jacrev,value_and_grad
import jax as jax

import jax.numpy as jnp
#import numpy as jnp

import ref_index as RF
import time
import timeit

def get_background_data():
    ## Read backgrounddata
    msis=xr.load_dataset(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/msis_cmam_climatology_200.nc')) #MSIS climatology
    sigma=np.load(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/o2Abandsigma100-600.npy')) #absoption cross-sections f(T,lambda)
    emission=np.load(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/o2Abandemission100-600.npy')) #normalized emission spectrum f(T,lambda)

    return msis,sigma,emission

def get_instrument_data():
    #Reads instrumentdata
    dftop = pd.read_pickle(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/verdec')) #VER data from MATS
    
    #Split data into channels
    ir1 = dftop[dftop['channel'] == 'IR1'].dropna().reset_index(drop=True)#[0:10]
    ir2 = dftop[dftop['channel'] == 'IR2'].dropna().reset_index(drop=True)#[0:10]
    ir3 = dftop[dftop['channel'] == 'IR3'].dropna().reset_index(drop=True)#[0:10]
    ir4 = dftop[dftop['channel'] == 'IR4'].dropna().reset_index(drop=True)#[0:10]
    uv1 = dftop[dftop['channel'] == 'UV1'].dropna().reset_index(drop=True)#[0:10]
    uv2 = dftop[dftop['channel'] == 'UV2'].dropna().reset_index(drop=True)#[0:10]

    #select part of orbit
    offset = 400
    num_profiles = 4 #use 50 profiles for inversion
    ir1 = ir1.loc[offset:offset+num_profiles].reset_index(drop=True)
    ir2 = ir2.loc[offset:offset+num_profiles].reset_index(drop=True)
    ir3 = ir3.loc[offset:offset+num_profiles].reset_index(drop=True)
    ir4 = ir4.loc[offset:offset+num_profiles].reset_index(drop=True)
    #offset=offset-414
    uv1 = uv1.loc[offset:offset+num_profiles].reset_index(drop=True)
    uv2 = uv2.loc[offset:offset+num_profiles].reset_index(drop=True)

    return ir1,ir2,ir3,ir4,uv1,uv2

def get_filtercurves():
    #%% Load filter curves

    IR1=np.loadtxt(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/F-N-IR1-ABandCenter_transmission_air_6degr.dat'),skiprows=1,unpack=True)
    IR2=np.loadtxt(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/F-N-IR2-ABandTotal_air_6degr.dat'),skiprows=1,unpack=True)
    IR3=np.loadtxt(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/F-N-IR3-BgShort_transmission_air_6degr.dat'),skiprows=1,unpack=True)
    IR4=np.loadtxt(expanduser('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/retrievals/Datafiles/F-N-IR4-BgLong_transmission_air_6degr.dat'),skiprows=1,unpack=True)

    #convert from wavelength in air to wavelength in vacuum for filtercurves
    grid = np.arange(12950, 13200, 0.002)
    IR1[0,:]=RF.air2vac(IR1[0,:])
    IR2[0,:]=RF.air2vac(IR2[0,:])
    IR1[1,:]/=100
    IR2[1,:]/=100
    IR3[1,:]/=100
    IR4[1,:]/=100
    filter1=np.interp(grid, 1e7/IR1[0,-1::-1],IR1[1,-1::-1],left=0,right=0)
    filter2=np.interp(grid, 1e7/IR2[0,-1::-1],IR2[1,-1::-1],left=0,right=0)
    filters=np.vstack([filter1,filter2])

    return filters
    
def prepare_profile(ch):
    # -*- coding: utf-8 -*-
    """Interpolates data in center column of input channel onto tangent height grid.
    Args:
        ch (MATS image): Data from a single Mats image as dataframe

    Returns:
        common_heights: common tangent height grid
        profile: profile on common tangent height grid
        heights: heights of each row of center column
    """

    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[0:-10, col-2:col+2].mean(axis=1)*1e13)
    profile = profile*1000/ch.TEXPMS #until files fixed
    common_heights = np.arange(60000,110250,1000)
    profile=np.interp(common_heights,heights[0:-10],profile)
    return common_heights, profile, heights


def prepare_measurment(ir1,ir2,ir3,ir4):
    # Subtracts background channels from IR1 and IR2
    # Converted to band strength
    common_heights,p1,_=prepare_profile(ir1)
    _,p2,_=prepare_profile(ir2)
    _,p3,_=prepare_profile(ir3)
    _,p4,_=prepare_profile(ir4)
    #plt.figure()
    #plt.semilogx(p1,z1,p2,z2,p3,z3,p4,z4)
    #If it is dark: Return altitudes and 
    if ir3.TPsza>98 : 
        raise ValueError('Nightglow not supported')
        #return common_heights,(p1-(p3+p4)/2)/(p2-(p3+p4)/2)*3.57/8.16
    
    p3=p3-p3[-4:].mean()/1.05
    p4=p4-p4[-4:].mean()/1.05
    p1=p1-(p3+p4)/2
    p2=p2-(p3+p4)/2
    IR1_BW = 3.57 #nm
    IR2_BW = 8.16 #nm
    ir1_photons = np.array(p1*IR1_BW)
    ir2_photons = np.array(p2*IR2_BW)
    return common_heights,ir1_photons,ir2_photons

def cart2sph(pos=np.array([[]])):
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    #radius =np.linalg.norm(pos,axis=1)
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / radius)
    return radius, longitude,latitude
    #return np.array([radius,longitude,latitude]).T

@jit
def cart2sph_jit(pos=jnp.array([[]])):
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    radius = jnp.sqrt(x**2 + y**2 + z**2)
    #radius =np.linalg.norm(pos,axis=1)
    longitude = jnp.arctan2(y, x)
    latitude = jnp.arcsin(z / radius)
    return radius, longitude,latitude
    #return np.array([radius,longitude,latitude]).T

def get_ecef_to_local(df,num_profiles):
    first = 0
    mid = int(num_profiles/2)
    last = num_profiles


    #Make local grid (select first, mid and last datapoint and generate grid along it)
    to_ecef = R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][first])))
    posecef_first = to_ecef.apply(df.afsTangentPointECI[first]).astype('float32')
    to_ecef = R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][mid])))
    posecef_mid = to_ecef.apply(df.afsTangentPointECI[mid]).astype('float32')
    to_ecef = R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][last])))
    posecef_last = to_ecef.apply(df.afsTangentPointECI[last]).astype('float32')
    observation_normal = np.cross(posecef_first,posecef_last)
    observation_normal = observation_normal/np.linalg.norm(observation_normal)
    posecef_mid_unit = posecef_mid/np.linalg.norm(posecef_mid)
    ecef_to_local = R.align_vectors([[1,0,0],[0,1,0]],[posecef_mid_unit,observation_normal])[0]

    return ecef_to_local

def get_ecivec(ch,quat):
    #Return a spline to compute ECI vector for a given row (ypixel)
    
    ypixels = np.linspace(0, ch['NROW']-1,5).astype('int')
    x, yv = pix_deg(ch, int(ch['NCOL']/2), ypixels)
    qp = R.from_quat(ch['qprime'])
    ecivec = np.zeros((3, len(yv)))
    for irow, y in enumerate(yv):
        los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
        ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
    cs_eci=CubicSpline(ypixels,ecivec.T) #cubic spline that give ecivec for rownumber

    return cs_eci,ecivec

#@njit(cache=True)
def get_ecivec_fast(quat,qprime,x,yv,ypixels):
    #Return a spline to compute ECI vector for a given row (ypixel)
    
    qp = qprime
    ecivec = np.zeros((3, len(yv)))
    for irow, y in enumerate(yv):
        los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
        ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
    cs_eci=CubicSpline(ypixels,ecivec.T) #cubic spline that give ecivec for rownumber

    return cs_eci,ecivec
    
@jit
def get_steps(sat_alt, theta, localR):
    # Estimate LOS distance to top of atmosphere (120 km)
    b = 2 * sat_alt * jnp.cos(theta) 
    root = jnp.sqrt(b**2 + 4 * ((120e3 + localR)**2 - sat_alt**2))
    s_120_1 = (-b - root) / 2  # LOS distance (nearest) to altitude of 120 km (spherical earth for local R)
    s_120_2 = (-b + root) / 2  # LOS distance (furthest) to altitude of 120 km (spherical earth for local R)

    return s_120_1,s_120_2

def get_steps_orig(sat_alt,theta,localR):

    #Estimate LOS distance to top of atmosphere (120 km)
    b=2*sat_alt*np.cos(theta) 
    root=np.sqrt(b**2+4*((120e3+localR)**2 - sat_alt**2))
    s_120_1 =(-b-root)/2 #LOS distance (nearest) to altitude of 120 km (spherical earth for local R)
    s_120_2 =(-b+root)/2 #LOS distance (furthest) to altitude of 120 km (spherical earth for local R)

    steps=2500 #m per step
    s_steps = np.arange(s_120_1,s_120_2,steps)

    return s_steps

def generate_grid(df):

    #Generate grid for mid measurement
    num_profiles = 4 #use 50 profiles for inversion
    mid = int(num_profiles/2)

    ecef_to_local = get_ecef_to_local(df,num_profiles) #Transform from ecef to local 

    ecipos= df['afsGnssStateJ2000'][mid][0:3] #position of satellite in ECI

    date = df['EXPDate'][mid]
    t = ts.from_datetime(date)
    localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat[mid], df.TPlon[mid], elevation_m=0).at(t).position.m) #get local R in mid of batch

    q = df['afsAttitudeState'][mid]
    quat = R.from_quat(np.roll(q, -1)) #quaterion of satellite pointing

    eci_to_ecef=R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][mid])))

    #calulate tangent geometries for limited number of rows and make a spline for it
    tanheights, _, _ = prepare_profile(df.iloc[mid])

    cs_eci_from_row,_ = get_ecivec(df.iloc[mid],quat)

    # #select row to get a LOS do define the local grid around
    irow = 90 #should prob be changed to (nrow/2 or a specific tanheight)
    ecivec=cs_eci_from_row(irow)
    ecefvec=eci_to_ecef.apply(ecivec)
    localvec=ecef_to_local.apply(ecefvec)


    sat_alt=np.linalg.norm(ecipos) #altitude of satellite (from earth center)
    theta=np.arccos(np.dot(ecipos,ecivec)/sat_alt) #angle between LOS and nadir

    steps = 2500  # m per step

    s_120_1,s_120_2 = get_steps(sat_alt,theta,localR)
    s_steps = jnp.arange(s_120_1,s_120_2,steps)   
    #hack to get around jax limitations
    #full_range = jnp.arange(1e6, 4e6, steps)
    #s_steps = full_range[jnp.logical_and(full_range >= s_120_1, full_range <= s_120_2)]
    
    #Get position in local coordinate system for all points along LOS
    satecefpos=eci_to_ecef.apply(ecipos)
    satlocalpos=ecef_to_local.apply(satecefpos)
    poslocal_i=(np.expand_dims(satlocalpos, axis=0).T+s_steps*np.expand_dims(localvec, axis=0).T).astype('float32')
    poslocal_i_sph = cart2sph(poslocal_i.T)   
    poslocal_i_sph=np.array(poslocal_i_sph).T   

    #Define retrieval model grid
    altitude_grid = np.arange(localR+30e3,localR+121e3,1e3)
    acrosstrack_grid = np.array([-0.3,0.3])
    #acrosstrack_grid = np.linspace(posecef_i_sph[:,1].min(),posecef_i_sph[:,1].max(),1)
    alongtrack_grid = np.linspace(poslocal_i_sph[:,2].min(),poslocal_i_sph[:,2].max(),100)
    alongtrack_grid[0] = alongtrack_grid[0]-0.5
    alongtrack_grid[-1] = alongtrack_grid[-1]+0.5

    #Define retrieval grid lat, lon, distance from earth center
    edges=[]
    edges.append(np.linspace(altitude_grid[0],altitude_grid[-1],110+1))
    edges.append(np.linspace(acrosstrack_grid[0],acrosstrack_grid[-1],1+1))
    edges.append(np.linspace(alongtrack_grid[0],alongtrack_grid[-1],120+1))
    radius=edges[0]#(edges[0][0:-1]+edges[0][1::])/2
    across_track=edges[1]#(edges[1][0:-1]+edges[1][1::])/2
    along_track=edges[2]#(edges[2][0:-1]+edges[2][1::])/2
    ret_ecef=ecef_to_local.inv().apply(np.array([radius[0]*np.cos(along_track),np.zeros(len(along_track)),radius[0]*np.sin(along_track)]).T)
    ret_lats=np.rad2deg(np.array(cart2sph(ret_ecef)).T[:,2])
    ret_lons=np.rad2deg(np.array(cart2sph(ret_ecef)).T[:,1])# %%

    return radius,across_track,along_track,ret_lats,ret_lons,date,edges,tanheights,ecef_to_local

def get_bg_atmosphere_on_grid(rs,lons,lats,ret_lats,ret_lons,d,msis):
    t = ts.from_datetime(d)
    Tarray=np.zeros([len(rs),len(lons)-1,len(lats)])
    VERarray=np.zeros_like(Tarray)
    o2array=np.zeros_like(Tarray)
    for i,retlat in enumerate(ret_lats):
        localR = np.linalg.norm(sfapi.wgs84.latlon(retlat, ret_lons[i], elevation_m=0).at(t).position.m)
        Tarray[:,0,i]=msis.T.sel(month=d.month).interp(lat=retlat,z=(rs-localR)/1000)
        o2array[:,0,i]=msis.o2.sel(month=d.month).interp(lat=retlat,z=(rs-localR)/1000)/1e6 # to cm-3
        VERarray[:,0,i]=2e3*1e6*stats.norm.pdf((rs-localR)/1000,88,4.5)+2e2*1e6 *np.exp(-((rs-localR)/1000-60)/20)
        
    Tarray = jnp.array(Tarray)
    o2array = jnp.array(o2array)
    VERarray = jnp.array(VERarray)

    return Tarray,o2array,VERarray

heights = []
ecipos = []
ecivecs = []
posecef_sph=[]
retrival_heights= np.arange(30,130,1)


#setup background atmosphere on retrieval grid




#%%

@jit
def interpT(x, xs, ys):
    ix = jnp.floor(x - 100).astype(jnp.int32)
    diff = xs[ix] - xs[ix - 1]
    return ((x - xs[ix - 1]) * ys[ix, :].T + (xs[ix] - x) * ys[ix - 1, :].T) / diff

@jit
def interppos_jit(radius, latitude, r_edges, along_edges):
    iz_values = jnp.floor((radius - r_edges[0]) / ((r_edges[1:] - r_edges[:-1]).sum() / (len(r_edges) - 1)))
    ix_values = jnp.floor((latitude - along_edges[0]) / ((along_edges[1:] - along_edges[:-1]).sum() / (len(along_edges) - 1)))

    # for i in I:
    #     iz_values = iz_values.at[i].set(0)
    # for j in J:
    #     iz_values = iz_values.at[j].set(len(r_edges))

    # I = jnp.where(ix_values<0)
    # J = jnp.where(ix_values>along_edges[-1])

    # for i in I:
    #     iz_values = iz_values.at[i].set(0)
    # for j in J:
    #     iz_values = iz_values.at[j].set(len(along_edges))

    iz = iz_values.astype(jnp.int32)
    ix = ix_values.astype(jnp.int32)

    return iz, ix
    #return((inArray[iz, 0, ix] + inArray[iz + 1, 0, ix + 1]) / 2)

def interppos(pos, edges):
    iz = jnp.array(jnp.floor((pos[:, 0] - edges[0][0]) / jnp.diff(edges[0]).mean()), dtype=int)
    ix = jnp.array(jnp.floor((pos[:, 2] - edges[2][0]) / jnp.diff(edges[2]).mean()), dtype=int)

    return iz,ix

##@jit
def cumsum_2d_axis_1(arr):
    out = jnp.empty_like(arr)
    out = out.at[:, 0].set(arr[:, 0])
    for i in range(1, arr.shape[1]):
        out = out.at[:, i].set(out[:, i - 1] + arr[:, i])
    return out

@jit
def ir1fun(radius, across,along,step,r_edges,across_edges,along_edges,filters,o2array,VER,Temps):
    #VER,Temps=atm # de-touple atmophere
    startT= jnp.linspace(100,600,501) #setting up temp grid
    iz,ix=interppos_jit(radius,along,r_edges,along_edges) #get temp at each point along path
    #pathtemps = jnp.empty_like(iz, dtype=jnp.float32)  # Assuming Temps dtype
    #for i in range(iz.shape[0]):
    #    pathtemps = pathtemps.at[i].set(Temps[iz[i], 0, ix[i]])
    pathtemps=Temps[iz,0,ix]
    pathtemps = jnp.array(pathtemps)
    sigmas=interpT(pathtemps,startT,sigma) #sigma for each point along the path
    emissions=interpT(pathtemps,startT,emission) #normalized emission spectrum for each point along the path
    
    o2=o2array[iz,0,ix]
    o2 = jnp.array(o2)
    # o2 = jnp.empty_like(iz, dtype=jnp.float32)  # Assuming Temps dtype
    # for i in range(iz.shape[0]):
    #     o2 = o2.at[i].set(o2array[iz[i], 0, ix[i]])

    tau = (sigmas*o2).cumsum(axis=1)*step*1e2
    #tau = cumsum_2d_axis_1(sigmas * o2) * s_steps * 1e2
    VERs=VER[iz,0,ix]
    VERs = jnp.array(VERs)
    #VERs = jnp.empty_like(iz, dtype=jnp.float32)  # Assuming Temps dtype
    #for i in range(iz.shape[0]):
    #    VERs = VERs.at[i].set(VER[iz[i], 0, ix[i]])

    VERs = VERs*step*1e2
    res=filters@(jnp.exp(-tau)*VERs*emissions) #get limb spectral radiance after filter at tanheight for VER (see other place)
    result = res[0].sum() #sum over wavelengths for filter 1

    return result

# ##@jit
# def compute_gradient_ir1(pos,path_step ,o2s,atm):
#     return grad(ir1fun,argnums=3)(pos,path_step ,o2s,atm) #take analytical (based on code) for argument #3 (which are VER and Temp (atm touple variable))

# ir1grad_non_compiled = compute_gradient_ir1
# #ir1grad_compiled = jit(compute_gradient_ir1)

##@jit
# def ir2fun(pos,path_step ,o2s,atm):
#        VER,Temps=atm
#     #    VER=jnp.array(VER)
#     #    Temps=jnp.array(Temps)
#     #    pos=jnp.array(pos)
#        startT= jnp.linspace(100,600,501)
#        pathtemps=interppos(pos,Temps)
#        #print(pathtemps)
#        sigmas=interpT(pathtemps,startT,sigma)
#        #print(sigmas.max())
#        emissions=interpT(pathtemps,startT,emission)
#        o2=interppos(pos,o2s)
#        tau = (sigmas*o2).cumsum(axis=1)*path_step * 1e2 #m -> cm
#        #print(tau)
#        VERs=interppos(pos,VER)*path_step*1e2 #m -> cm
#        res=filters@(jnp.exp(-tau)*VERs*emissions)
#        return res[1].sum()

#calculate jacobian for one pixel
#@jit
def calc_intensity(ecipos,ecivec,eci_to_ecef,ecef_to_local,r_edges,across_edges,along_edges,filters,localR,o2array,atm):
    VERarray,Tarray = atm

    ecefvec = jnp.matmul(eci_to_ecef,ecivec)
    localvec = jnp.matmul(ecef_to_local,ecefvec)

    satecefpos = jnp.matmul(eci_to_ecef,ecipos)
    satlocalpos = jnp.matmul(ecef_to_local,satecefpos)

    sat_alt=jnp.linalg.norm(ecipos) #altitude of satellite (from earth center)
    theta=jnp.arccos(jnp.dot(ecipos,ecivec)/sat_alt) #angle between LOS and nadir

    s_120_1,s_120_2 = get_steps(sat_alt,theta,localR)
    steps=2500 #m per step
    #hack to get around jax limitations
    #full_range = jnp.arange(1e6, 4e6, steps)
    #s_steps = full_range[jnp.logical_and(full_range >= s_120_1, full_range <= s_120_2)]
    s_steps = jnp.arange(s_120_1,s_120_2,steps)
    #s_steps = jnp.arange(1.5e6,4e6,steps)

    poslocal_i=jnp.expand_dims(satlocalpos, axis=0).T+s_steps*jnp.expand_dims(localvec, axis=0).T#.astype('float32')
    radius, across,along = cart2sph_jit(poslocal_i.T)   
    #poslocal_i_sph=jnp.array(np.asarray(np.array(poslocal_i_sph).T))
    #edges
    ir1calc=ir1fun(radius,across,along,steps,r_edges,across_edges,along_edges,filters,o2array,VERarray,Tarray)
    return ir1calc

#ir2grad_non_compiled = compute_gradient_ir2
#ir2grad_compiled = jit(compute_gradient_ir2)
#%%
## Run code ##

ts = sfapi.load.timescale()
ir1, ir2,_,_,_,_ = get_instrument_data()
filters = jnp.array(get_filtercurves())

radius,across_track,along_track,ret_lats,ret_lons,date,edges,tanheights,ecef_to_local = generate_grid(ir1)
msis, sigma,emission = get_background_data()
Tarray,o2array,VERarray = get_bg_atmosphere_on_grid(radius,across_track,along_track,ret_lats,ret_lons,date,msis)

Tarray = jnp.array(Tarray)
o2array = jnp.array(o2array)
VERarray = jnp.array(VERarray)
#for one tangenheigh

i = 2 #choose image
t = ts.from_datetime(date) #should perhaps use datetime of [i] (but should not matter alot)
ecipos= jnp.array(ir1['afsGnssStateJ2000'][i][0:3]) #position of satellite in ECI
q = ir1['afsAttitudeState'][i]
quat = R.from_quat(np.roll(q, -1)) #quaterion of satellite pointing
qprime = R.from_quat(ir1['qprime'][i]) #quaterion of channel pointing

localR = np.linalg.norm(sfapi.wgs84.latlon(ir1.TPlat.iloc[i], ir1.TPlon.iloc[i], elevation_m=0).at(t).position.m)
eci_to_ecef=R.from_matrix(itrs.rotation_at(ts.from_datetime(ir1['EXPDate'][i])))


# ypixels = np.linspace(0, ir1['NROW'][i]-1,5).astype('int')
# x, yv = pix_deg(ir1.iloc[i], int(ir1['NCOL'][i]/2), ypixels)



#%%
row = 90

cs_eci_from_row,_ = get_ecivec(ir1.iloc[i],quat)
ecivec = jnp.array(cs_eci_from_row(row)) #LOS of row

eci_to_ecef_mat = jnp.array(eci_to_ecef.as_matrix())
ecef_to_local_mat = jnp.array(ecef_to_local.as_matrix())
r_edges = jnp.array(edges[0])
across_edges = jnp.array(edges[1])
along_edges = jnp.array(edges[2])

# Run the operations to be profiled

start = time.time()
ir1calc = calc_intensity(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat,r_edges,across_edges,along_edges,filters,localR,o2array,[VERarray,Tarray])
end = time.time()
print('time:' + str(end - start))

start = time.time()
ir1calc = calc_intensity(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat,r_edges,across_edges,along_edges,filters,localR,o2array,[VERarray,Tarray])
end = time.time()
print('time:' + str(end - start))


start = time.time()
ir1calc, [vergrad1,tgrad1] = value_and_grad(calc_intensity,argnums=10)(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat,r_edges,across_edges,along_edges,filters,localR,o2array,[VERarray,Tarray])
end = time.time()
print('time:' + str(end - start))



#number = 2
#execution_time_1 = timeit.repeat(lambda: calc_VER(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat,r_edges,across_edges,along_edges,filters,localR,o2array,VERarray,Tarray), number=number, repeat=2)
#print(np.array(execution_time_1)/number)

#execution_time_2 = timeit.repeat(lambda: calc_VER_orig(ecipos,ecivec,eci_to_ecef,ecef_to_local), number=10, repeat=10)
#print(execution_time_2)

#print(np.mean(np.array(execution_time_1)/np.array(execution_time_2)))

#calc_VER(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat)
#%timeit calc_VER(row,ecipos,eci_to_ecef,ecef_to_local,edges,filters,localR,quat,qp,x,yv)


#%%
#calculate jacobian for all measurements


# ks = sp.lil_array((len(tanheights)*len(df)*2,Tarray.size+ VERarray.size)) #determine shape of K

# ir1profiles = []
# ir2profiles = []
# heights = []
# ecipos = []
# ecivecs = []
# posecef_sph=[]
# ir1calcs=[]
# ir2calcs=[]
# k_row = 0


# for i in range(0,len(df)):
#     start_long = time.time()
#     zs, ir1m,ir2m = prepare_measurment(ir1.iloc[i],ir2.iloc[i],ir3.iloc[i],ir4.iloc[i])
#     ir1profiles.append(ir1m)
#     ir2profiles.append(ir2m)
#     heights.append(zs)

#     ecipos.append(df.iloc[i]['afsGnssStateJ2000'][0:3])
#     d = df.iloc[i]['EXPDate']
#     t = ts.from_datetime(d)
#     localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat.iloc[i], df.TPlon.iloc[i], elevation_m=0).at(t).position.m)
#     q = df['afsAttitudeState'].iloc[i]
#     quat = R.from_quat(np.roll(q, -1))

#     #calculate tangent geometries for centre column limited rows and make a spline for it
#     ypixels = np.linspace(0, df['NROW'].iloc[i]-1, 5).astype('int')
#     x, yv = pix_deg(df.iloc[i], int(df['NCOL'].iloc[i]/2), ypixels)
#     qp = R.from_quat(df['qprime'].iloc[i])
#     ecivec = np.zeros((3, len(yv)))
#     k=np.zeros((df['NROW'].iloc[i]-10,len(retrival_heights)))
#     for irow, y in enumerate(yv):
#         los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
#         ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
#     cs_eci=CubicSpline(ypixels,ecivec.T)
#     tanheights2ecivec=CubicSpline(orig_heights[ypixels],ecivec.T)

#     to_ecef=R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][i])))
#     satecefpos=to_ecef.apply(ecipos[-1])
#     satlocalpos=ecef_to_local.apply(satecefpos)
    
#     #calculate jacobian for each row
#     #for irow in range(df['NROW'].iloc[i]-10):
#     stop_2 = time.time()

#     for irow in range(len(tanheights)):

#         calc_VER(tanheights[i],ecef_to_local,edges,filters)
#         start_short = time.time()
#         #ecivec=cs_eci(irow)
#         ecivec=tanheights2ecivec(tanheights[irow])
#         ecefvec=to_ecef.apply(ecivec)
#         localvec=ecef_to_local.apply(ecefvec)
#         zs=np.linalg.norm(ecipos[-1])
#         theta=np.arccos(np.dot(ecipos[-1],ecivec)/zs)
#         b=2*zs*np.cos(theta)
#         root=np.sqrt(b**2+4*((120e3+localR)**2 - zs**2))
#         s_120_1 =(-b-root)/2
#         s_120_2 =(-b+root)/2
#         s_steps = np.arange(s_120_1,s_120_2,steps)
#         stop_3 = time.time()
#         poslocal_i=(np.expand_dims(satlocalpos, axis=0).T+s_steps*np.expand_dims(localvec, axis=0).T).astype('float32')
#         poslocal_i_sph = cart2sph(poslocal_i.T)   
#         poslocal_i_sph=jnp.array(np.asarray(np.array(poslocal_i_sph).T))   
#         stop_4 = time.time()
#         ir1calc=ir1fun(poslocal_i_sph,steps,o2array,[VERarray,Tarray])
#         ir2calc=ir2fun(poslocal_i_sph,steps,o2array,[VERarray,Tarray])
#         stop_5 = time.time()
#         [vergrad1,tgrad1]=ir1grad_non_compiled(poslocal_i_sph,steps,o2array,[VERarray,Tarray])
#         [vergrad2,tgrad2]=ir2grad_non_compiled(poslocal_i_sph,steps,o2array,[VERarray,Tarray])
#         stop_6 = time.time()
#        # ir1calcs.append(ir1calc.item())
#        # ir2calcs.append(ir2calc.item())
#         stop_7 = time.time()
#         hist=np.hstack([vergrad1[:,0,:],tgrad1[:,0,:]])
#         k_ir1 = hist.reshape(-1)
#         hist=np.hstack([vergrad2[:,0,:],tgrad2[:,0,:]])
#         k_ir2 = hist.reshape(-1)
#         stop_8 = time.time()
#         ks[k_row,:] = k_ir1
#         ks[k_row+len(tanheights)*len(df),:] = k_ir2
#         k_row = k_row+1
#         stop_short = time.time()
#         print("Time for one pixel: " + str(stop_short - start_short))
#         print("Time for one pixel, step 3: " + str(stop_3 - start_short))
#         print("Time for one pixel, step 4: " + str(stop_4 - stop_3))
#         print("Time for one pixel, step 5: " + str(stop_5 - stop_4))
#         print("Time for one pixel, step 6: " + str(stop_6 - stop_5))
#         print("Time for one pixel, step 7: " + str(stop_7 - stop_6))
#         print("Time for one pixel, step 8: " + str(stop_8 - stop_7))


#     stop_long = time.time()
#     print("Time for one profile: " + str(stop_long - start_long))
#     print("Time for one profile, step 1: " + str(stop_1 - start_long))
#     print("Time for one profile, step 2: " + str(stop_2 - stop_1))


# # with open("runningfile_2", "wb") as file:
# #     pickle.dump((ir1calcs,ir2calcs,ir1profiles,ir2profiles,ks,Tarray,VERarray), file)
    
        
    
# # # ecivecs= np.reshape(ecivecs,(len(df),-1,3))
# # #posecef_sph= np.reshape(posecef_sph,(-1,3))
# # #%%

# # z= np.array(heights).mean(axis=0)
# # inputdata = xr.Dataset({
# #     'time': (['time'], df.EXPDate),
# #     'channel': (['time'], df.channel),
# #     'satlat': (['time'], df.satlat, {'long_name': "MATS' Latitude", 'units': 'deg'}),
# #     'satlon': (['time'], df.satlon),
# #     'TPlat': (['time'], df.TPlat, {'long_name': "TP Latitude", 'units': 'deg'}),
# #     'TPlon': (['time'], df.TPlon),
# #     'TPsza': (['time'], df.TPsza),
# #     'TPssa': (['time'], df.TPssa),
# #     'ecipos': (['time', 'xyz'], ecipos),   
# #     'z':  (['z'], z, {'long_name': 'Approx Altitude', 'units': 'm'}),
# #     'ir1profile': (['time', 'z'], np.asarray(ir1profiles), {'long_name': 'LOS  intensity', 'units': 'Photons m-2 nm-1 sr-1 s-1'}),
# #     'ir2profile': (['time', 'z'], np.asarray(ir2profiles), {'long_name': 'LOS  intensity', 'units': 'Photons m-2 nm-1 sr-1 s-1'}),
# #     'heights': (['time', 'z'],  heights),
# #     'ret_grid_z': (['z_r'], edges[0]),
# #     'ret_grid_lon': (['lon_r'], edges[1]),
# #     'ret_grid_lat': (['lat_r'], edges[2]),
# # })

# # inputdata.to_netcdf('IR1IR2test_400-520.nc')
# # # %%
# # filename = "jacobianIR1IR2_400-520.pkl"
# # with open(filename, "wb") as file:
# #     pickle.dump((edges, ks, ecef_to_local), file)
# # #%%
# # filename = "intensityIR1IR2_400-520.pkl"
# # with open(filename, "wb") as file:
# #     pickle.dump((ir1calcs,ir2calcs), file)

# # # %%
# # z1,p1,orig_heights=prepare_profile(ir1.iloc[i])
# # z2,p2,_=prepare_profile(ir2.iloc[i])
# # z3,p3,_=prepare_profile(ir3.iloc[i])
# # z4,p4,_=prepare_profile(ir4.iloc[i])
# # plt.figure()
# # plt.semilogx(p1,z1,p2,z2,p3,z3,p4,z4)

# # p3=p3-p3[-4:].mean()/1.05
# # p4=p4-p4[-4:].mean()/1.05
# # plt.semilogx(p3,z3,p4,z4)
# # sc=2
# # p1=p1-(p3+p4)/2/sc
# # p2=p2-(p3+p4)/2/sc
# # #p1=p1-p1[-4:-1].mean()/1.1
# # #p2=p2-p2[-4:-1].mean()/1.1
# # plt.semilogx(p1,z3,p2,z4)
# # plt.legend(['IR1','IR2','IR3','IR4','IR3c','IR4c','IR1c','IR2c',])
# # plt.xlim([5e12,1e15])
# # plt.show()
# # plt.figure()
# # plt.plot(3.57/8.16*p1/p2,z1)
# # plt.xlim([0.1,0.8])
# #  # %%
# # with open("runningfile_1", "rb") as file:
# #     [i,irow,ir1calcs,ir2calcs,profiles,ks] = pickle.load(file)
# # # %%
# # plt.figure()
# # plt.pcolor(vergrad1[:,0,:])
# # plt.colorbar()
# # #%%
# # plt.figure()
# # plt.plot(ir1calcs)
# # #%%
# %%
