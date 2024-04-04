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


# def prepare_measurment(ir1,ir2,ir3,ir4):
#     # Subtracts background channels from IR1 and IR2
#     # Converted to band strength
#     common_heights,p1,_=prepare_profile(ir1)
#     _,p2,_=prepare_profile(ir2)
#     _,p3,_=prepare_profile(ir3)
#     _,p4,_=prepare_profile(ir4)
#     #plt.figure()
#     #plt.semilogx(p1,z1,p2,z2,p3,z3,p4,z4)
#     #If it is dark: Return altitudes and 
#     if ir3.TPsza>98 : 
#         raise ValueError('Nightglow not supported')
#         #return common_heights,(p1-(p3+p4)/2)/(p2-(p3+p4)/2)*3.57/8.16
    
#     p3=p3-p3[-4:].mean()/1.05
#     p4=p4-p4[-4:].mean()/1.05
#     p1=p1-(p3+p4)/2
#     p2=p2-(p3+p4)/2
#     IR1_BW = 3.57 #nm
#     IR2_BW = 8.16 #nm
#     ir1_photons = np.array(p1*IR1_BW)
#     ir2_photons = np.array(p2*IR2_BW)
#     return common_heights,ir1_photons,ir2_photons

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

@jit
def interpT(x, xs, ys):
    ix = jnp.floor(x - 100).astype(jnp.int32)
    diff = xs[ix] - xs[ix - 1]
    return ((x - xs[ix - 1]) * ys[ix, :].T + (xs[ix] - x) * ys[ix - 1, :].T) / diff

@jit
def interppos_jit(radius, latitude, r_edges, along_edges):
    iz_values = jnp.floor((radius - r_edges[0]) / ((r_edges[1:] - r_edges[:-1]).sum() / (len(r_edges) - 1)))
    ix_values = jnp.floor((latitude - along_edges[0]) / ((along_edges[1:] - along_edges[:-1]).sum() / (len(along_edges) - 1)))

    iz = iz_values.astype(jnp.int32)
    ix = ix_values.astype(jnp.int32)

    return iz, ix

def interppos(pos, edges):
    iz = jnp.array(jnp.floor((pos[:, 0] - edges[0][0]) / jnp.diff(edges[0]).mean()), dtype=int)
    ix = jnp.array(jnp.floor((pos[:, 2] - edges[2][0]) / jnp.diff(edges[2]).mean()), dtype=int)

    return iz,ix

@jit
def ir1fun(radius, across,along,step,r_edges,across_edges,along_edges,filters,o2array,VER,Temps):
    startT= jnp.linspace(100,600,501) #setting up temp grid
    iz,ix=interppos_jit(radius,along,r_edges,along_edges) #get temp at each point along path
    pathtemps=Temps[iz,0,ix]
    pathtemps = jnp.array(pathtemps)
    sigmas=interpT(pathtemps,startT,sigma) #sigma for each point along the path
    emissions=interpT(pathtemps,startT,emission) #normalized emission spectrum for each point along the path
    
    o2=o2array[iz,0,ix]
    o2 = jnp.array(o2)

    tau = (sigmas*o2).cumsum(axis=1)*step*1e2

    VERs=VER[iz,0,ix]
    VERs = jnp.array(VERs)

    VERs = VERs*step*1e2
    res=filters@(jnp.exp(-tau)*VERs*emissions) #get limb spectral radiance after filter at tanheight for VER (see other place)
    result = res[0].sum() #sum over wavelengths for filter 1

    return result

@jit
def ir2fun(radius, across,along,step,r_edges,across_edges,along_edges,filters,o2array,VER,Temps):
    startT= jnp.linspace(100,600,501) #setting up temp grid
    iz,ix=interppos_jit(radius,along,r_edges,along_edges) #get temp at each point along path
    pathtemps=Temps[iz,0,ix]
    pathtemps = jnp.array(pathtemps)
    sigmas=interpT(pathtemps,startT,sigma) #sigma for each point along the path
    emissions=interpT(pathtemps,startT,emission) #normalized emission spectrum for each point along the path
    
    o2=o2array[iz,0,ix]
    o2 = jnp.array(o2)

    tau = (sigmas*o2).cumsum(axis=1)*step*1e2

    VERs=VER[iz,0,ix]
    VERs = jnp.array(VERs)

    VERs = VERs*step*1e2
    res=filters@(jnp.exp(-tau)*VERs*emissions) #get limb spectral radiance after filter at tanheight for VER (see other place)
    result = res[1].sum() #sum over wavelengths for filter 1

    return result


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

def get_geometry_for_image(df_row):
    t = ts.from_datetime(df_row['EXPDate']) #should perhaps use datetime of [i] (but should not matter alot)
    ecipos= jnp.array(df_row['afsGnssStateJ2000'][0:3]) #position of satellite in ECI
    q = df_row['afsAttitudeState']
    quat = R.from_quat(np.roll(q, -1)) #quaterion of satellite pointing
    qprime = R.from_quat(df_row['qprime']) #quaterion of channel pointing

    localR = np.linalg.norm(sfapi.wgs84.latlon(df_row.TPlat, df_row.TPlon, elevation_m=0).at(t).position.m)
    eci_to_ecef=R.from_matrix(itrs.rotation_at(t))

    ypixels = np.linspace(0, df_row['NROW']-1,5).astype('int')
    x, yv = pix_deg(df_row, int(df_row['NCOL']/2), ypixels)

    return quat,qprime,x,yv,ypixels,eci_to_ecef,ecipos,localR

#%%
## Run code ##

ts = sfapi.load.timescale()
ir1, ir2,_,_,_,_ = get_instrument_data()
filters = jnp.array(get_filtercurves())

radius,across_track,along_track,ret_lats,ret_lons,date,edges,tanheights,ecef_to_local = generate_grid(ir1)
r_edges = jnp.array(edges[0])
across_edges = jnp.array(edges[1])
along_edges = jnp.array(edges[2])
ecef_to_local_mat = jnp.array(ecef_to_local.as_matrix())

msis, sigma,emission = get_background_data()

Tarray,o2array,VERarray = get_bg_atmosphere_on_grid(radius,across_track,along_track,ret_lats,ret_lons,date,msis)
Tarray = jnp.array(Tarray)
o2array = jnp.array(o2array)
VERarray = jnp.array(VERarray)

#for one image
i=2
quat,qprime,x,yv,ypixels,eci_to_ecef,ecipos,localR = get_geometry_for_image(ir1.iloc[i])
eci_to_ecef_mat = jnp.array(eci_to_ecef.as_matrix())


#for one row(pixel)
row = 90
cs_eci_from_row,_ = get_ecivec_fast(quat,qprime,x,yv,ypixels)
ecivec = jnp.array(cs_eci_from_row(row)) #LOS of row

# calc jacobian
start = time.time()
ir1calc, [vergrad1,tgrad1] = value_and_grad(calc_intensity,argnums=10)(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat,r_edges,across_edges,along_edges,filters,localR,o2array,[VERarray,Tarray])
end = time.time()
print('time:' + str(end - start))

start = time.time()
ir1calc, [vergrad1,tgrad1] = value_and_grad(calc_intensity,argnums=10)(ecipos,ecivec,eci_to_ecef_mat,ecef_to_local_mat,r_edges,across_edges,along_edges,filters,localR,o2array,[VERarray,Tarray])
end = time.time()
print('time:' + str(end - start))