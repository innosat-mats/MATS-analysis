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

#%%
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

def prepare_measurment(ch,ir3,ir4,subtract=True,endcut=-25):

    z1,p1=prepare_profile(ch)
    z3,p3=prepare_profile(ir3)
    z4,p4=prepare_profile(ir4)

    if subtract:

        p3 = p3[10:endcut]
        p4 = p4[10:endcut]

        p3=p3-p3[-4:].mean()/1.05
        p4=p4-p4[-4:].mean()/1.05

        p1=p1[10:endcut]-(p3+p4)/2

    return np.array(z1),np.array(p1)


def prepare_profile(ch):
    # TBD comment
    # calculates mean profile and tan heights
    image = image = np.stack(ch.ImageCalibrated)/ch.TEXPMS*1000
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))

    profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e13)
    common_heights = np.arange(60000,110250,250)
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile

def generate_ks(df,ir3,ir4,retrival_heights, s_steps, splinedfactor,abs,subtract,endcut):
    # TBD comment
    # calculates ks from mean profile
    ts = sfapi.load.timescale()
    profiles = []
    heights = []
    ecipos = []
    ecivecs = []
    ks=[]

    for i in range(len(df)):
        print(f'{i+1}/{len(df)}')
        zs, p = prepare_measurment(df.iloc[i], ir3.iloc[i], ir4.iloc[i],subtract,endcut=-25)
        #zs, p = prepare_profile(df.iloc[i])

        # remove some pixels top and bottom
        if df.iloc[i]['channel'] == 'IR1':
            zs = zs[10:endcut]
            p = p*3.57
            nrows = len(zs)
        
        if df.iloc[i]['channel'] == 'IR2':
            zs = zs[10:endcut]
            p = p*8.17
            nrows = len(zs)

        profiles.append(p)
        heights.append(zs)
        ecipos.append(df['afsGnssStateJ2000'][i][0:3])
        d = df['EXPDate'][i]
        t = ts.from_datetime(d)
        localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat[i], df.TPlon[i], elevation_m=0).at(t).position.m)
        q = df['afsAttitudeState'][i]
        
        quat = R.from_quat(np.roll(q, -1))
        ypixels = np.linspace(0, nrows, 5)
        x, yv = pix_deg(df.iloc[i], int(df['NCOL'][i]/2), ypixels)
        qp = R.from_quat(df['qprime'][i])
        ecivec = np.zeros((3, len(yv)))
        k=np.zeros((nrows,len(retrival_heights)))
        for irow, y in enumerate(yv):
            los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
            ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
        cs_eci=CubicSpline(ypixels,ecivec.T)
        for irow in range(nrows):
            ecivec=cs_eci(irow)
            pos=np.expand_dims(ecipos[-1], axis=0).T+s_steps*np.expand_dims(ecivec, axis=0).T
            point_height=np.linalg.norm(pos,axis=0)-localR
            
            # absorption weight
            if abs:
                minarg=point_height[:].argmin() # pos of lowest tan height
                target_tangent=zs[irow]/1000
                distances=(s_steps-s_steps[minarg])/1000 #to km for intepolation
                lowertan=bisect_left(tanz,target_tangent)-1
                uppertan=lowertan+1
                lowerfactor=splinedfactor[lowertan](distances)
                upperfactor=splinedfactor[uppertan](distances)
                newfactor=interp1d([tanz[lowertan],tanz[uppertan]],np.array([lowerfactor,upperfactor]).T)
                weight=np.exp(newfactor(target_tangent))

                #if irow > 150:
                #    plt.plot(weight)
                #    print(np.sum(weight)())

                # calculate k
                counts,bins=np.histogram(point_height/1000,np.hstack((retrival_heights,retrival_heights[-1]+1)),weights=weight)
            else:
                counts,bins=np.histogram(point_height/1000,np.hstack((retrival_heights,retrival_heights[-1]+1)))


            k[irow,:]=counts*steps
            ecivecs.append(ecivec)

        ks.append(k)
    ecivecs= np.reshape(ecivecs,(len(df),-1,3))

    return ks, profiles, ecipos, ecivecs, heights

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
# load images
channels=['IR1', 'IR2']
l1b_version="0.5"
oned_version="0.1"
# dayglow stuff
ascending=True
# end day
start_day=6
end_day=7
# dayglow (ALL SZA)
dmin,dmax = 0, 95
tplat0, tplat1 = -70, 70

for channel in channels:
    for day in range(start_day,end_day):

        print(f'------------------------------------------------------------')
        print(f'DAY: {day} -- CHANNEL: {channel} -- STARTING DAY {start_day} -- ENDING DAY {end_day}')
        print(f'------------------------------------------------------------')
        starttime=datetime(2023,3,day,0,0)
        stoptime=datetime(2023,3,day,23,59)
        #stoptime=datetime(2023,3,day,2,0) ########### NOT 24 !!!!!!
        
        try:

            dftop=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1]})

            if ascending:
                midday=DT.time(12, 0, 0)
                dftop['TPlocaltime'] = pd.to_datetime(dftop['TPlocaltime'])
                dftop = dftop[(dftop['TPlocaltime'].dt.time > midday)]

            df = dftop[dftop['channel'] == channel].dropna().reset_index()#[0:10]
            ir3 = dftop[dftop['channel'] == 'IR3'].dropna().reset_index()#[0:10]
            ir4 = dftop[dftop['channel'] == 'IR4'].dropna().reset_index()#[0:10]

            # absorption weights (TBD: call generation of these)
            abs=True #WRONG DATE  ---------------------------------------
            tanz, splinedfactor=np.load(f"splinedlogfactors{channel}_sep1730_new.npy",allow_pickle=True)

            # select part of orbit
            offset = 10
            num_profiles = len(df)-50 #use 50 profiles for inversion
            print(num_profiles)
            df = df.loc[offset:offset+num_profiles]
            df = df.reset_index(drop=True)

            ir3 = ir3.loc[offset:offset+num_profiles]
            ir3 = ir3.reset_index(drop=True)
            ir4 = ir4.loc[offset:offset+num_profiles]
            ir4 = ir4.reset_index(drop=True)

            # retrieval settings
            retrival_heights= np.arange(60,100,1)
            retrival_heights= np.arange(60,105,1)
            s_140=1700e3
            steps=100 #m steps
            s_steps = np.arange(s_140,s_140 + 2e6,steps)
            abs_bool=True

            # Generate ks // save netcdf
            ks, profiles, ecipos, ecivecs, heights = generate_ks(df,ir3,ir4,retrival_heights,
                                                                s_steps, splinedfactor,
                                                                abs=abs_bool,subtract=True,
                                                                endcut=-25)
            z=np.array(heights).mean(axis=0)

            ch = xr.Dataset({
                'time': (['time'], df.EXPDate),
                'channel': (['time'], df.channel),
                'satlat': (['time'], df.satlat, {'long_name': "MATS' Latitude", 'units': 'deg'}),
                'satlon': (['time'], df.satlon),
                'TPlat': (['time'], df.TPlat, {'long_name': "TP Latitude", 'units': 'deg'}),
                'TPlon': (['time'], df.TPlon),
                'TPsza': (['time'], df.TPsza),
                'TPssa': (['time'], df.TPssa),
                'z':  (['z'], z, {'long_name': 'Approx Altitude', 'units': 'm'}),
                'ecipos': (['time', 'xyz'], ecipos),
                'ecivec': (['time', 'z', 'xyz'], ecivecs),
                'profile': (['time', 'z'], profiles, {'long_name': 'LOS  intensity', 'units': 'Photons m-2 nm-1 sr-1 s-1'}),
                'heights': (['time', 'z'], heights),
                'z_r': (['z_r'], retrival_heights),
                'k': (['time','z','z_r'],ks),
            })

            # save the results
            ch.to_netcdf(f'/media/waves/AVAGO/data/MATS/1D_inversions/v{oned_version}/IR_comparison/{channel}_Mar{day}_l1bv{l1b_version}_sepabs.nc')
            
        except KeyboardInterrupt:
            sys.exit()  
        
        except Warning:
            print('No data!') 


# %%
