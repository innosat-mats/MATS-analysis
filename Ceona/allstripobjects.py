# %% create list of all strips.
import datetime as DT
import pandas as pd
from Keogram import CenterStrip
import numpy as np
from aurora_analysis import set_aurora_spec, save_strips
from skyfield.api import load
from skyfield.units import Distance
from skyfield.toposlib import wgs84
from skyfield.positionlib import Geocentric
from aacgmv2 import get_aacgm_coord
import scipy
import time

start_time = DT.datetime(2023,3,15,0,0,0)
stop_time = DT.datetime(2023,3,22,0,0,0)
numdays = stop_time-start_time #number of days
items = pd.read_pickle('15to21febIR1')

# %%
def TPpos(ccditem):
    """Function giving the GPS position in mlat, mlon, mlt"""
    ecipos= ccditem['afsGnssStateJ2000'][0: 3]
    d = ccditem['EXPDate']
    ts= load.timescale()
    t = ts.from_datetime(d)
    satpo = Geocentric(position_au=Distance(
        m=ecipos).au, t=t)
    position = wgs84.geographic_position_of(satpo)
    TPalt = position.elevation.km
    TPlat = position.latitude.degrees
    TPlon = position.longitude.degrees
    mlat, mlon, mlt = get_aacgm_coord(TPlat,TPlon,TPalt,ccditem.EXPDate, method='ALLOWTRACE')
    return mlat, mlon, mlt

# %%
def all_strips(items):
    "returns all strips, saves as panda object list"
    airglowlim = 160
    auroramean = 50
    strips = []
    stripsNH = []
    stripsSH = []
    
    for i in range(0, len(items)-1):
        ccd = items.iloc[i]
        new_strip = CenterStrip(ccd)
        new_strip.makeVerticalStrip()
        ccd_strip = new_strip.strip 
        print(i)
        if ccd.TPlat > 0: #north hemisphere
            auroraintensity = 55
            
            #finds the row of the max intensity value of each strip, above airglow limit
            row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
            top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

            #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
            top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                    
            if ccd_strip.item(top_max) >= auroraintensity: #check so we have aurora above row.
                if top_mean > auroramean:
                    #sets the position coordinates of the max intensity point of strips with aurora
                    set_aurora_spec(new_strip,ccd,row)

            strips.append(new_strip)
            stripsNH.append(new_strip)  

        elif ccd.TPlat < 0: #south hemisphere
            auroraintensity = 70

            #finds the row of the max intensity value of each strip, above airglow limit
            row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
            top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

            #gives the row of the maximum 10 rows above the limit to check for aurora there.
            top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
            
            if ccd_strip.item(top_max) > auroraintensity:  #satlat to avoid SAA
                if items.iloc[i].satlat > -50 and items.iloc[i].satlon > -90 and items.iloc[i].satlon < 40:
                    pass
                elif top_mean > auroramean:    
                    set_aurora_spec(new_strip,ccd,row)
            
            strips.append(new_strip)
            stripsSH.append(new_strip)  
    save_strips(strips,'allstrips2weekmar.mat','allstrips')
    #save_strips(stripsNH,'allstripsNH2weekfeb.mat','allstripsNH')
    #save_strips(stripsSH,'allstripsSH2weekfeb.mat','allstripsSH')
    return

def sattelite_MLT(items):
    MLT = []
    Mlat = []
    for k, ccd in items.iterrows():
        st = time.time()
        mlat, mlon, mlt = TPpos(ccd)
        print(time.time()-st)
        MLT.append(mlt)
        Mlat.append(mlat)
    
    scipy.io.savemat('MLT.mat',{'MLT': MLT, 'label':'mlt'}) #saves to matlabfile
    scipy.io.savemat('MLat.mat',{'Mlat': Mlat, 'label':'mlat'}) #saves to matlabfile
    return

# %%
