# %% create list of strips.
import datetime as DT
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
from Keogram import CenterStrip
import numpy as np
from aurora_analysis import set_aurora_spec, save_strips

start_time = DT.datetime(2023,2,15,0,0,0)
stop_time = DT.datetime(2023,2,17,0,0,0)
numdays = stop_time-start_time #number of days
Tperiod = timedelta(minutes=100)
items = pd.read_pickle('15to16febIR1')


# %%
def all_strips(items, numdays, Tperiod):
    "returns all strips, saves as panda list"
    n = 0
    orb = 0
    centercol = 22
    airglowlim = 160
    auroramean = 50
    strips = []
    
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        
        #this for loop goes through the images starting from the end of previous day
        startday = items.iloc[orb].EXPDate
        for i in range(n, len(items)-1):
            
            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6 and i < len(items)-2:
                continue
            else:  #if this is True, next image will belong to next orbit.                         
                #print(i,'orbitstart', items.iloc[i].EXPDate)
                if items.iloc[i].TPlat > 0: #north hemisphere
                    auroraintensity = 55
                    #creates orbit from index n to i
                    #print('NH',orb,i)
                    NH = items.iloc[orb:i+1]
                    if len(NH) == 0 : #if empty, go to next hemisphere
                        continue
                    
                    # This for loop goes through the images belonging to NH
                    for k, ccd in NH.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip = ccdimage[:,centercol]

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim

                        #if row >= airglowlim + 5 and ccd_strip.item[row] > auroraintensity:
                        top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

                        #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        new_strip = CenterStrip(ccd)
                        new_strip.makeVerticalStrip() 

                        if ccd_strip.item(top_max) > auroraintensity: #check so we have aurora above row.
                            if top_mean > auroramean:
                                #sets the position coordinates of the max intensity point of strips with aurora
                                set_aurora_spec(new_strip,ccd,row,centercol)
                        strips.append(new_strip)     

                elif items.iloc[i].TPlat < 0: #south hemisphere
                    auroraintensity = 70
                    #print('SH',orb,i)
                    SH = items.iloc[orb:i+1]
                    if len(SH) == 0 :
                        continue
                    for k, ccd in SH.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip = ccdimage[:,centercol]

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
                        
                        top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

                        #gives the row of the maximum 10 rows above the limit to check for aurora there.
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        new_strip = CenterStrip(ccd)
                        new_strip.makeVerticalStrip()

                        if ccd_strip.item(top_max) > auroraintensity:
                            if items.iloc[i].TPlat > -50 and items.iloc[i].TPlon > -90 and items.iloc[i].TPlon < 40:
                                pass
                            elif top_mean > auroramean:    
                                set_aurora_spec(new_strip,ccd,row,centercol)
                        strips.append(new_strip) 
               
                orb = i+1   #start number for next orbit
                nextday_startdate = items.iloc[orb].EXPDate
                #print('nextday', nextday_startdate)
                #comparing the day at start of the new orbit with the active orbits start.
                if startday.day != nextday_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day", startday)
                    n = orb
                    break
    save_strips(strips,'allstrips.mat','allstrips')
    return
    

def get_stripRow(strips):
    "Get all rows from given list of strip objects"
    allrows = strips['row']
    return allrows

# %%
