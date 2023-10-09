# %% create list of strips.
import datetime as DT
import pandas as pd
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
from datetime import timedelta
from sklearn.linear_model import LinearRegression
# setting path
sys.path.append(r'C:\Users\ceona\Documents\GitHub\MATS-analysis\MATS-analysis\Ceona')
from Keogram import makeStripMatrix
from aurora_analysis import set_aurora_spec, save_strips

# %% 
def gradientmatrix(items,airglowlim):
    """Removes the gradient background from airglowlim and up, updates the strips corresponding to that keogram
    Returns updated keogram matrix and list of strips"""
    matrix, striplist = makeStripMatrix(items)
    newmatrix = matrix.copy()
    for rowindex in range(airglowlim,matrix.shape[0]):
        y = matrix[rowindex,:] 
        x = np.arange(0,len(y)).reshape((-1, 1))
        
        model = LinearRegression().fit(x,y)
        y_lin = model.predict(x)
        newmatrix[rowindex,:] = y-y_lin

    #Update the strips from the new matrix
    for col, strip in enumerate(striplist):
        strip.strip = newmatrix[:,col]

    return newmatrix, striplist
    
def strips(items, numdays):
    Tperiod = timedelta(minutes=100)
    airglowlim = 140
    auroramean = 20
    n = 0
    strips = []
    stripsNH = []
    stripsSH = []
    aurorastrips = []

    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        #this for loop goes through the images starting from the end of previous day
        for i in range(n, len(items)-1):
            
            orbit_startdate = items.iloc[n].EXPDate
            
            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6 and i < len(items)-2: 
                continue
            else:                          
                #creates orbit from index n to i
                if items.iloc[i].TPlat > 0: #north hemisphere
                    auroraintensity = 20
                    NH = items.iloc[n:i+1]
                    if len(NH) == 0 :
                        continue
                   
                    #gets the removed gradient matrix corresponding to that hemisphere
                    matrix, striplist = gradientmatrix(NH, airglowlim)
                    for m, strip in enumerate(striplist):
                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(strip.strip[airglowlim:]) + airglowlim
                        
                        top_mean = np.sum(strip.strip[airglowlim+10:])/len(strip.strip[airglowlim+10:])
                        #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
                        top_max = np.argmax(strip.strip[airglowlim+10:]) + airglowlim + 10        
 
                        if strip.strip.item(top_max) >= auroraintensity: #check so we have aurora above row.
                            if top_mean > auroramean:
                                #sets the position coordinates of the max intensity point of strips with aurora
                                print('Row',row,'RowI',strip.strip.item(row),'Topmax',strip.strip.item(top_max),'Mean',top_mean,strip.time)

                                ccd = NH.iloc[m]
                                set_aurora_spec(strip,ccd,row)
                                aurorastrips.append(strip)
                        strips.append(strip)
                        stripsNH.append(strip)  
                elif items.iloc[i].TPlat < 0: #south hemisphere
                    auroraintensity = 20
                    SH = items.iloc[n:i+1]
                    if len(SH) == 0 :
                        continue
                    #gets the removed gradient matrix corresponding to that hemisphere
                    matrix, striplist = gradientmatrix(SH,airglowlim)
                    for m, strip in enumerate(striplist):

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(strip.strip[airglowlim:]) + airglowlim
                        top_mean = np.sum(strip.strip[airglowlim+10:])/len(strip.strip[airglowlim+10:])

                        #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
                        top_max = np.argmax(strip.strip[airglowlim+10:]) + airglowlim + 10   
                        
                        if SH.iloc[m].TPlat >= -60 and SH.iloc[m].TPlon >= 0 and SH.iloc[m].TPlon <= 40 or (SH.iloc[m].TPlat >= -60 and SH.iloc[m].TPlon >= 300):
                            print(strip.latitude, SH.iloc[m].TPlon, strip.time)
                            pass     
                        elif strip.strip.item(top_max) > auroraintensity:  #satlat to avoid SAA 
                            if top_mean > auroramean:
                                #print('Row',row,'RowI',strip.strip.item(row),'Topmax',strip.strip.item(top_max),'Mean',top_mean,strip.time)
  
                                ccd = SH.iloc[m]
                                set_aurora_spec(strip,ccd,row)
                                #sets the position coordinates of the max intensity point of strips with aurora
                                aurorastrips.append(strip)
                        strips.append(strip)
                        stripsSH.append(strip)  
                n = i+1 #start number for next orbit, to avoid unnecessary search
                nextorbit_startdate = items.iloc[n].EXPDate
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day")
                    print(orbit_startdate)       
                    break
    save_strips(strips,'allstripstest.mat','allstrips')
    save_strips(aurorastrips,'aurorastripstest.mat','aurorastrips')

    return 
# %%
#linear regression test for short interval (only one orbit)
testitems = pd.read_pickle(r'C:\Users\ceona\Documents\GitHub\MATS-analysis\MATS-analysis\Ceona\3marstest')

def test(testitems):
    airglowlim = 140
    matrix, striplist = makeStripMatrix(testitems)
    newmatrix = matrix.copy()
    for rowindex in range(airglowlim,matrix.shape[0]):
        y = matrix[rowindex,:] 
        x = np.arange(0,len(y)).reshape((-1, 1))
        
        model = LinearRegression().fit(x,y)
        y_lin = model.predict(x)
        newmatrix[rowindex,:] = y-y_lin

        #if rowindex == airglowlim:
            #plt.plot(y)
            #plt.plot(y_lin)
            #plt.plot(y-y_lin)

    #Update the strips from the new matrix
    for col, strip in enumerate(striplist):
        strip.strip = newmatrix[:,col]
    
    scipy.io.savemat('keogramgrad.mat',{'keogramgrad': newmatrix, 'label':'pixel'}) #saves to matlabfile
    #plt.pcolormesh(newmatrix)
    print(newmatrix[:,0] -striplist[0].strip )
    return

# %%
def Main():
    start_time = DT.datetime(2023,2,15,00,0,0)
    stop_time = DT.datetime(2023,2,17,00,0,0)
    numdays = stop_time-start_time
    items = pd.read_pickle(r'C:\Users\ceona\Documents\GitHub\MATS-analysis\MATS-analysis\Ceona\15to16febIR1')
    strips(items,numdays)
    return
