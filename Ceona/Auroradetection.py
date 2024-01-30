# %% detect aurora using keograms with removed gradient background.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
import datetime as DT
#from sklearn.linear_model import LinearRegression
import scipy
from Keogram import makeStripMatrix
from aurora_analysis import set_aurora_spec, save_strips, get_aurora_max, set_strip_spec, save_TPMLT

# %% 
def gradientmatrix(items,airglowlim, auroralim):
    """Removes the gradient background from airglowlim and up, updates the strips corresponding to that keogram
    Returns updated keogram matrix and its list of strip objects"""
    matrix, striplist = makeStripMatrix(items)
    newmatrix = matrix.copy()

    for rowindex in range(airglowlim,matrix.shape[0]):
        y = matrix[rowindex,:] 
        x = np.arange(0,len(y))

        #Polynomial regression model
        if rowindex >= auroralim:
            #weird intense spikes above auroralim, like the moon pass, will be ignored in the fit
            #for April week use y < 400.
            valid_indices = np.where(y < 400)[0]
            y_filt = y[valid_indices]
            x_filt = np.arange(0,len(y_filt))
            model = np.poly1d(np.polyfit(x_filt, y_filt, 3))
        else:
            valid_indices = np.where(y < 600)[0]
            y_filt = y[valid_indices]
            x_filt = np.arange(0,len(y_filt))
            model = np.poly1d(np.polyfit(x_filt, y_filt, 3))
        
        y_reg = model(x) #new y_points
        newmatrix[rowindex,:] = y-y_reg

    #Update the strips from the new matrix
    for col, strip in enumerate(striplist):
        strip.strip = newmatrix[:,col]

    return newmatrix, striplist
    
#Because all the strips are needed for the polynomial regression of each hemisphere keogram matrix 
#to get the gradient removed matrix
def get_strips(items, numdays,filedate):
    """Returns a list of all strips and list of only aurorastrips, using linear regression and set aurora conditions"""
    Tperiod = timedelta(minutes=100)
    airglowlim = 130
    n = 0
    allstripsNH = []
    allstripsSH = []
    allstrips = []   #used for plotting red dots on onverview
    aurorastripsNH = []
    aurorastripsSH = []
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
                    auroralim = 150
                    auroramean = 15
                    auroraintensity = 18
                    NH = items.iloc[n:i+1]  #needs the plus 1 for iloc function
                    if len(NH) == 0 :
                        continue
                   
                    #gets the removed gradient matrix corresponding to that hemisphere
                    matrix, striplist = gradientmatrix(NH, airglowlim,auroralim)
                    for m, strip in enumerate(striplist):
                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(strip.strip[auroralim:]) + auroralim
                        
                        #top_mean = np.sum(strip.strip[airglowlim+10:])/len(strip.strip[airglowlim+10:])
                        mean = np.sum(strip.strip[auroralim:])/len(strip.strip[auroralim:])

                        #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
                        top_max = np.argmax(strip.strip[auroralim+10:]) + auroralim + 10        
                        ccd = NH.iloc[m]
                        if strip.strip.item(top_max) >= auroraintensity and mean > auroramean: #check so we have aurora above row.
                            #sets the position coordinates of the max intensity point of strips with aurora
                            #print(strip.latitude, NH.iloc[m].TPlon, strip.time)  
                            #print('Row',row,'RowI',strip.strip.item(row),'Topmax',strip.strip.item(top_max),'Mean',top_mean,strip.time)
                            set_aurora_spec(strip,ccd,row)
                            aurorastrips.append(strip)
                            aurorastripsNH.append(strip)
                        else:
                            set_strip_spec(strip,ccd)

                        allstripsNH.append(strip)
                        allstrips.append(strip)

                elif items.iloc[i].TPlat < 0: #south hemisphere
                    auroralim = 120    #150 for pol-reg, 120 for normal keogram
                    auroramean = 30     #15 for pol-reg, 30 for normal keogram
                    auroraintensity = 40   #18 used for pol-reg, 40 for normal
                    SH = items.iloc[n:i+1]
                    if len(SH) == 0 :
                        continue
                    #gets the removed gradient matrix corresponding to that hemisphere
                    #matrix, striplist = gradientmatrix(SH,airglowlim,auroralim)
                    matrix, striplist = makeStripMatrix(SH)   #Used for april week 2-4 and may, when background is more linear

                    for m, strip in enumerate(striplist):

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(strip.strip[auroralim:]) + auroralim
                        #top_mean = np.sum(strip.strip[airglowlim+10:])/len(strip.strip[airglowlim+10:])
                        mean = np.sum(strip.strip[auroralim:])/len(strip.strip[auroralim:])

                        #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
                        top_max = np.argmax(strip.strip[auroralim+10:]) + auroralim + 10   
                        ccd = SH.iloc[m]
                        #to avoid strips from SAA
                        if SH.iloc[m].TPlat >= -60 and SH.iloc[m].TPlon >= 0 and SH.iloc[m].TPlon <= 40 or (SH.iloc[m].TPlat >= -60 and SH.iloc[m].TPlon >= 290):
                            set_strip_spec(strip,ccd)    
                        elif strip.strip.item(top_max) > auroraintensity and mean > auroramean:
                                #print(strip.latitude, SH.iloc[m].TPlon, strip.time)
                                print('Row',row,'RowI',strip.strip.item(row),'Topmax',strip.strip.item(top_max),'Mean',mean,strip.time)
                                set_aurora_spec(strip,ccd,row)
                                #sets the position coordinates of the max intensity point of strips with aurora
                                aurorastrips.append(strip)
                                aurorastripsSH.append(strip)
                        else:
                            #if not aurora set magnetic coordinates of TP
                            set_strip_spec(strip,ccd)

                        allstripsSH.append(strip)
                        allstrips.append(strip)
                n = i+1 #start number for next orbit, to avoid unnecessary search
                nextorbit_startdate = items.iloc[n].EXPDate
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day", orbit_startdate)
                    break
    save_strips(allstripsNH,filedate +'allstripsNH.mat',filedate +'allstripsNH')
    save_strips(allstripsSH,filedate +'allstripsSH.mat',filedate +'allstripsSH')
    save_strips(allstrips,filedate +'allstrips.mat',filedate + 'allstrips')
    save_strips(aurorastrips,filedate +'aurorastrips.mat', filedate +'aurorastrips')

    return aurorastrips
# %%
def testLinreg():
    testitems = pd.read_pickle(r'C:\Users\ceona\Documents\GitHub\MATS-analysis\MATS-analysis\Ceona\MatsData\26aprorb3')
    """polynomial regression test for short interval (only one orbit), to check keogram as well
    """
    airglowlim = 130
    auroralim = 150
    matrix, striplist = makeStripMatrix(testitems)
    newmatrix = matrix.copy()
    fig, axs = plt.subplots(2,1)

    for rowindex in range(airglowlim,matrix.shape[0]):
        y = matrix[rowindex,:]
        x = np.arange(0,len(y)) 
        #Polynomial regression model
        if rowindex >= auroralim:
            #weird intense spikes above auroralim, like the moon pass, will be ignored in the fit
            valid_indices = np.where(y < 400)[0]
            y_filtered = y[valid_indices]
            x_filt = np.arange(0,len(y_filtered))
            model = np.poly1d(np.polyfit(x_filt, y_filtered, 1))
        else:
            valid_indices = np.where(y < 600)[0]
            y_filtered = y[valid_indices]
            x_filt = np.arange(0,len(y_filtered))
            model = np.poly1d(np.polyfit(x_filt, y_filtered, 1))
        y_reg = model(x) #new y_points
        
        newmatrix[rowindex,:] = y-y_reg
        if rowindex == airglowlim or rowindex == auroralim+20: 
            axs[0].plot(y, label = f"Row {rowindex}")
            axs[0].plot(y_reg, label = f"y_reg, Row {rowindex}")
            axs[0].plot(y-y_reg, label = f"Row {rowindex} Diff")
            axs[0].set_ylim(-30,400)
            #scipy.io.savemat('y.mat',{'y': y, 'label':'I'}) #saves to matlabfile
            #scipy.io.savemat('y_reg.mat',{'y_reg': y_reg, 'label':'I'}) #saves to matlabfile

    #Update the strips from the new matrix
    for col, strip in enumerate(striplist):
        strip.strip = newmatrix[:,col]
    
    #scipy.io.savemat('keogramgrad.mat',{'keogramgrad': newmatrix, 'label':'pixel'}) #saves to matlabfile
    #scipy.io.savemat('keogram.mat',{'keogram': matrix, 'label':'pixel'}) #saves to matlabfile
    save_strips(striplist,'test' +'allstrips.mat', 'test' +'allstrips')
    
    axs[1].pcolormesh(newmatrix,vmin=-10, vmax=300)
    fig.legend()
    #print(newmatrix[:,0] -striplist[0].strip )
    plt.show()
    return

# %%
def Main():
    start_time = DT.datetime(2023,4,26,00,0,0)
    stop_time = DT.datetime(2023,4,28,00,0,0)
    numdays = stop_time-start_time
    filedate = 'aprtest'
    items = pd.read_pickle(r'MatsData\26to27aprIR1')
    #save_TPMLT(items,filedate)
    aurorastrips = get_strips(items,numdays,filedate) 
    get_aurora_max(aurorastrips,filedate)
    return

# %%
