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

#%%
def plot_background(p3, p4, label_suffix):
    axs[bgr_index].plot(p3, dstest.altitude.values[10:endcut], color='blue', label=f'IR3{label_suffix}')
    axs[bgr_index].plot(p4, dstest.altitude.values[10:endcut], color='purple', label=f'IR4{label_suffix}')

def flatten_and_reverse(array, flip):
    """Flatten a nested array and reverse it if flip is True."""
    flat_list = [x for xs in array for x in xs]  # Flatten the array
    return flat_list[::-1] if flip else flat_list

def process_channel(df, channel):
    """Filter, drop duplicates, and set index for the given channel."""
    filtered_df = df[df['channel'] == channel].dropna().reset_index()
    filtered_df = filtered_df.drop_duplicates(subset=['EXPDate']).set_index('EXPDate')
    return filtered_df

def calc_distance(coord1, coord2):
    '''function to calculate distance between the two satellite positions'''
    distance = geodesic(coord1, coord2).kilometers
    return distance # km

def prepare_measurment(ch,ir3,ir4,offset,subtract=True,endcut=-25):

    z1,p1=prepare_profile(ch,offset)
    _,p3=prepare_profile(ir3,offset=0)
    _,p4=prepare_profile(ir4,offset=0)

    p3 = p3[10:endcut]
    p4 = p4[10:endcut]

    p3=p3-p3[-4:].mean()/1.05
    p4=p4-p4[-4:].mean()/1.05

    if subtract:

        p1=p1[10:endcut]-(p3+p4)/2

    return np.array(z1),np.array(p1),np.array(p3),np.array(p4)

def prepare_profile(ch,offset):
    # This function averages some columns and
    # calculates tangent heights
    
    image = image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2) - offset
    #cs = col_heights(ch, col, 10, spline=True) # TEST BETTER 
    cs = col_heights(ch, col, spline=True)
    heights = np.array(cs(range(ch['NROW'])))

    # multiply with factor to get right values (depends on version?)
    #profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e12)
    profile = np.array(image[:, col]*1e12)


    # set heights
    common_heights = np.arange(60,110,0.25)*1000
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile

#%%

# NLC
NLC=False
DAYGLOW=True
NIGHTGLOW=False

# CHANNEL
chan_str='IR1' # DO NOT CHANGE

OSIRIS_stray=True
#bgr_removal=False # FOR BOTH OSIRIS AND MATS
plt_advanced=True # plot
plt_bgr_channels=False

if NLC:
    data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/Odin/OSIRIS_NLCscans.mat')
    savestr='NLC'
# AIRGLOW
else:    
    data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/Odin/sm_ABand_2022_2023_v3.mat')
    if DAYGLOW:
        savestr='DAYGLOW'
        txt_fpath=f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/{chan_str}/DAYGLOW/'
    if NIGHTGLOW:
        savestr='NIGHTGLOW'
        txt_fpath=f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/{chan_str}/NIGHTGLOW/'

# %%
# sattime lat lon
if NLC:
    lat = np.array(data['lat']).flatten()
    lon = np.array(data['lon']).flatten()
    timeO = np.array(data['t']).flatten()
else:
    lat = np.array(data['sm']['lat90'])
    lon = np.array(data['sm']['lon90'])
    timeO = np.array(data['sm']['time90'])
    time = np.array(timeO)

jd = timeO + 2400000.5
t = Time(jd, format='jd')
timeO = t.to_datetime()

# SZA at 90
sza = np.zeros(len(timeO))
for i in range(len(sza)):
    sza[i] = 90 - pysolar.solar.get_altitude(
        lat[i],
        lon[i],
        pd.Timestamp(timeO[i], tzinfo=timezone.utc).to_pydatetime()
    )

#%%
#  GENERATE MEANS FROM THE ISOLATED DATA
# OSIRIS NEEDS TO BE PUT ON A COMMON GRID
#OSheights = np.arange(61.25,111.25,0.25) # SHIFT
OSheights = np.arange(60,110,0.25)

# messy OSIRIS stuff
# Initialize OSIRIS output arrays
n_mjd = len(data['sm']['mjd'])
OSinter = np.zeros((n_mjd, len(OSheights)))
OSinter3 = np.zeros((n_mjd, len(OSheights)))
OSinter4 = np.zeros((n_mjd, len(OSheights)))

# Process each time step
for i in range(n_mjd):
    zt = data['sm']['zt'][i]
    flip = False  # For sorting

    # Flatten zt and check if reversal is needed
    zt = flatten_and_reverse(zt, flip)
    if zt[0] > zt[1]:
        flip = True
        zt = zt[::-1]

    # SELECT CHANNEL FROM ODIN
    Li = data['sm'][f'L{1 if chan_str == "IR1" else 2}'][i]
    Li = flatten_and_reverse(Li, flip)

    # BACKGROUND CHANNELS
    Li3 = flatten_and_reverse(data['sm']['L3'][i], flip)
    Li4 = flatten_and_reverse(data['sm']['L4'][i], flip)

    # Interpolate to OSheights
    OSinter[i, :] = np.interp(OSheights, zt, Li)
    OSinter3[i, :] = np.interp(OSheights, zt, Li3)
    OSinter4[i, :] = np.interp(OSheights, zt, Li4)

    # Apply stray correction if needed
    if OSIRIS_stray:
        OSinter3[i, :] -= OSinter3[i, -4:].mean()
        OSinter4[i, :] -= OSinter4[i, -4:].mean()

# SUBTRACT BACKGROUND
#OSheights = OSheights - 1.25 # SHIFT

ds = xr.Dataset(
    data_vars=dict(
        lat=(["time"], lat),
        lon=(["time"], lon),
        sza=(["time"], sza),
        Li=(["time", "altitude"], OSinter),
        L3=(["time", "altitude"], OSinter3),
        L4=(["time", "altitude"], OSinter4),

    ),
    coords=dict(
        time=timeO,
        altitude=OSheights
    ),
    attrs=dict(description="ODIN"),
)

#%%
# Apply conditions for NIGHTGLOW and DAYGLOW
if NIGHTGLOW:
    ds = ds.where(ds.sza > 98)
if DAYGLOW:
    ds = ds.where(ds.sza < 90)

# Select data within specified date range and drop NaNs
ds = ds.sel(time=slice("2022-02-01", "2023-04-01")).dropna(dim='time')
starttime=datetime(2023,1,1,0,0)
stoptime=datetime(2023,2,1,0,0)
l1b_version="0.6"
download=False
run=True

dftop_full=pd.read_pickle(f'/media/waves/AVAGO/data/MATS/MATS_ODIN_DATA/{savestr}/{str(starttime)}_{str(stoptime)}.pkl')
dftop_full=dftop_full.drop(labels=['SID','schedule_description_short'],axis=1)
dftop_full['EXPDate'] = pd.to_datetime(dftop_full['EXPDate'])
#%%
# Creates a new file
file_path=f'{txt_fpath}{starttime.strftime("%Y-%m-%d")}_{stoptime.strftime("%Y-%m-%d")}.txt'
with open(file_path, 'w') as fp:
    pass

# read in BGR channels
ccdsel0, ccdsel1 = 2, 3

# Process IR3 and IR4 channels
ir2 = process_channel(dftop_full, 'IR2')
ir3 = process_channel(dftop_full, 'IR3')
ir4 = process_channel(dftop_full, 'IR4')

# Select channel for analysis
dftop_all = process_channel(dftop_full, chan_str)

#%%

##    dftop_all = dftop_all.set_index('EXPDate')
# RUN DAY-BY-DAY
polar=True
dtime =15 # minutes before and after Odin MEAS
mindis=300
endcut=-25

#CCDs = CCDs[(CCDs['TPlocaltime'].dt.time > midday)]
for dayi in range(1,30): # for every day 
    MATS_times, ODIN_times = [], []
    
    day0 = starttime + DT.timedelta(dayi)
    day1 = day0 + DT.timedelta(1) 
    ds_day=ds.sel(time=slice(day0,day1))
    dftop=dftop_all.loc[str(day0)[0:-9]:str(day1)[0:-9]]

    if (len(dftop.index) > 0) and (len(ds_day.time) > 0):

        # convert to common (epoch) time
        test_EXPdate = np.zeros(len(dftop.index.values))
        test_odint = np.zeros(len(ds_day.time))
        for i in range(0,len(dftop.index.values)):
            test_EXPdate[i]=int(dftop.index.values[i])
        for i in range(0,len(ds_day.time)):
            test_odint[i]=int(ds_day.time.values[i])

        # save str versions of the time also
        test_ttot=test_EXPdate
        test_EXPdate_str = []
        for i in range(0,len(dftop.index.values)):
            test_EXPdate_str.append(dftop.index[i])

        # lat, lon
        TPlati = dftop.TPlat.values
        TPloni = dftop.TPlon.values

        print(test_EXPdate_str)

        # MATS xarray
        MATSint=[]
        MATSint = xr.Dataset(
            data_vars=dict(
                TPlati=(["time"], TPlati),
                TPlatilon=(["time"], TPloni),
                time_str=(["time"], test_EXPdate_str),
                #sza=(["time"], sza)
            ),
            coords=dict(
                time=test_ttot,
            ),
            attrs=dict(description="MATS"),
        )

        # GET RID OF DUBLICATED TIMES (instrument anomalies?)
        MATSint = MATSint.drop_duplicates(dim = 'time')

        # ISOLATE TEN MINUTES BEFORE AND AFTER OSIRIS
        time0=((ds_day.time[0] - np.timedelta64(dtime, 'm')).values)
        time1=((ds_day.time[0] + np.timedelta64(dtime, 'm')).values)
        MATSi=MATSint.sel(time=slice(int(time0), int(time1)),drop=True)
        for i in range(0, len(test_odint)):
            time0=((ds_day.time[i] - np.timedelta64(dtime, 'm')).values)
            time1=((ds_day.time[i] + np.timedelta64(dtime, 'm')).values)
            MATSi=xr.concat([MATSi,MATSint.sel(time=slice(int(time0), int(time1)),drop=True)],dim='time')
        MATSi = MATSi.drop_duplicates(dim = 'time').dropna(dim='time') # perhaps not necessary

        if len(MATSi.time) > 0: # if any near measurements in time

            # check distance
            distance_bool = np.zeros([len(ds_day.time),len(MATSi.time)]) # matrix to flag near
            for i in range(0,len(ds_day.time)): # for every odin point
                
                #isolate odin
                coord_ODIN = (ds_day.lat.values[i], ds_day.lon.values[i])

                #isolate MATS (nearby in time)
                time0=((ds_day.time[i] - np.timedelta64(dtime, 'm')).values)
                time1=((ds_day.time[i] + np.timedelta64(dtime, 'm')).values)
                MATS_dt=MATSi.where(MATSi.time > int(time0))
                MATS_dt=MATS_dt.where(MATS_dt.time < int(time1))

                for j in range(0,len(MATS_dt.time)): # for each of these points
                    if MATS_dt.TPlati[j].notnull():
                        coord_MATS = (MATS_dt.TPlati.values[j],MATS_dt.TPlatilon.values[j])
                        dis=calc_distance(coord_ODIN, coord_MATS) # calculate distance
                        if dis < mindis: # if below limit
                            distance_bool[i,j] = 1

            MATS_meas = MATSi.where(np.sum(distance_bool,axis=0) > 0) # within any 
            ODIN_meas = ds_day.time.where(np.sum(distance_bool,axis=1) > 0) # within any
            ODIN_meas = ds_day.where(ds_day.time == ODIN_meas).dropna(dim='time')
            # note that distance_bool has structure for ds_day and MATS_meas 

            # if any points close in time and distance
            if (len(ODIN_meas.dropna('time').lat) > 0) and (len(MATS_meas.dropna('time').time) > 0):

                # plot
                if NLC or polar:
                    
                    if DAYGLOW or NLC:
                        fig, axs = pplt.subplots(figwidth='20cm',ncols=1, nrows=1,abc='a.',sharex=0,proj='splaea')
                        axs.format(coast=True, boundinglat=-60,landzorder=5)
                    else:
                        fig, axs = pplt.subplots(figwidth='20cm',ncols=1, nrows=1,abc='a.',sharex=0,proj='nplaea')
                        axs.format(coast=True, boundinglat=60,landzorder=5)
                # plot
                else:
                    pplt.rc.axesfacecolor = 'gray4'
                    fig, axs = pplt.subplots(figwidth='15cm',ncols=1, nrows=1,abc='a.',sharex=0,proj='cyl')
                    axs.format(coast=True,landzorder=5)

                # plot
                axs.scatter(ODIN_meas.dropna(dim='time').lon,ODIN_meas.dropna(dim='time').lat, s=15, c='blue',marker='x')

                # plot distances
                for j in range(0,len(ds_day.time)): # for every odin point
                    first=True
                    append_ODIN=False
                    if np.sum(distance_bool[j,:],axis=0) > 0: # continue if close to any MATS
                        append_ODIN=True 
                        for i in range(0,len(MATSi.time)):  # for every MATS point
                            if distance_bool[j,i] == 1: # close to odin point


                                MATS_times.append(MATSi.time_str.values[i]) # save this MATS point
                                
                                # plot distance between
                                axs.plot([MATSi.TPlatilon[i], ds_day.lon[j]],
                                            [MATSi.TPlati[i],ds_day.lat[j]],
                                            c='grey6', alpha=0.6, linestyle='--',
                                            linewidth=0.5)

                                # plot MATS point
                                m=axs.scatter(MATSi.TPlatilon[i], MATSi.TPlati[i],
                                            c='red',
                                            s=5, vmin=-4, vmax=4,cmap='coral')
                                
                                if first: # if first time considering this ODIN point

                                    # calculate time difference
                                    timediff= np.around((MATSi.time.values[i] - int(ds_day.time.values[j]))/(1e9*60),
                                                        decimals=2)
                                    axs.text(s=f'      dt: {str(timediff)} min',
                                                x=MATSi.TPlatilon.values[i],
                                                y=MATSi.TPlati.values[i],
                                                transform='map',fontsize=6,
                                                c='black', weight="bold",bbox=False,
                                                bboxalpha=0.1,bboxcolor='red')
                                    
                                    # plot arrow where dt was calculated
                                    axs.plot([MATSi.TPlatilon[i], ds_day.lon[j]],
                                                [MATSi.TPlati[i], ds_day.lat[j]], c='black',
                                                alpha=0.75, linestyle='--',linewidth=0.75)
                                    first=False

                        if append_ODIN:
                            ODIN_times.append(ds_day.time.values[j])

                    # plot texts for ODIN points
                    for i in range(0,len(ODIN_meas.dropna(dim='time').time)):
                        axs.text(s='      '+str(ODIN_meas.dropna(dim='time').time.values[i])[0:-10],
                                    x=ODIN_meas.dropna(dim='time').lon.values[i],
                                    y=ODIN_meas.dropna(dim='time').lat.values[i],transform='map',fontsize=6,
                                    c='green', weight="bold",bbox=False,
                                    bboxalpha=0.1,bboxcolor='green')

                fig.format(suptitle="ODIN / MATS: "+day0.strftime('%Y-%m-%d')+" -- "+day1.strftime('%Y-%m-%d'))
                axs.format(title='red: MATS - blue: Odin //')

                if not NLC:
                    if DAYGLOW:
                        directory=f"/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/{chan_str}/DAYGLOW/{day0.strftime('%Y-%m-%d')}"
                    if NIGHTGLOW:
                        directory=f"/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/{chan_str}/NIGHTGLOW/{day0.strftime('%Y-%m-%d')}"

                    if not os.path.exists(directory):
                        # Create the directory
                        os.makedirs(directory)
                        print("Directory created successfully!")

                    fig.savefig(f'{directory}/map_mindis{mindis}_mintime{dtime}.png',format='png')

                else:   
                    print('SAVING IMAGE NOT YET IMPLEMENTED FOR NLC')
                    #fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/NLC/map_mindis{mindis}_mintime{dtime}_{str(day0)}.png',format='png')

            #  GENERATE MEANS FROM THE ISOLATED DATA
            # OSIRIS NEEDS TO BE PUT ON A COMMON GRID

            dstest=ds.sel(time=ODIN_times) # xarray

            # compute means encounter by encounter (MATS)
            pass_no=0
            print(f"({len(ODIN_times)} encounters on {day0.strftime('%Y-%m-%d')})")
            #with open(file_path, 'a') as file:
            #    file.write(f"--{str(day0)}: ({len(ODIN_times)}) encounters--\n")
            for j in range(0,len(distance_bool[:,0])):

                MATS_times_tmp=MATSi.time_str.where(distance_bool[j,:] > 0).dropna(dim='time')
                MATS_times_int=MATSi.time.where(distance_bool[j,:] > 0).dropna(dim='time')
                if len(MATS_times_tmp) > 0:
                    pass_no=pass_no+1
                    # mean of all MATS measurements done at MATS_times (epoch)
                    MATS_df=dftop_all.loc[MATS_times_tmp] # pandas
                    MATS_df['EXPDate'] = MATS_df.index


                    #MALi = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 
                    MALi = np.zeros([len(MATS_df), len(common_heights),3]) # if no bgr rem
                    Ma2 = np.zeros([len(MATS_df), len(common_heights),3]) # if no bgr rem
                    Ma3 = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 
                    Ma4 = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 

                    for i in range(0,len(MATS_df)):
                        #heights, profile = prepare_profile(MATS_df.iloc[i])

                        # nearest time
                        ir2_p=ir2.iloc[ir2.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
                        ir3_p=ir3.iloc[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
                        ir4_p=ir4.iloc[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')

                        # set EXPDate accordingly
                        ir2_p['EXPDate'] = ir2.index[ir2.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                        ir3_p['EXPDate'] = ir3.index[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                        ir4_p['EXPDate'] = ir4.index[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                        
                        k = 0
                        for offs in [-15,0,15]:
                            zsi,MALi[i,:,k],Ma3[i,:,k],Ma4[i,:,k] = prepare_measurment(MATS_df.iloc[i], ir3_p, ir4_p,subtract=False,endcut=endcut,offset=offs)
                            zsi2,Ma2[i,:,k],Ma3[i,:,k],Ma4[i,:,k] = prepare_measurment(ir2_p, ir3_p, ir4_p,subtract=False,endcut=endcut,offset=offs)
                            k=k+1

                    if plt_advanced:
                        fig, axs = pplt.subplots(figwidth='25cm',ncols=3, nrows=3,abc='a.',sharey=0, sharex=0)
                        k_sza=2
                    else:
                        fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=1,abc='a.',sharey=0)
                        k_sza=1
                    
                    fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))
                    
                    labels=['MATS L','MATS M','MATS R']
                    colos=['red','orange','green']
                    for k in range(0,3):
                        for i in range(0,len(MALi[:,0,k])):
                            if i == 0:
                                axs[0].plot(MALi[i,:,k][10:endcut], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                                            color=colos[k])
                                #axs[0].plot(MALi[i,:,k], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                                #            color=colos[k])
                            else:
                                axs[0].plot(MALi[i,:,k][10:endcut], zsi[10:endcut]/1000, color=colos[k], alpha=0.2)
                                #axs[0].plot(MALi[i,:,k], zsi[10:endcut]/1000, color=colos[k], alpha=0.2)
                    axs[0].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
                    axs[0].legend(fontsize=1,ncol=2)
                    
                    if DAYGLOW:
                        xlims_main=[5e13,1.4e15]
                    if NIGHTGLOW:
                        xlims_main=[1e12,1.7e14]

                    axs[0].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                                    xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)

                    axs[k_sza].format(title=str(dstest.time[pass_no-1].values))
                    axs[k_sza].plot(MATS_df.EXPDate,MATS_df.TPsza.values, color='red')
                    
                    axs[k_sza].axhline(dstest.sza[pass_no-1], color='blue',linestyle='--',alpha=0.5)
                    axs[k_sza].axvline(dstest.time[pass_no-1].values, color='blue',linestyle='--',alpha=0.5)
                    axs[k_sza].scatter(dstest.time[pass_no-1].values,dstest.sza[pass_no-1], s=15, marker='x', color='blue')


                    timediff_first= np.around((MATS_times_int.values[0] - int(dstest.time[pass_no-1].values))/(1e9*60),
                                         decimals=2)
                    
                    timediff_last= np.around((MATS_times_int.values[-1] - int(dstest.time[pass_no-1].values))/(1e9*60),
                                         decimals=2)
                    
                    #with open(file_path, 'a') as file:
                    #    file.write(f"{str(pass_no)}/{len(ODIN_times)} (ODIN time {str(dstest.time[pass_no-1].values)}): {len(MATS_times_tmp)} profiles " + 
                    #                f"[t1: {timediff_first} min, tn: {timediff_last} min] \n")

                    # Creates a new file
                    fpath_odin=f'{directory}/pass{pass_no}_times_ODIN.txt'
                    with open(fpath_odin, 'a') as fp:
                        fp.write(str(dstest.time[pass_no-1].values))

                    # Creates a new file
                    fpath_mats=f'{directory}/pass{pass_no}_times_MATS.txt'
                    with open(fpath_mats, 'a') as fp:
                        for mtime in MATS_df.EXPDate:
                            fp.write(f'{str(mtime)} \n')

                    if not plt_advanced:
                        fig.savefig(f'{directory}/profiles_mindis{mindis}_mintime{dtime}_pass{pass_no}.png',format='png')

                    #### BACKGROUND CHANNELS

                    if plt_advanced:
                        bgr_index = 3 if plt_bgr_channels else 6
                        p3 = dstest.L3.values[pass_no-1, :]
                        p4 = dstest.L4.values[pass_no-1, :]
                        if OSIRIS_stray:
                            p3, p4 = p3[10:endcut], p4[10:endcut]
                            plot_background(p3, p4, 'o')
                        else:
                            plot_background(p3, p4, '')

                        for k in range(3):
                            for i in range(len(MALi[:, 0, k])):
                                alpha_val = 0.2
                                label = labels[k] if i == 0 else None
                                axs[1].plot(dstest.L1.values[pass_no-1, 10:endcut] / MALi[i, :, k][10:endcut], zsi[10:endcut] / 1000, alpha=alpha_val, label=label, color=colos[k])
                                if not plt_bgr_channels:
                                    axs[3].plot(Ma2[i, :, k][10:endcut], zsi[10:endcut] / 1000, color='pink', alpha=alpha_val)
                                    axs[4].plot(dstest.L2.values[pass_no-1, 10:endcut] / Ma2[i, :, k][10:endcut], zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                                if k == 0:
                                    axs[bgr_index].plot(Ma3[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='red', label='IR3m')
                                    axs[bgr_index].plot(Ma4[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='green', label='IR4m')
                                    if OSIRIS_stray:
                                        axs[bgr_index + 1].plot(p3 / Ma3[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='red', label='IR3')
                                        axs[bgr_index + 1].plot(p4 / Ma4[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='green', label='IR4')
                                if plt_bgr_channels:
                                    osi = dstest.L1.values[pass_no-1, 10:endcut] - (p3 + p4) / 2
                                    matsi = MALi[i, :, k][10:endcut] - (Ma3[i, :, k] + Ma4[i, :, k]) / 2
                                    axs[6].plot(osi, zsi[10:endcut] / 1000, alpha=alpha_val, color='blue')
                                    axs[6].plot(matsi, zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                                    axs[7].plot(osi / matsi, zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                                    axs[7].format(title='OSIRIS/MATS bgr sub', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])

                        if not plt_bgr_channels:
                            axs[3].format(title='Limb radiance (IR2)', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=xlims_main)
                            axs[3].format(xscale='log')

                        axs[bgr_index].format(title='Limb radiance (bgr)', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=[1e12, 3e14])
                        axs[bgr_index].format(xscale='log')
                        axs[1].format(title='OSIRIS/MATS', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                        axs[4].format(title='OSIRIS/MATS', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                        axs[7].format(title='OSIRIS/MATS bgr sub', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                        if plt_bgr_channels:
                            axs[6].format(title=f'Limb radiance ({chan_str})  bgr sub', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=xlims_main)
                            axs[6].format(xscale='log')
                        axs[0].format(xscale='log')
                        axs[bgr_index].format(title='Limb radiance (bgr)', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=[1e12, 4e14])
                        axs[bgr_index].format(xscale='log')

                        fig.savefig(f'{directory}/profiles_mindis{mindis}_mintime{dtime}_pass{pass_no}.png', format='png')

run=False

 # %%
