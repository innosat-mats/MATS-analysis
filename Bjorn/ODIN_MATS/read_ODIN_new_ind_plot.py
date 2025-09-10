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
import copy

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

def prepare_measurment(ch,ir3,ir4,offset,subtract=True, subtract_dark=True, endcut=-25):

    z1,p1,dark=prepare_profile(ch,offset)
    _,p3,_=prepare_profile(ir3,offset=0)
    _,p4,_=prepare_profile(ir4,offset=0)

    p3 = p3[10:endcut]
    p4 = p4[10:endcut]

    p3=p3-p3[-4:].mean()#/1.05
    p4=p4-p4[-4:].mean()#/1.05

    if subtract:

        p1=p1[10:endcut]-(p3+p4)/2

    return np.array(z1),np.array(p1),np.array(p3),np.array(p4),dark

def prepare_profile(ch,offset,subtract_dark=True):
    # This function averages some columns and
    # calculates tangent heights
    
    image = image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2) - offset
    #cs = col_heights(ch, col, 10, spline=True) # TEST BETTER 
    cs = col_heights(ch, col, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    # multiply with factor to get right values (depends on version?)
    profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e12)
    
    dark=0

    if NIGHTGLOW:
        if ch.channel == 'IR1':
            dark = profile[-25:-15].mean()
            profile = profile - dark
            #print(dark)

        if ch.channel == 'IR2':
            dark = profile[-25:-15].mean()

            #plt.plot(profile.T,heights,color='black')
            #plt.plot(profile[-25:-15].T,heights[-25:-15],color='red')
            #plt.show()

            profile = profile - dark


    #profile = np.array(image[:, col]*1e12)


    # set heights
    common_heights = np.arange(60,100,0.25)*1000
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile, dark

#%%

# NLC
NLC=False
DAYGLOW=False
NIGHTGLOW=True

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
#OSheights = np.arange(59.25,109.25,0.25) # SHIFT
OSheights = np.arange(60,100,0.25)


# messy OSIRIS stuff
# Initialize OSIRIS output arrays
n_mjd = len(data['sm']['mjd'])
OSinter = np.zeros((n_mjd, len(OSheights)))
OSinter2 = np.zeros((n_mjd, len(OSheights)))
OSinter3 = np.zeros((n_mjd, len(OSheights)))
OSinter4 = np.zeros((n_mjd, len(OSheights)))
dz = np.zeros(n_mjd)
zref = np.zeros(n_mjd)
exptime = np.zeros(n_mjd)
zmax = np.zeros(n_mjd)

# Process each time step
for i in range(n_mjd):
    zt = data['sm']['zt'][i]
    if type(data['sm']['exposure_time'][i]) is list:
        exptime[i] = exptime[i] 
    else: 
        exptime[i] = data['sm']['exposure_time'][i]
    flip = False  # For sorting

    # Flatten zt and check if reversal is needed
    zt = flatten_and_reverse(zt, flip)
    if zt[0] > zt[1]:
        flip = True
        zt = zt[::-1]

    # approx dz at 90 (sometimes)
    dz[i] = zt[-15] - zt[-16]
    zref[i] = zt[-15]

    # SELECT CHANNEL FROM ODIN
    Li1 = flatten_and_reverse(data['sm']['L1'][i], flip)
    Li2 = flatten_and_reverse(data['sm']['L2'][i], flip)

    # remove mean of top altitudes if NIGHTGLOW (straylight):
    if NIGHTGLOW:
        #print(np.mean(Li2[-3:])/10**12)
        #Li1 = Li1 - np.mean(Li1[-3:])
        #Li2 = Li2 - np.mean(Li2[-3:])
        Li1 = np.array(Li1) - 0.65*1e13
        Li2 = np.array(Li2) - 0.65*1e13

    zmax[i] = zt[-1]        

    # BACKGROUND CHANNELS
    Li3 = flatten_and_reverse(data['sm']['L3'][i], flip)
    Li4 = flatten_and_reverse(data['sm']['L4'][i], flip)

    # Interpolate to OSheights
    OSinter[i, :] = np.interp(OSheights, zt, Li1)
    OSinter2[i, :] = np.interp(OSheights, zt, Li2)
    OSinter3[i, :] = np.interp(OSheights, zt, Li3)
    OSinter4[i, :] = np.interp(OSheights, zt, Li4)

    # Apply stray correction if needed
    if OSIRIS_stray:
        OSinter3[i, :] -= OSinter3[i, -4:].mean()
        OSinter4[i, :] -= OSinter4[i, -4:].mean()

# SUBTRACT BACKGROUND
#OSheights = OSheights + 0.75 # SHIFT
common_heights = OSheights*1000

ds_sza = xr.Dataset(
    data_vars=dict(
        dz = (["time"], dz),
        zref = (["time"], zref),
        zmaxs = (["time"], zmax),
        exptime = (["time"], exptime),
        lat=(["time"], lat),
        lon=(["time"], lon),
        sza=(["time"], sza),
        L1=(["time", "altitude"], OSinter),
        L2=(["time", "altitude"], OSinter2),
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
    ds = ds_sza.where(ds_sza.sza > 100)
if DAYGLOW:
    ds = ds_sza.where(ds_sza.sza < 90)
#%%

# Select data within specified date range and drop NaNs
ds = ds.sel(time=slice("2022-02-01", "2023-05-01")).dropna(dim='time')
starttime=datetime(2023,1,1,0,0)
stoptime=datetime(2023,2,1,0,0)
#l1b_version="0.6"
download=False
run=True

# with updated data version v09
dftop_full=pd.read_pickle(f'/media/waves/AVAGO3/data/MATS/MATS_ODIN_DATA/{savestr}/v09/{str(starttime)}_{str(stoptime)}.pkl')
dftop_full=dftop_full.drop(labels=['SID','schedule_description_short'],axis=1)
dftop_full['EXPDate'] = pd.to_datetime(dftop_full['EXPDate'])

if NIGHTGLOW:
    dftop_additional = pd.read_pickle('/media/waves/AVAGO3/data/MATS/MATS_ODIN_DATA/NIGHTGLOW/v09/2022-12-20 00:00:00_2023-01-01 00:00:00.pkl')
    dftop_additional = dftop_additional.drop(labels=['SID','schedule_description_short'],axis=1)
    dftop_additional['EXPDate'] = pd.to_datetime(dftop_additional['EXPDate'])
    dftop_full = pd.concat([dftop_full, dftop_additional])
    # sort by EXPDate
    dftop_full = dftop_full.sort_values(by='EXPDate')

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

# save epoch times
test_EXPdate = np.zeros(len(dftop_all.index.values))

for i in range(0, len(dftop_all.index.values)):
    test_EXPdate[i] = int(dftop_all.index.values[i])

# save string times
test_EXPdate_str = []
for i in range(0, len(dftop_all.index.values)):
    test_EXPdate_str.append(dftop_all.index[i])

dftop_xr = xr.Dataset(
    data_vars=dict(
        TPlati=(["time"], dftop_all.TPlat.values),
        TPlatilon=(["time"], dftop_all.TPlon.values),
        time_str=(["time"], test_EXPdate_str),
    ),
    coords=dict(
        time=test_EXPdate,
    ),
    attrs=dict(description="MATS"))

dftop_xr = dftop_xr.drop_duplicates(dim='time')


#%%

##    dftop_all = dftop_all.set_index('EXPDate')
# RUN DAY-BY-DAY
polar=True
dtime = 30 # minutes before and after Odin MEAS
mindis= 150
endcut=-25

smallest = mindis
# Initialize lists to collect data
ODIN_meas_list = []
MATS_meas_list = []
distance_bool_list = []
ds_day_list = []
MATSi_list = []
days = []
nearest_list = []

#for starttime in [datetime(2022,12,20,0,0), datetime(2023,1,1,0,0)]

if NIGHTGLOW:
    startnum = -1
else:
    startnum = 0

for dayi in range(startnum, 31):  # for every day

    # isolate data of the day
    MATS_times, ODIN_times = [], []

    if dayi == -1:
        day0 = datetime(2022,12,24,0,0)
        day1 = datetime(2022,12,24,0,0) + DT.timedelta(1)
    else:
        day0 = starttime + DT.timedelta(dayi)
        day1 = day0 + DT.timedelta(1)

    ds_day = ds.sel(time=slice(day0, day1))
    #dftop = dftop_all.loc[str(day0)[0:-9]:str(day1)[0:-9]]

    MATSint = dftop_xr.where(dftop_xr.time > int(day0.timestamp() * 1e9))
    MATSint = MATSint.where(MATSint.time < int(day1.timestamp() * 1e9)).dropna(dim='time')

    # if there is data
    if (len(MATSint.time) > 0) and (len(ds_day.time) > 0):

        # save epoch times
        test_odint = np.zeros(len(ds_day.time))
        for i in range(0, len(ds_day.time)):
            test_odint[i] = int(ds_day.time.values[i])

        # for keeping track of conjunctions
        distance_bool = np.zeros([len(ds_day.time), len(MATSint.time)])

        # for every ODIN data this day
        nearest = []
        for i in range(0, len(test_odint)):

            # isolate MATS within dt
            time0 = ((ds_day.time[i] - np.timedelta64(dtime, 'm')).values)
            time1 = ((ds_day.time[i] + np.timedelta64(dtime, 'm')).values)
            MATS_dt = MATSint.where(MATSint.time > int(time0))
            MATS_dt = MATS_dt.where(MATS_dt.time < int(time1))

            coord_ODIN = (ds_day.lat.values[i], ds_day.lon.values[i])

            # check if distances are within mindis
            smallest = mindis
            for j in range(0, len(MATSint.time)):
                if MATS_dt.TPlati[j].notnull():
                    coord_MATS = (MATS_dt.TPlati.values[j], MATS_dt.TPlatilon.values[j])
                    dis = calc_distance(coord_ODIN, coord_MATS)
                    
                    if dis < mindis:
                        distance_bool[i, j] = 1
                        if dis < smallest:
                            smallest = dis
                            nearest_mats = j # not in use yet
            
            if smallest < mindis:
                nearest.append(nearest_mats)
                        
            
        MATS_meas = MATSint.where(np.sum(distance_bool, axis=0) > 0)
        ODIN_meas = ds_day.time.where(np.sum(distance_bool, axis=1) > 0)
        ODIN_meas = ds_day.where(ds_day.time == ODIN_meas).dropna(dim='time')
        
        nearest_list.append(nearest)
        ODIN_meas_list.append(ODIN_meas)
        MATS_meas_list.append(MATS_meas)
        distance_bool_list.append(distance_bool)
        days.append(day0)

        print(f"({len(ODIN_meas.time)} encounters on {day0.strftime('%Y-%m-%d')})")

#%%
# LOOP TO PLOT THE MAP OF ENCOUNTERS

# technically same as day loop
for ddd in range(0,len(distance_bool_list)):
    
    ODIN_meas = ODIN_meas_list[ddd]
    distance_bool = distance_bool_list[ddd]
    MATS_meas = MATS_meas_list[ddd]
    day0 = days[ddd]


    #day0 = starttime + DT.timedelta(dayi)
    day1 = day0 + DT.timedelta(1)

    ds_day = ds.sel(time=slice(day0, day1))
    MATSint = dftop_xr.where(dftop_xr.time > int(day0.timestamp() * 1e9))
    MATSint = MATSint.where(MATSint.time < int(day1.timestamp() * 1e9)).dropna(dim='time')

    if len(ODIN_meas.dropna(dim='time').time) > 0:
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
                for i in range(0,len(MATSint.time)):  # for every MATS point
                    if distance_bool[j,i] == 1: # close to odin point

                        MATS_times.append(MATSint.time_str.values[i]) # save this MATS point
                        
                        # plot distance between
                        axs.plot([MATSint.TPlatilon[i], ds_day.lon[j]],
                                    [MATSint.TPlati[i],ds_day.lat[j]],
                                    c='grey6', alpha=0.6, linestyle='--',
                                    linewidth=0.5)

                        # plot MATS point
                        m=axs.scatter(MATSint.TPlatilon[i], MATSint.TPlati[i],
                                    c='red',
                                    s=5, vmin=-4, vmax=4,cmap='coral')
                        
                        if first: # if first time considering this ODIN point

                            # calculate time difference
                            timediff= np.around((MATSint.time.values[i] - int(ds_day.time.values[j]))/(1e9*60),
                                                decimals=2)
                            axs.text(s=f'      dt: {str(timediff)} min',
                                        x=MATSint.TPlatilon.values[i],
                                        y=MATSint.TPlati.values[i],
                                        transform='map',fontsize=6,
                                        c='black', weight="bold",bbox=False,
                                        bboxalpha=0.1,bboxcolor='red')
                            
                            # plot arrow where dt was calculated
                            axs.plot([MATSint.TPlatilon[i], ds_day.lon[j]],
                                        [MATSint.TPlati[i], ds_day.lat[j]], c='black',
                                        alpha=0.75, linestyle='--',linewidth=0.75)
                            first=False

                #if append_ODIN:
                #    ODIN_times.append(ds_day.time.values[j])

            # plot texts for ODIN points
            for i in range(0,len(ODIN_meas.dropna(dim='time').time)):
                axs.text(s='      '+str(ODIN_meas.dropna(dim='time').time.values[i])[0:-10],
                            x=ODIN_meas.dropna(dim='time').lon.values[i],
                            y=ODIN_meas.dropna(dim='time').lat.values[i],transform='map',fontsize=6,
                            c='green', weight="bold",bbox=False,
                            bboxalpha=0.1,bboxcolor='green')

        fig.format(suptitle="ODIN / MATS: "+day0.strftime('%Y-%m-%d')+" -- "+day1.strftime('%Y-%m-%d'))
        axs.format(title='red: MATS - blue: Odin //')

#%%
"""
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
"""
#%%%

# Initialize a dictionary to store results for each day and pass
results_dict = {}

# Loop through each day (`ddd`)
for ddd in range(len(distance_bool_list)):

    ODIN_meas = ODIN_meas_list[ddd]
    distance_bool = distance_bool_list[ddd]
    MATS_meas = MATS_meas_list[ddd]
    day0 = days[ddd]
    nearest = nearest_list[ddd]
    day1 = day0 + DT.timedelta(1)

    #ds_day = ds.sel(time=slice(day0, day1))
    MATSint = dftop_xr.where(dftop_xr.time > int(day0.timestamp() * 1e9))
    MATSint = MATSint.where(MATSint.time < int(day1.timestamp() * 1e9)).dropna(dim='time')

    #  GENERATE MEANS FROM THE ISOLATED DATA
    # OSIRIS NEEDS TO BE PUT ON A COMMON GRID
    dstest=ds.sel(time=ODIN_meas.time) # xarray 

    # Dictionary to store Ma fields and ds for each day (ddd)
    daily_data = {
        'Ma1': [],
        'Ma2': [],
        'Ma3': [],
        'Ma4': [],
        'Ma1_sza': [],
        'Ma2_sza': [],
        'Ma2_coords': [],
        'MaExps': [],
        'days': [],
        'ds': [],  # Store the daily ds for reference
        'dark1': [],
        'dark2': []
    }

    # compute means encounter by encounter
    pass_no=0

    print(f"({len(ODIN_meas.time)} encounters on {day0.strftime('%Y-%m-%d')})")
    with open(file_path, 'a') as file:
        file.write(f"--{str(day0)}: ({len(ODIN_times)}) encounters--\n")
    for j in range(0,len(distance_bool[:,0])):

        MATS_times_tmp=MATS_meas.time_str.where(distance_bool[j,:] > 0).dropna(dim='time')

        if len(MATS_times_tmp) > 0:
            pass_no=pass_no+1

            # mean of all MATS measurements done at MATS_times (epoch)
            MATS_df=dftop_all.loc[MATS_times_tmp] # pandas
            MATS_df['EXPDate'] = MATS_df.index

            Ma1 = np.zeros([len(MATS_df), len(common_heights),3])
            Ma2 = np.zeros([len(MATS_df), len(common_heights),3])
            Ma3 = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 
            Ma4 = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 
            Ma1_sza = np.zeros(len(MATS_df))
            Ma2_sza = np.zeros(len(MATS_df))
            Ma2_coords = np.zeros([len(MATS_df),2])
            MaExps = []
            dark1 = np.zeros([len(MATS_df),3])
            dark2 = np.zeros([len(MATS_df),3])
            

            for i in range(0,len(MATS_df)):

                # nearest time
                ir2_p=ir2.iloc[ir2.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
                ir3_p=ir3.iloc[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
                ir4_p=ir4.iloc[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
                #print(ir2.index[ir2.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')])

                # set EXPDate accordingly
                ir2_p['EXPDate'] = ir2.index[ir2.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                ir3_p['EXPDate'] = ir3.index[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                ir4_p['EXPDate'] = ir4.index[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                
                Ma1_sza[i] = MATS_df.TPsza.values[i]
                Ma2_sza[i] = ir2_p.TPsza
                Ma2_coords[i] = [ir2_p.TPlat, ir2_p.TPlon]
                MaExps.append(MATS_df.EXPDate.values[i])
                #print(MaExps[i])

                k = 0
                for offs in [-15,0,15]:
                    zsi, Ma1[i, :, k], Ma3[i, :, k], Ma4[i, :, k],dark1[i,k] = prepare_measurment(MATS_df.iloc[i], ir3_p, ir4_p, subtract=False, endcut=endcut, offset=offs)
                    zsi2,Ma2[i,:,k], _, _,dark2[i,k] = prepare_measurment(ir2_p, ir3_p, ir4_p,subtract=False,endcut=endcut,offset=offs)
                    
                    # Apply rolling mean to smooth out the data
                    window_size = 6  # Adjust the window size as needed
                    Ma1[i, :, k] = pd.Series(Ma1[i, :, k]).rolling(window=window_size, min_periods=1).mean().values
                    Ma2[i, :, k] = pd.Series(Ma2[i, :, k]).rolling(window=window_size, min_periods=1).mean().values
                    Ma3[i, :, k] = pd.Series(Ma3[i, :, k]).rolling(window=window_size, min_periods=1).mean().values
                    Ma4[i, :, k] = pd.Series(Ma4[i, :, k]).rolling(window=window_size, min_periods=1).mean().values
                    
                    # apply calibration factors to convert to latest version (not if v09)
                    Ma1[i, :, k] = Ma1[i, :, k] #* (10.2/10.5)
                    Ma2[i, :, k] = Ma2[i, :, k] #* (2.99/3.36)
                    Ma3[i, :, k] = Ma3[i, :, k] #* (18.3/19.1)
                    Ma4[i, :, k] = Ma4[i, :, k] #* (25.5/25.4)
                    dark1[i,k] = dark1[i,k] #* (10.2/10.5)
                    dark2[i,k] = dark2[i,k] #* (2.99/3.36)

                    k=k+1

            # Store Ma fields and additional info in the daily_data dict
            daily_data['Ma1'].append(Ma1)
            daily_data['Ma2'].append(Ma2)
            daily_data['Ma3'].append(Ma3)
            daily_data['Ma4'].append(Ma4)
            daily_data['Ma1_sza'].append(Ma1_sza)
            daily_data['Ma2_sza'].append(Ma2_sza)
            daily_data['Ma2_coords'].append(Ma2_coords)
            daily_data['MaExps'].append(MaExps)
            daily_data['dark1'].append(dark1)
            daily_data['dark2'].append(dark2)

            # rolling mean odin
            window_size = 1 # no rolling mean
            dstest.L1.values[pass_no-1, :] = pd.Series(dstest.L1.values[pass_no-1, :]).rolling(window=window_size, min_periods=1).mean().values
            dstest.L2.values[pass_no-1, :] = pd.Series(dstest.L2.values[pass_no-1, :]).rolling(window=window_size, min_periods=1).mean().values

        daily_data['ds'].append(dstest)

    daily_data['days'].append(day0)

    # Store all daily data in the results_dict
    results_dict[ddd] = daily_data


#%%

for results in results_dict.values(): # for every day
                
    day0 = results['days'][0]
    dstest = results['ds'][0]


    if len(results['Ma1']) > 0: # if there are any encounters

        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])

            if plt_advanced:
                fig, axs = pplt.subplots(figwidth='25cm',ncols=3, nrows=3,abc='a.',sharey=0, sharex=0)
                k_sza=2
                k_sza2=5
            else:
                fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=1,abc='a.',sharey=0)
                k_sza=1

            fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

            labels=['MATS L','MATS M','MATS R']
            colos=['red','orange','green']
            for k in range(0,3):
                for i in range(0,len(Ma1[:,0,k])):
                    if i == 0:
                        axs[0].plot(Ma1[i,:,k][10:endcut], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                                    color=colos[k])
                        #axs[0].plot(Ma1[i,:,k], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                        #            color=colos[k])
                    else:
                        axs[0].plot(Ma1[i,:,k][10:endcut], zsi[10:endcut]/1000, color=colos[k], alpha=0.2)
                        #axs[0].plot(Ma1[i,:,k], zsi[10:endcut]/1000, color=colos[k], alpha=0.2)
            axs[0].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            axs[0].legend(fontsize=1,ncol=2)

            if DAYGLOW:
                xlims_main=[5e13,1.4e15]
            if NIGHTGLOW:
                xlims_main=[1e12,1.7e14]

            axs[0].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            axs[0].format(xscale='log')

            axs[k_sza].format(title=str(dstest.time[pass_no-1].values))
            axs[k_sza].plot(MaExps,Ma1_sza, color='red')

            axs[k_sza].axhline(dstest.sza[pass_no-1], color='blue',linestyle='--',alpha=0.5)
            axs[k_sza].axvline(dstest.time[pass_no-1].values, color='blue',linestyle='--',alpha=0.5)
            axs[k_sza].scatter(dstest.time[pass_no-1].values,dstest.sza[pass_no-1], s=15, marker='x', color='blue')


            #with open(file_path, 'a') as file:
            #    file.write(f"{str(pass_no)}/{len(ODIN_times)} (ODIN time {str(dstest.time[pass_no-1].values)}): {len(MATS_times_tmp)} profiles " + 
            #                f"[t1: {timediff_first} min, tn: {timediff_last} min] \n")

            # Creates a new file
            #fpath_odin=f'{directory}/pass{pass_no}_times_ODIN.txt'
            #with open(fpath_odin, 'a') as fp:
            #    fp.write(str(dstest.time[pass_no-1].values))

            # Creates a new file
            #fpath_mats=f'{directory}/pass{pass_no}_times_MATS.txt'
            #with open(fpath_mats, 'a') as fp:
            #    for mtime in MATS_df.EXPDate:
            #        fp.write(f'{str(mtime)} \n')

            #if not plt_advanced:
            #    fig.savefig(f'{directory}/profiles_mindis{mindis}_mintime{dtime}_pass{pass_no}.png',format='png')

            #### BACKGROUND CHANNELS

            if plt_advanced:
                bgr_index = 3 if plt_bgr_channels else 6
                p3 = dstest.L3.values[pass_no-1, :]
                p4 = dstest.L4.values[pass_no-1, :]
                
                if DAYGLOW:
                    if OSIRIS_stray:
                        p3, p4 = p3[10:endcut], p4[10:endcut]
                        plot_background(p3, p4, 'o')
                    else:
                        plot_background(p3, p4, '')

                for k in range(0,3):
                    for i in range(len(Ma1[:, 0, k])):
                        alpha_val = 0.2
                        label = labels[k] if i == 0 else None
                        axs[1].plot(dstest.L1.values[pass_no-1, 10:endcut] / Ma1[i, :, k][10:endcut], zsi[10:endcut] / 1000, alpha=alpha_val, label=label, color=colos[k])
                        if not plt_bgr_channels:
                            axs[3].plot(Ma2[i, :, k][10:endcut], zsi[10:endcut] / 1000, color=colos[k], alpha=alpha_val)
                            axs[3].plot(dstest.L2.values[pass_no-1, 10:endcut], zsi[10:endcut] / 1000, color='blue')
                            axs[4].plot(dstest.L2.values[pass_no-1, 10:endcut] / Ma2[i, :, k][10:endcut], zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                        if k == 0:
                            if DAYGLOW:
                                axs[bgr_index].plot(Ma3[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='red', label='IR3m')
                                axs[bgr_index].plot(Ma4[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='green', label='IR4m')
                                if OSIRIS_stray:
                                    axs[bgr_index + 1].plot(p3 / Ma3[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='red', label='IR3')
                                    axs[bgr_index + 1].plot(p4 / Ma4[i, :, k], zsi[10:endcut] / 1000, alpha=alpha_val, color='green', label='IR4')
                        if DAYGLOW:
                            if plt_bgr_channels:
                                osi = dstest.L1.values[pass_no-1, 10:endcut] - (p3 + p4) / 2
                                matsi = Ma1[i, :, k][10:endcut] - (Ma3[i, :, k] + Ma4[i, :, k]) / 2
                                axs[6].plot(osi, zsi[10:endcut] / 1000, alpha=alpha_val, color='blue')
                                axs[6].plot(matsi, zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                                axs[7].plot(osi / matsi, zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                                axs[7].format(title='OSIRIS/MATS bgr sub', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                
                        if NIGHTGLOW:
                            mats_t = Ma1[i, :, k][10:endcut]/Ma2[i, :, k][10:endcut]
                            osiris_t = dstest.L1.values[pass_no-1][10:endcut]/dstest.L2.values[pass_no-1][10:endcut]
                            axs[6].plot(mats_t, zsi[10:endcut] / 1000, alpha=alpha_val, color='blue')
                            axs[6].plot(osiris_t, zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                            axs[7].plot(osiris_t / mats_t, zsi[10:endcut] / 1000, alpha=alpha_val, color=colos[k])
                            axs[7].format(title='OSIRIS/MATS bgr sub', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                
                
                axs[k_sza2].format(title=str(dstest.time[pass_no-1].values))
                axs[k_sza2].scatter(MaExps,Ma2_sza, color='red')
                
                if not plt_bgr_channels:
                    axs[3].format(title='Limb radiance (IR2)', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=xlims_main)
                    axs[3].format(xscale='log')

                axs[1].format(title='OSIRIS/MATS', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                axs[4].format(title='OSIRIS/MATS', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                axs[7].format(title='OSIRIS/MATS', ylabel='Tangent altitude [km]', xlabel='Fraction', xlim=[0.25, 1.5])
                if DAYGLOW:
                    axs[bgr_index].format(title='Limb radiance (bgr)', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=[1e12, 3e14])
                    axs[bgr_index].format(xscale='log')
                    if plt_bgr_channels:
                        axs[6].format(title=f'Limb radiance ({chan_str})  bgr sub', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=xlims_main)
                        axs[6].format(xscale='log')
                        axs[bgr_index].format(title='Limb radiance (bgr)', ylabel='Tangent altitude [km]', xlabel='[m-2 s-1 str-1 nm-1]', xlim=[1e12, 4e14])
                        axs[bgr_index].format(xscale='log')
                if NIGHTGLOW:
                    axs[bgr_index].format(title='IR1/IR2', ylabel='Tangent altitude [km]', xlim=[0.25, 1.75])
                    
            


            #fig.savefig(f'{directory}/profiles_mindis{mindis}_mintime{dtime}_pass{pass_no}.png', format='png')


 # %%

for results in results_dict.values(): # for every day
                
    day0 = results['days'][0]
    dstest = results['ds'][0]

    if len(results['Ma1']) > 0: # if there are any encounters

        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
        

            #Ma1 = Ma1 - np.mean(Ma1[:, -10::, :], axis=1)



            fig, axs = pplt.subplots(figwidth='16cm',ncols=2, nrows=1,abc='a.',sharey=1, sharex=0,
                                     includepanels=True)
            pxs=axs.panel('r', space=0, width='7em')

            fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

            labels=['MATS L','MATS M','MATS R']
            colos=['red','orange','green']
            for k in range(0,3):
                for i in range(0,len(Ma1[:,0,k])):
                    if i == 0:
                        
                        # IR1 
                        #axs[0].plot(Ma1[i, :, k][10:endcut], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                        #            color=colos[k])

                        # IR2
                        #axs[1].plot(Ma2[i, :, k][10:endcut], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                        #            color=colos[k])
                        pass
                        #pxs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut]/Ma2[i, :, k][10:endcut], zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
                        #            color=colos[k])
                                    
                    else:
                        # IR1 
                        #axs[0].plot(Ma1[i, :, k][10:endcut], zsi[10:endcut]/1000, alpha=0.1, color='red')
                        #pxs[0].plot(dstest.L1.values[pass_no-1,:][10:endcut]/Ma1[i, :, k][10:endcut], zsi[10:endcut]/1000, 
                        #            alpha=0.2,color=colos[k])
                        
                        # IR2
                        #axs[1].plot(Ma2[i, :, k][10:endcut], zsi[10:endcut]/1000, alpha=0.1,
                        #            color='red')
                        pass
                        #pxs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut]/Ma2[i, :, k][10:endcut], zsi[10:endcut]/1000, alpha=0.2,
                        #            color=colos[k])
                
                meanz = np.mean(Ma1[:, :, k],axis=0)
                axs[0].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                    color=colos[k],linewidth=0.75)
                pxs[0].plot(dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                            color=colos[k],linewidth=0.75)
                
                meanz = np.mean(Ma2[:, :, k],axis=0)
                 
                axs[1].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                    color=colos[k],linewidth=0.75)
                pxs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                            color=colos[k],linewidth=0.75)

            fractions_IR1 = np.zeros([len(Ma1[:,0,0]),len(Ma1[0,:,0]),3])
            fractions_IR2 = np.zeros([len(Ma2[:,0,0]),len(Ma2[0,:,0]),3])
            for k in range(0,3):
                for i in range(0,len(Ma1[:,0,k])):
                    fractions_IR1[i,:,k] = dstest.L1.values[pass_no-1,:]/Ma1[i,:,k]
                    fractions_IR2[i,:,k] = dstest.L2.values[pass_no-1,:]/Ma2[i,:,k]

            # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
            axs[0].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            axs[1].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma2[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

            axs[0].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            axs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

            # fill pxs between min and max fractions
            pxs[0].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
            pxs[1].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[0].legend(fontsize=1,ncol=2)
            #axs[3].legend(fontsize=1,ncol=2)

            if DAYGLOW:
                xlims_main=[5e13,1.4e15]
            if NIGHTGLOW:
                xlims_main=[3e12,1.3e14]

            axs[0].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            axs[1].format(title=f'Limb radiance (IR2)',ylabel='Tangent altitude [km]',
                         xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            pxs[0].format(title='OSIRIS/MATS',ylabel='Tangent altitude [km]',xlabel='Fraction',xlim=[0.4, 1.6])
            pxs[1].format(title='OSIRIS/MATS',ylabel='Tangent altitude [km]',xlabel='Fraction',xlim=[0.5, 1.5])
            
            axs[0].format(xscale='log')
            axs[1].format(xscale='log')
            #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
            
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)

            directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
            if NIGHTGLOW:
                savstr = f'NG_{day0.strftime("%Y-%m-%d")}'
            else:
                savstr = f'DG_{day0.strftime("%Y-%m-%d")}'
            fig.savefig(f'{directory}/{savstr}_profiles_mindis{mindis}_mintime{dtime}_pass{pass_no}_v09.png',format='png')

#%%
# DARK CURRENT REMOVAL SETTINGS

if NIGHTGLOW:
    zind0_ir1 = 0
    zind1_ir1 = 0
    print(f'zind0_ir1: {zsi[zind0_ir1]}, zind1_ir1: {zsi[zind1_ir1]}')

else: 
    zind0_ir1 = 0
    zind1_ir1 = 0
    print(f'zind0_ir1: {zsi[zind0_ir1]}, zind1_ir1: {zsi[zind1_ir1]}')

endcut=-1
#%%
import matplotlib.pyplot as plt
# plot all mats l1 profiles to see how highest values look
for results in results_dict.values(): # 

    dstest = results['ds'][0]
    for i in range(0,len(results['Ma1'])):
        for j in range(0,len(dstest.time)): # for every pass
            Ma1 = np.array(results['Ma1'][j])
        plt.plot(Ma1[i,:,1]/1e13, zsi/1000, alpha=0.2)
    plt.title('MATS L1 profiles')
#%%
# plot all osiris l1 profiles to see hwo highest values look
for results in results_dict.values(): # for every day
    dstest = results['ds'][0]
    for i in range(0,len(dstest.time)):
        plt.plot(dstest.L1.values[i,:], dstest.altitude.values, alpha=0.2)
    plt.title('OSIRIS L1 profiles')

#%%
# plot all mats l1 profiles to see how highest values look
for results in results_dict.values(): # 

    dstest = results['ds'][0]
    for i in range(0,len(results['Ma2'])):
        for j in range(0,len(dstest.time)): # for every pass
            Ma2 = np.array(results['Ma2'][j])
        plt.plot(Ma2[i,:,1]/1e13, zsi/1000, alpha=0.2)
    plt.title('MATS L2 profiles')
#%%
# plot all osiris l2 profiles to see hwo highest values look
for results in results_dict.values(): # for every day
    dstest = results['ds'][0]
    for i in range(0,len(dstest.time)):
        plt.plot(dstest.L2.values[i,:], dstest.altitude.values, alpha=0.2)
    plt.title('OSIRIS L2 profiles')
# %%

# statistics

diffs_ir1 = np.zeros([0,len(zsi[10:endcut])])
diffs_ir2 = np.zeros([0,len(zsi[10:endcut])])

refs_ir1 = np.zeros([0,len(zsi[10:endcut])])
refs_ir2 = np.zeros([0,len(zsi[10:endcut])])

tot_pass = 0

import copy

for results in results_dict.values(): # for every day

    results_temp = copy.deepcopy(results)            
    day0 = results_temp['days'][0]
    dstest = results_temp['ds'][0]

    if NIGHTGLOW:
        dstest.L1.values = dstest.L1.values/10**13
        dstest.L2.values = dstest.L2.values/10**13
    else:
        dstest.L1.values = dstest.L1.values/10**14
        dstest.L2.values = dstest.L2.values/10**14

    if len(results['Ma1']) > 0: # if there are any encounters


        diff_ir1 = np.zeros([0,len(zsi[10:endcut])])
        diff_ir2 = np.zeros([0,len(zsi[10:endcut])])

        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1


            run=True

            # remove mean from highest altitudes
            #dstest.L1.values[pass_no-1,:] = dstest.L1.values[pass_no-1,:] - np.mean(dstest.L1.values[pass_no-1,-30:endcut])
            #dstest.L2.values[pass_no-1,:] = dstest.L2.values[pass_no-1,:] - np.mean(dstest.L2.values[pass_no-1,-30:endcut])


            if NIGHTGLOW:
                Ma1 = np.array(results['Ma1'][i])/10**13
                Ma2 = (np.array(results['Ma2'][i]))/10**13 
            else:
                Ma1 = np.array(results['Ma1'][i])/10**14
                Ma2 = (np.array(results['Ma2'][i]))/10**14

            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)

            if NIGHTGLOW:
                if MaExps[0].day == 3:
                    run = False
                
                if (MaExps[0].day == 24 and pass_no == 1):
                    print('skip')
                    run = False

                #if MaExps[0].day == 11:
                #    run = False

                if (MaExps[0].day == 9) and pass_no == 2:
                    run = False

                if MaExps[0].day == 15 and pass_no == 2:
                    run = False

            #else:
            #    if (MaExps[0].day == 11 and pass_no == 2):
            #        run = False


            if run:

                #diff_ir1 = np.zeros([3,len(Ma1[:,0,0]),len(zsi[10:endcut])])

                #k=2
                for k in range(0,3):
                    for j in range(0,len(Ma1[:,0,k])):

                        # statistics to compare dstest.L1 and Ma1
                        # differences

                        # remove mean from highest altitudes
                        if NIGHTGLOW:
                            Ma1[j,:,k] = Ma1[j,:,k]  #- np.mean(Ma1[j,zind0_ir1:zind1_ir1,k],axis=0)
                            Ma2[j,:,k] = Ma2[j,:,k]  #- np.mean(Ma2[j,zind0_ir1:zind1_ir1,k],axis=0)

                        dv_ir1 = (Ma1[j,:,k][10:endcut] - dstest.L1.values[pass_no-1,:][10:endcut])
                        diff_ir1 = np.append(diff_ir1, [dv_ir1],axis=0)

                        dv_ir2 = (Ma2[j,:,k][10:endcut] - dstest.L2.values[pass_no-1,:][10:endcut])
                        diff_ir2 = np.append(diff_ir2, [dv_ir2],axis=0)
                        
                        #diff_ir1 = np.append(diff_ir1, dstest.L1.values[pass_no-1,:][10:endcut] - Ma1[j,:,k][10:endcut],axis=0)
                        #diff_ir1= dstest.L1.values[pass_no-1,:][10:endcut] - Ma1[j,:,k][10:endcut]

                    refs_ir1 = np.append(refs_ir1, [dstest.L1.values[pass_no-1,:][10:endcut]],axis=0)
                    refs_ir2 = np.append(refs_ir2, [dstest.L2.values[pass_no-1,:][10:endcut]],axis=0)

                    # append all diffs_ir1 to one array
                    diffs_ir1 = np.append(diffs_ir1,diff_ir1,axis=0)
                    diffs_ir2 = np.append(diffs_ir2,diff_ir2,axis=0)

                tot_pass = tot_pass + 1
                #print(dstest.time[pass_no-1].values)


mean_diff = np.mean(diffs_ir1,axis=0)
std_diff = np.std(diffs_ir1,axis=0)
#std_diff_norm = np.std(diffs_ir1,axis=0)
median_diff = np.median(diffs_ir1,axis=0)
max_diff = np.max(diffs_ir1,axis=0)
min_diff = np.min(diffs_ir1,axis=0)

mean_diff2 = np.mean(diffs_ir2,axis=0)
std_diff2 = np.std(diffs_ir2,axis=0)
median_diff2 = np.median(diffs_ir2,axis=0)
max_diff2 = np.max(diffs_ir2,axis=0)
min_diff2 = np.min(diffs_ir2,axis=0)

# means of ref
mean_ref = np.mean(refs_ir1,axis=0)
mean_ref2 = np.mean(refs_ir2,axis=0)

# plot
#if i == 0:
#    axs[0].plot(diff_ir1, zsi[10:endcut]/1000, alpha=0.2,label=labels[k],
#                color=colos[k])
#else:
#    axs[0].plot(diff_ir1, zsi[10:endcut]/1000, alpha=0.2,color=colos[k])
fig, axs = pplt.subplots(figwidth='16cm',ncols=2, nrows=1,abc='a.',sharey=1, sharex=0,
            includepanels=True)
pxs=axs.panel('r', space=0, width='7em')
axs[0].format(suptitle=f'Mean difference between OSIRIS and MATS')

# plot statistics
axs[0].plot(mean_diff, zsi[10:endcut]/1000, alpha=1,label='IR1',
            color='red',linewidth=1)
#axs[0].plot(median_diff, zsi[10:endcut]/1000, alpha=1,label='median',
#            color='red',linewidth=0.75,linestyle='--')

axs[0].plot(mean_diff2, zsi[10:endcut]/1000, alpha=1,label='IR2',
            color='blue',linewidth=1)

pxs[0].plot(mean_diff/mean_ref, zsi[10:endcut]/1000, color='red',linewidth=1)
pxs[0].plot(mean_diff2/mean_ref2, zsi[10:endcut]/1000, color='blue',linewidth=1)



# fill between min and max
axs[0].fill_betweenx(zsi[10:endcut]/1000, mean_diff - 1*std_diff, mean_diff + 1*std_diff, alpha=0.1, color='red')
axs[0].fill_betweenx(zsi[10:endcut]/1000, mean_diff2 - 1*std_diff2, mean_diff2 + 1*std_diff2, alpha=0.1, color='blue')

# fill between diff+std/ref and diff-std/ref
#pxs[0].fill_betweenx(zsi[10:endcut]/1000, (mean_diff - 1*std_diff)/mean_ref, (mean_diff + 1*std_diff)/mean_ref, alpha=0.1, color='red')
#pxs[0].fill_betweenx(zsi[10:endcut]/1000, (mean_diff2 - 1*std_diff2)/mean_ref2, (mean_diff2 + 1*std_diff2)/mean_ref2, alpha=0.1, color='blue')

std_diff = std_diff/np.sqrt(tot_pass)
std_diff2 = std_diff2/np.sqrt(tot_pass)

# fill between min and max
axs[0].fill_betweenx(zsi[10:endcut]/1000, mean_diff - 1*std_diff, mean_diff + 1*std_diff, alpha=0.3, color='red')
axs[0].fill_betweenx(zsi[10:endcut]/1000, mean_diff2 - 1*std_diff2, mean_diff2 + 1*std_diff2, alpha=0.3, color='blue')

## fill between diff+std/ref and diff-std/ref
pxs[0].fill_betweenx(zsi[10:endcut]/1000, (mean_diff - 1*std_diff)/mean_ref, (mean_diff + 1*std_diff)/mean_ref, alpha=0.3, color='red')
pxs[0].fill_betweenx(zsi[10:endcut]/1000, (mean_diff2 - 1*std_diff2)/mean_ref2, (mean_diff2 + 1*std_diff2)/mean_ref2, alpha=0.3, color='blue')








#axs[0].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
#axs[0].format(xscale='log')
#axs[0].legend(fontsize=1,ncol=1)
axs[1].legend(fontsize=1,ncol=1)
pxs[0].format(title='OSIRIS/MATS',ylabel='Tangent altitude [km]',xlim=[-0.6, 0.6])
pxs[0].axvline(0, color='black', linestyle='--', alpha=1, linewidth=0.75)

axs[0].axvline(0, color='black', linestyle='--', alpha=1, linewidth=0.75)



xlims_main = [-2.1, 1.8]

savstr = 'Nightglow'
if DAYGLOW:
    savstr = 'Dayglow'
    xlims_main = [-1.2, 1.9]
    axs[0].format(title=f'{savstr}',ylabel='Tangent altitude [km]',
                    xlabel=r'$10^{14}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main)
    axs[1].format(title=f'{savstr}',ylabel='Tangent altitude [km]',
                    xlabel=r'$10^{14}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main)

    #axs[0].legend(fontsize=1,ncol=1, loc='upper left')
else:
    #axs[0].legend(fontsize=0.5,ncol=1, loc='lower left')
    axs[0].format(title=f'{savstr}',ylabel='Tangent altitude [km]',
                    xlabel=r'$10^{13}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main)
    axs[1].format(title=f'{savstr}',ylabel='Tangent altitude [km]',
                    xlabel=r'$10^{14}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main)

axs[0].legend(fontsize=1,ncol=1, loc='upper left')
axs[0].format(title='Nightglow')
axs[1].format(title='Dayglow')


directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/statistics'
fig.savefig(f'{directory}/{savstr}_mindis{mindis}_dtime{dtime}_totpass{tot_pass}_{len(diffs_ir1[:,0])}profiles_v09.png', format='png')

# %%


#%% PAPER SCRIPT (2 EXAMPLES AND A PASS)

array = [  # the "picture" (0 == nothing, 1 == subplot A, 2 == subplot B, etc.)
    [1, 1, 2, 2,3,3, 3,3,3],
    [1, 1, 2,2,3,3, 3,3,3],
    [4, 4, 5,5,3,3, 3,3,3],
    [4, 4, 5, 5,3,3, 3,3,3],
]



fig, axs = pplt.subplots(array,figwidth='16cm',abc='a.',sharey=1, sharex=0,
                            includepanels=True,figheight='10cm', proj=(None, None, 'ortho', None, None),
                            proj_kw={'central_latitude': 82,'central_longitude': 130})
#pxs=axs.panel('r', space=0, width='7em')

#remove panel from pxs[2]
#pxs[2].remove()


#increase width of axs[2] plot
axs[2].format(coast=True, reso='med', land=True,
              landcolor='white',landzorder=2,facecolor='lightsteelblue',
              lonlim=(70, 150), latlim=(70, 89), latlines=10, loninline=True, lonlabels=True)

# zoom in in axs[2]
#axs[2].set_boundary(radius=0.4)

# add gridlines to axs[2]
axs[2].gridlines(alpha=0.6,zorder=3)

axs[2].text(s=r' $85^\circ$N',x=80,y=85,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $80^\circ$N',x=80,y=80,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $75^\circ$N',x=80,y=75,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $70^\circ$N',x=80,y=70,transform='map',fontsize=7, c='red')

# longitude lines
axs[2].text(s=r' $120^\circ$E',x=120,y=76,transform='map',fontsize=7, c='red')
axs[2].text(s=r'     $60^\circ$E',x=61,y=75,transform='map',fontsize=7, c='red')


# make axs[2] a polar plot
#axs[2].format(boundinglat=-60, coast=True)

fig.format(toplabels=('        2023-01-11 03:13:16', '          2023-01-15 02:09:05', None))

plot=False
import copy

for results in results_dict.values(): # for every day

    results_temp = copy.deepcopy(results)            
    day0 = results_temp['days'][0]
    dstest = results_temp['ds'][0]

    dstest.L1.values = dstest.L1.values/10**14
    dstest.L2.values = dstest.L2.values/10**14

    if len(results['Ma1']) > 0: # if there are any encounters

        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])/10**14
            Ma2 = np.array(results['Ma2'][i])/10**14
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])

            # convert MaExps[0] to datetime from numpy datetime64
            MaExps = pd.to_datetime(MaExps)

            # check if MaEXPs contain a given datetime day
            if MaExps[0].day == 11 and pass_no == 1:
                plot=True
                plt_ind=0
            # elseif for the other day
            elif MaExps[0].day == 15 and pass_no == 1:
                plot=True
                plt_ind=1

            if plot:
                #print(MaExps[0])
                labels=['MATS L','MATS M','MATS R']
                colos=['red','black','green']

                pxs0=axs[plt_ind+0].panel('r', space=0, width='5em')
                pxs1=axs[plt_ind+3].panel('r', space=0, width='5em')

                for k in range(0,3):
                    
                    meanz = np.mean(Ma1[:, :, k],axis=0) 
                    # remove mean from highest altitudes
                    meanz = meanz #- np.mean(meanz[zind0_ir1:zind1_ir1])
                    
                    axs[plt_ind+0].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                        color=colos[k],linewidth=0.75)
                    pxs0.plot((dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut])**-1, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                color=colos[k],linewidth=0.75)
                    
                    meanz = np.mean(Ma2[:, :, k],axis=0)
                    meanz = meanz #- np.mean(meanz[zind0_ir1:zind1_ir1])
                    axs[plt_ind+3].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                        color=colos[k],linewidth=0.75)
                    pxs1.plot((dstest.L2.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut])**-1, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                color=colos[k],linewidth=0.75)

                fractions_IR1 = np.zeros([len(Ma1[:,0,0]),len(Ma1[0,:,0]),3])
                fractions_IR2 = np.zeros([len(Ma2[:,0,0]),len(Ma2[0,:,0]),3])
                for k in range(0,3):
                    for i in range(0,len(Ma1[:,0,k])):
                        Ma1[i,:,k] = Ma1[i,:,k] #- np.mean(Ma1[i,zind0_ir1:zind1_ir1,k])
                        Ma2[i,:,k] = Ma2[i,:,k] #- np.mean(Ma2[i,zind0_ir1:zind1_ir1,k])
                        fractions_IR1[i,:,k] = (dstest.L1.values[pass_no-1,:]/(Ma1[i,:,k]))**(-1)
                        fractions_IR2[i,:,k] = (dstest.L2.values[pass_no-1,:]/(Ma2[i,:,k]))**(-1)


                # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
                axs[plt_ind+0].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
                axs[plt_ind+3].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma2[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

                axs[plt_ind+0].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
                axs[plt_ind+3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

                # plot ma2 positions on map

                # plot osisris positions on map
                osi_coords = np.array([dstest.lat[pass_no-1],dstest.lon[pass_no-1]])
                #print(osi_coords)

                if MaExps[0].day == 11:
                    axs[2].scatter(Ma2_coords[:,1],Ma2_coords[:,0],color='red',label='MATS',s=3,alpha=1, zorder=10)
                    axs[2].scatter(osi_coords[1],osi_coords[0],color='blue',label='OSIRIS',marker='x', s=12,alpha=1, zorder=10)

                else:
                    axs[2].scatter(Ma2_coords[:,1],Ma2_coords[:,0],color='red',s=3,alpha=1, zorder=10)
                    axs[2].scatter(osi_coords[1],osi_coords[0],color='blue',marker='x', s=12,alpha=1, zorder=10)



                # fill pxs between min and max fractions
                pxs0.fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
                pxs1.fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
                #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
                #axs[0].legend(fontsize=1,ncol=2)
                #axs[3].legend(fontsize=1,ncol=2)

                if DAYGLOW:
                    xlims_main=[5e13,1.4e15]
                if NIGHTGLOW:
                    xlims_main=[3e12/10**14,1.3e14/10**14]

                axs[plt_ind+0].format(title=f'{chan_str}',ylabel='Tangent altitude [km]',
                                xlabel=r'$10^{14}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main,fontsize=8)
                axs[plt_ind+3].format(title=f'IR2',ylabel='Tangent altitude [km]',
                            xlabel=r'$10^{14}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main,fontsize=8)
                pxs0.format(xlim=[0.3, 1.7], title='Ratios', fontsize=8,facecolor='gray3')
                pxs1.format(xlim=[0.3, 1.7],title='Ratios', fontsize=8,facecolor='gray3')
                pxs0.axvline(1.0, color='black', linestyle='--', alpha=0.5, linewidth=0.5)
                pxs1.axvline(1.0, color='black', linestyle='--', alpha=0.5, linewidth=0.5)

                # add texts to axs
                #axs[plt_ind+0].text(0.5, 0.9, 'IR1')
                #axs[plt_ind+0].text(0.5, 0.9, 'IR2')

                axs[plt_ind+0].format(xscale='log')
                axs[plt_ind+3].format(xscale='log')


                lgd=axs[2].legend(loc='lower right',fontsize=1)
                axs[2].text(s=f'      2023-01-{MaExps[0].day}',
                                x=osi_coords[1],
                                y=osi_coords[0],
                                transform='map',fontsize=6,
                                c='black', weight="bold",bbox=False,
                                bboxalpha=0.1,bboxcolor='red')

                #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
                #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
                
                #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)

                directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
                if NIGHTGLOW:
                    savstr = f'NG_{day0.strftime("%Y-%m-%d")}'
                else:
                    savstr = f'DG_{day0.strftime("%Y-%m-%d")}'

                plot=False

lgd.legendHandles[0]._sizes = [10]
lgd.legendHandles[1]._sizes = [10]

# change font size of legend items
for text in lgd.get_texts():
    text.set_fontsize(8)

axs[2].format(title='Tangent point locations',fontsize=8)
axs[0].set_xlabel('')
axs[1].set_xlabel('')

fig.savefig(f'{directory}/NG_EXAMPLES_v1.png',format='png')

# %%
#%% PAPER SCRIPT (2 x 2 DAYGLOW IMAGE)

array = [  # the "picture" (0 == nothing, 1 == subplot A, 2 == subplot B, etc.)
    [1, 1, 2, 2,3,3, 3,3,3],
    [1, 1, 2,2,3,3, 3,3,3],
    [4, 4, 5,5,3,3, 3,3,3],
    [4, 4, 5, 5,3,3, 3,3,3],
]


fig, axs = pplt.subplots(array,figwidth='16cm',abc='a.',sharey=1, sharex=0,
                            includepanels=True,figheight='10cm', proj=(None, None, 'ortho',None, None),
                            proj_kw={'central_latitude': -65,'central_longitude': 0}) #110
#pxs=axs.panel('r', space=0, width='7em')

#remove panel from pxs[2]
#pxs[2].remove()

fig.format(toplabels=('      2023-02-11 12:46:35', '        2023-02-21 14:03:28', None))
#increase width of axs[2] plot
axs[2].format(coast=True, reso='med', land=True, lonlim=(0, 220), 
              landcolor='white',landzorder=2,facecolor='lightsteelblue',
              latlim=(-70, -90), latlines=20, inlinelabels=True, rotatelabels=True, latlabels=True)


axs[2].text(s=r' $90^\circ$S',x=60,y=-90,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $85^\circ$S',x=60,y=-85,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $80^\circ$S',x=60,y=-80,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $75^\circ$S',x=60,y=-75,transform='map',fontsize=7, c='red')

# longitude lines
axs[2].text(s=r' $0^\circ$E',x=0,y=-79,transform='map',fontsize=7, c='red')
axs[2].text(s=r' $60^\circ$W',x=-60,y=-77,transform='map',fontsize=7, c='red')
# zoom in in axs[2]
#axs[2].set_boundary(radius=0.4)

# add gridlines to axs[2]
axs[2].gridlines(alpha=0.6,zorder=3)


# make axs[2] a polar plot
#axs[2].format(boundinglat=-60, coast=True)

plot=False
import copy

for results in results_dict.values(): # for every day

    results_temp = copy.deepcopy(results)            
    day0 = results_temp['days'][0]
    dstest = results_temp['ds'][0]

    dstest.L1.values = dstest.L1.values/10**15
    dstest.L2.values = dstest.L2.values/10**15

    #print(dstest.L1.values)

    if len(results['Ma1']) > 0: # if there are any encounters

        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])/10**15
            Ma2 = np.array(results['Ma2'][i])/10**15
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])

            # convert MaExps[0] to datetime from numpy datetime64
            MaExps = pd.to_datetime(MaExps)

            # check if MaEXPs contain a given datetime day
            if MaExps[0].day == 11 and pass_no == 2:
                plot=True
                plt_ind=0

            if MaExps[0].day == 21 and pass_no == 1:
                plot=True
                plt_ind=1

            if plot:
                #print(MaExps[0])
                labels=['MATS L','MATS M','MATS R']
                colos=['red','black','green']

                pxs0=axs[plt_ind+0].panel('r', space=0, width='5em')
                pxs1=axs[plt_ind+3].panel('r', space=0, width='5em')

                for k in range(0,3):
                    
                    meanz = np.mean(Ma1[:, :, k],axis=0)
                    axs[plt_ind+0].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                        color=colos[k],linewidth=0.75)
                    pxs0.plot((dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut])**-1, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                color=colos[k],linewidth=0.75)
                    
                    meanz = np.mean(Ma2[:, :, k],axis=0)
                    axs[plt_ind+3].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                        color=colos[k],linewidth=0.75)
                    pxs1.plot((dstest.L2.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut])**-1, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                color=colos[k],linewidth=0.75)

                fractions_IR1 = np.zeros([len(Ma1[:,0,0]),len(Ma1[0,:,0]),3])
                fractions_IR2 = np.zeros([len(Ma2[:,0,0]),len(Ma2[0,:,0]),3])
                for k in range(0,3):
                    for i in range(0,len(Ma1[:,0,k])):
                        fractions_IR1[i,:,k] = (dstest.L1.values[pass_no-1,:]/Ma1[i,:,k])**(-1)
                        fractions_IR2[i,:,k] = (dstest.L2.values[pass_no-1,:]/Ma2[i,:,k])**(-1)

                # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
                axs[plt_ind+0].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
                axs[plt_ind+3].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma2[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

                axs[plt_ind+0].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
                axs[plt_ind+3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

                # plot ma2 positions on map

                # plot osisris positions on map
                osi_coords = np.array([dstest.lat[pass_no-1],dstest.lon[pass_no-1]])
                #print(osi_coords)

                if pass_no == 1:
                    axs[2].scatter(Ma2_coords[:,1],Ma2_coords[:,0],color='red',label='MATS',s=3,alpha=1,zorder=3)
                    axs[2].scatter(osi_coords[1],osi_coords[0],color='blue',label='OSIRIS',marker='x', s=12,alpha=1,zorder=3)

                else:
                    axs[2].scatter(Ma2_coords[:,1],Ma2_coords[:,0],color='red',s=3,alpha=1,zorder=3)
                    axs[2].scatter(osi_coords[1],osi_coords[0],color='blue',marker='x', s=12,alpha=1,zorder=3)



                # fill pxs between min and max fractions
                pxs0.fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
                pxs1.fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
                #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
                #axs[0].legend(fontsize=1,ncol=2)
                #axs[3].legend(fontsize=1,ncol=2)

                if DAYGLOW:
                    xlims_main=[5e13/10**15,1.4e15/10**15]
                if NIGHTGLOW:
                    xlims_main=[3e12,1.3e14]

                axs[plt_ind+0].format(title=f'{chan_str}',ylabel='Tangent altitude [km]',
                                xlabel=r'$10^{15}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main,fontsize=8)
                axs[plt_ind+3].format(title=f'IR2',ylabel='Tangent altitude [km]',
                            xlabel=r'$10^{15}$ $m^{-2} s^{-1} str^{-1} nm^{-1}$',xlim=xlims_main, fontsize=8)
                pxs0.format(xlim=[0.3, 1.7], title='Ratios', fontsize=8,facecolor='gray3')
                pxs1.format(xlim=[0.3, 1.7],title='Ratios', fontsize=8,facecolor='gray3')
                pxs0.axvline(1.0, color='black', linestyle='--', alpha=0.5, linewidth=0.5)
                pxs1.axvline(1.0, color='black', linestyle='--', alpha=0.5, linewidth=0.5)

                # add texts to axs
                #axs[plt_ind+0].text(0.5, 0.9, 'IR1')
                #axs[plt_ind+0].text(0.5, 0.9, 'IR2')

                axs[plt_ind+0].format(xscale='log')
                axs[plt_ind+3].format(xscale='log')


                lgd=axs[2].legend(loc='upper right',fontsize=1)

                axs[2].text(s=f'      2023-02-{MaExps[0].day}',
                                x=osi_coords[1],
                                y=osi_coords[0],
                                transform='map',fontsize=6,
                                c='black', weight="bold",bbox=False,
                                bboxalpha=0.1,bboxcolor='red')

                #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
                #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
                
                #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
                #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)

                directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
                if NIGHTGLOW:
                    savstr = f'NG_{day0.strftime("%Y-%m-%d")}'
                else:
                    savstr = f'DG_{day0.strftime("%Y-%m-%d")}'

                plot=False

lgd.legendHandles[0]._sizes = [10]
lgd.legendHandles[1]._sizes = [10]

# change font size of legend items
for text in lgd.get_texts():
    text.set_fontsize(8)

axs[2].format(title='Tangent point locations', fontsize=8)
axs[0].set_xlabel('')
axs[1].set_xlabel('')

fig.savefig(f'{directory}/DG_EXAMPLES_v1.png',format='png')
# %%

### APPENDIX FIGURE

if DAYGLOW:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=4, nrows=5,abc='a.',sharey=1, sharex=0,
                                includepanels=True)
    pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Dayglow IR1')
else:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=4, nrows=4,abc='a.',sharey=1, sharex=0,
                                includepanels=True)
    pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Nightglow IR1')

#fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

o = 0

for results in results_dict.values(): # for every day
    results_temp = copy.deepcopy(results)            
    day0 = results_temp['days'][0]
    dstest = results_temp['ds'][0]

    if NIGHTGLOW:
        dstest.L1.values = dstest.L1.values/10**14
        dstest.L2.values = dstest.L2.values/10**14
    else:
        dstest.L1.values = dstest.L1.values/10**15
        dstest.L2.values = dstest.L2.values/10**15

    if len(results['Ma1']) > 0: # if there are any encounters

        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1



            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)

            if NIGHTGLOW:
                Ma1 = Ma1/10**14
                Ma2 = Ma2/10**14
            else:
                Ma1 = Ma1/10**15
                Ma2 = Ma2/10**15

            labels=['MATS L','MATS M','MATS R']
            colos=['red','black','green']

            for k in range(0,3):
                
                meanz = np.mean(Ma1[:, :, k],axis=0)

                # remove mean of highest altitudes
                #if NIGHTGLOW:
                    #meanz = np.abs(meanz - np.mean(meanz[zind0_ir1:zind1_ir1]))

                axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                                    color=colos[k],linewidth=0.75)
                pxs[o].plot((dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut])**-1, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                            color=colos[k],linewidth=0.75)
            
                #meanz = np.mean(Ma2[:, :, k],axis=0)
                #axs[1].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #    color=colos[k],linewidth=0.75)
                #pxs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)

            fractions_IR1 = np.zeros([len(Ma1[:,0,0]),len(Ma1[0,:,0]),3])
            fractions_IR2 = np.zeros([len(Ma2[:,0,0]),len(Ma2[0,:,0]),3])
            for k in range(0,3):
                for j in range(0,len(Ma1[:,0,k])):
                    #if NIGHTGLOW:
                    #    fractions_IR1[j,:,k] = dstest.L1.values[pass_no-1,:]/(Ma1[j,:,k])
                    #    fractions_IR2[j,:,k] = dstest.L2.values[pass_no-1,:]/(Ma2[j,:,k])
                    #else:
                    fractions_IR1[j,:,k] = (dstest.L1.values[pass_no-1,:]/Ma1[j,:,k])**-1
                    fractions_IR2[j,:,k] = (dstest.L2.values[pass_no-1,:]/Ma2[j,:,k])**-1

            # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines

            axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[1].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma2[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

            axs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

            # fill pxs between min and max fractions
            pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
            #pxs[1].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[0].legend(fontsize=1,ncol=2)
            #axs[3].legend(fontsize=1,ncol=2)



            if DAYGLOW:
                xlims_main=[0.05,2]
            if NIGHTGLOW:
                xlims_main=[3e12,1.3e14]

                # inv rem
                xlims_main=[0.01,2]

            axs[o].format(xlim=xlims_main)
            #axs[1].format(title=f'Limb radiance (IR2)',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            pxs[o].format(xlim=[0.4, 1.6],title='Ratios',facecolor='gray3')
            #pxs[1].format(title='OSIRIS/MATS',ylabel='Tangent altitude [km]',xlabel='Fraction',xlim=[0.5, 1.5])
            
            axs[o].format(xscale='log')

            # time as axs title
            if NIGHTGLOW:
                axs[o].format(title=f'   {MaExps[0].day}/1 ({pass_no})',fontsize=8, ylabel='Tangent altitude [km]')

            else:
                axs[o].format(title=f'   {MaExps[0].day}/2 ({pass_no})',fontsize=8, ylabel='Tangent altitude [km]')

            #axs[1].format(xscale='log')
            #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
            
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)


            

            o = o + 1

#remove last axs

if NIGHTGLOW:
    axs[0].format(title=f'   24/12 (1)',fontsize=8, ylabel='Tangent altitude [km]')
    axs[1].format(title=f'   24/12 (2)',fontsize=8, ylabel='Tangent altitude [km]')


directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
if NIGHTGLOW:
    savstr = f'NG'
else:
    savstr = f'DG'
fig.savefig(f'{directory}/ALL_{savstr}_IR1_v2.png',format='png')
# %%
# %%

### APPENDIX FIGURE
if DAYGLOW:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=4, nrows=5,abc='a.',sharey=1, sharex=0,
                                includepanels=True)
    pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Dayglow IR2')
else:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=4, nrows=4,abc='a.',sharey=1, sharex=0,
                                includepanels=True)
    pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Nightglow IR2')

#fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

o = 0

for results in results_dict.values(): # for every day

    results_temp = copy.deepcopy(results)            
    day0 = results_temp['days'][0]
    dstest = results_temp['ds'][0]

    if NIGHTGLOW:
        dstest.L1.values = dstest.L1.values/10**14
        dstest.L2.values = dstest.L2.values/10**14
    else:
        dstest.L1.values = dstest.L1.values/10**15
        dstest.L2.values = dstest.L2.values/10**15


    if len(results['Ma1']) > 0: # if there are any encounters
        
        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            # remove mean of highest altitudes
            #dstest.L2.values[pass_no-1,:][10:endcut] = dstest.L2.values[pass_no-1,:][10:endcut] - np.mean(dstest.L2.values[pass_no-1,:][-35:endcut])

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i]) # + 1*4*10**12
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)

            if NIGHTGLOW:
                Ma1 = Ma1/10**14
                Ma2 = Ma2/10**14
            else:
                Ma1 = Ma1/10**15
                Ma2 = Ma2/10**15


            labels=['MATS L','MATS M','MATS R']
            colos=['red','black','green']
            for k in range(0,3):
                
                #meanz = np.mean(Ma1[:, :, k],axis=0)
                #axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #                    color=colos[k],linewidth=0.75)
                #pxs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)
                
                meanz = np.mean(Ma2[:, :, k],axis=0)

                if NIGHTGLOW:
                    # remove mean of highest altitudes
                    meanz = meanz

                axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                    color=colos[k],linewidth=0.75)
                pxs[o].plot((dstest.L2.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut])**-1, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                            color=colos[k],linewidth=0.75)

            fractions_IR1 = np.zeros([len(Ma1[:,0,0]),len(Ma1[0,:,0]),3])
            fractions_IR2 = np.zeros([len(Ma2[:,0,0]),len(Ma2[0,:,0]),3])
            for k in range(0,3):
                for j in range(0,len(Ma1[:,0,k])):
                    #if NIGHTGLOW:
                    #    fractions_IR1[j,:,k] = dstest.L1.values[pass_no-1,:]/(Ma1[j,:,k] - np.mean(Ma1[j,zind0_ir1:zind1_ir1,k]))
                    #    fractions_IR2[j,:,k] = dstest.L2.values[pass_no-1,:]/(Ma2[j,:,k] - np.mean(Ma2[j,zind0_ir1:zind1_ir1,k]))
                    #else:
                    fractions_IR1[j,:,k] = (dstest.L1.values[pass_no-1,:]/Ma1[j,:,k])**-1
                    fractions_IR2[j,:,k] = (dstest.L2.values[pass_no-1,:]/Ma2[j,:,k])**-1

            # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

            axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma2[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

            #axs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            axs[o].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

            # fill pxs between min and max fractions
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
            pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[0].legend(fontsize=1,ncol=2)
            #axs[3].legend(fontsize=1,ncol=2)

            if NIGHTGLOW:
                axs[o].format(title=f'   {MaExps[0].day}/1 ({pass_no})',fontsize=8, ylabel='Tangent altitude [km]')
            else:
                axs[o].format(title=f'   {MaExps[0].day}/2 ({pass_no})',fontsize=8, ylabel='Tangent altitude [km]')

            if DAYGLOW:
                xlims_main=[0.05,2]
            if NIGHTGLOW:
                xlims_main=[0.01,2]

            #axs[o].format(xlim=xlims_main)
            axs[o].format(xlim=xlims_main)
            #pxs[o].format(xlim=[0.4, 1.6])
            pxs[o].format(xlim=[0.4, 1.6],title='Ratios',facecolor='gray3')
            
            #axs[o].format(xscale='log')
            axs[o].format(xscale='log')
            #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
            
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)



            

            o = o + 1

if NIGHTGLOW:
    axs[0].format(title=f'   24/12 (1)',fontsize=8, ylabel='Tangent altitude [km]')
    axs[1].format(title=f'   24/12 (2)',fontsize=8, ylabel='Tangent altitude [km]')

directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
if NIGHTGLOW:
    savstr = f'NG'
else:
    savstr = f'DG'
fig.savefig(f'{directory}/ALL_{savstr}_IR2_v2.png',format='png')

# %%
## FRACTIONS
### APPENDIX FIGURE

if DAYGLOW:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=4, nrows=5,abc='a.',sharey=1, sharex=1,
                                includepanels=True)
    pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Dayglow IR1/IR2')
else:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=4, nrows=5,abc='a.',sharey=1, sharex=1,
                                includepanels=True)
    pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Nightglow IR1/IR2')

#fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

pxs.format(facecolor='gray3')

o = 0

for results in results_dict.values(): # for every day
                
    results_temp = copy.deepcopy(results)            
    day0 = results_temp['days'][0]
    dstest = results_temp['ds'][0]


    if len(results['Ma1']) > 0: # if there are any encounters
        
        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            # remove mean of highest altitudes
            #dstest.L1.values[pass_no-1,:][10:endcut] = dstest.L1.values[pass_no-1,:][10:endcut] - np.mean(dstest.L1.values[pass_no-1,:][-35:endcut])
            #dstest.L2.values[pass_no-1,:][10:endcut] = dstest.L2.values[pass_no-1,:][10:endcut] - np.mean(dstest.L2.values[pass_no-1,:][-35:endcut])

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)


            labels=['MATS L','MATS M','MATS R']
            colos=['red','orange','green']
            for k in range(0,3):
                
                #meanz = np.mean(Ma1[:, :, k],axis=0)
                #axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #                    color=colos[k],linewidth=0.75)
                #pxs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)
                
                meanz1 = np.mean(Ma1[:, :, k],axis=0)
                meanz2 = np.mean(Ma2[:, :, k],axis=0)

                # remove mean of highest altitudes
                #meanz1 = meanz1 - np.mean(meanz1[zind0_ir1:zind1_ir1])
                #meanz2 = meanz2 - np.mean(meanz2[zind0_ir1:zind1_ir1])

                temper=meanz1[10:endcut]/meanz2[10:endcut]

                axs[o].plot(temper, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                    color=colos[k],linewidth=0.75)
                
                axs[o].axvline(1.0, color='black', linestyle='--', alpha=0.5, linewidth=0.5)
                
                osiris_temper=dstest.L1.values[pass_no-1,:][10:endcut]/dstest.L2.values[pass_no-1,:][10:endcut]

                pxs[o].plot(osiris_temper/temper, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                            color=colos[k],linewidth=0.75)


            # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma2[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')

            #axs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            axs[o].plot(osiris_temper, dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

            # fill pxs between min and max fractions
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[0].legend(fontsize=1,ncol=2)
            #axs[3].legend(fontsize=1,ncol=2)

            axs[o].format(title=f'   {MaExps[0].day}/1 ({pass_no})',fontsize=8)

            #if DAYGLOW:
            xlims_main=[0,2.4]
            #if NIGHTGLOW:
            #    xlims_main=[3e12,1.3e14]

            axs[o].format(xlim=xlims_main)
            #axs[o].format(xlim=xlims_main)
            pxs[o].format(xlim=[0.30, 1.1])
            #pxs[o].format(xlim=[0.4, 1.6])
            
            #axs[o].format(xscale='log')
            #axs[o].format(xscale='log')
            #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
            
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)



            

            o = o + 1
directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
if NIGHTGLOW:
    savstr = f'NG'
else:
    savstr = f'DG'
fig.savefig(f'{directory}/ALL_{savstr}_IR1divIR2_v09.png',format='png')

# %%
# %%
# SCATTER PLOT
### APPENDIX FIGURE


fig, axs = pplt.subplots(figwidth='16cm',figheight='18cm', ncols=1, nrows=2,abc='a.',sharey=1, sharex=1,
                        includepanels=True)

o = 0

for results in results_dict.values(): # for every day
                
    day0 = results['days'][0]
    dstest = results['ds'][0]

    if len(results['Ma1']) > 0: # if there are any encounters
        
        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)

            labels=['MATS L','MATS M','MATS R']
            colos=['red','orange','green']
            
            for k in [0,2]:
            #k = 0
                meanz1 = np.mean(Ma1[:, :, k],axis=0)
                meanz2 = np.mean(Ma2[:, :, k],axis=0)
                for alt in range(10,len(meanz1)):
                    axs[0].scatter(meanz1[alt], dstest.L1.values[pass_no-1,alt],
                                    color='red',linewidth=0.75,s=1,alpha=0.050)
                    axs[1].scatter(meanz2[alt], dstest.L2.values[pass_no-1,alt], 
                                color='red',linewidth=0.75,s=1,alpha=0.050)
                    
                    # plot straight line
                    axs[0].plot([0,1.5e15],[0,1.5e15],color='black',linewidth=0.5)
                    axs[1].plot([0,1.5e15],[0,1.5e15],color='black',linewidth=0.5)
            #axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
            #                    color=colos[k],linewidth=0.75)

            
            #meanz = np.mean(Ma2[:, :, k],axis=0)
            #axs[1].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
            #    color=colos[k],linewidth=0.75)
            #pxs[1].plot(dstest.L2.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
            #            color=colos[k],linewidth=0.75)

axs[0].format(title=f'IR1',ylabel='OSIRIS',xlabel='MATS',fontsize=8, xlim=[0,0.5e14], ylim=[0,0.6e14])
axs[1].format(title=f'IR2',ylabel='OSIRIS',xlabel='MATS',fontsize=8, xlim=[0,0.5e14],ylim=[0,0.6e14])
#remove last axs


directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
if NIGHTGLOW:
    savstr = f'NG'
else:
    savstr = f'DG'
fig.savefig(f'{directory}/ALL_{savstr}_SCATTER_v09.png',format='png')
# %%

#### BGR CHANNELS
### APPENDIX FIGURE
if DAYGLOW:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='24cm', ncols=4, nrows=5,abc='a.',sharey=1, sharex=1,
                                includepanels=False)
    #pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Dayglow IR3')
else:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='24cm', ncols=4, nrows=4,abc='a.',sharey=1, sharex=1,
                                includepanels=False)
    #pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Nightglow IR3')

#fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

o = 0

for results in results_dict.values(): # for every day
                
    day0 = results['days'][0]
    dstest = results['ds'][0]

    if len(results['Ma1']) > 0: # if there are any encounters
        
        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)


            labels=['MATS L','MATS M','MATS R']
            colos=['red','orange','green']
            for k in range(0,3):
                
                #meanz = np.mean(Ma1[:, :, k],axis=0)
                #axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #                    color=colos[k],linewidth=0.75)
                #pxs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)
                
                meanz = np.mean(Ma3[:, :, k],axis=0)
                axs[o].plot(meanz, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                    color=colos[k],linewidth=0.75)
                #pxs[o].plot(dstest.L3.values[pass_no-1,:][10:endcut]/meanz[:], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)

            fractions_IR1 = np.zeros([len(Ma1[:,0,0]),len(Ma1[0,:,0]),3])
            fractions_IR2 = np.zeros([len(Ma2[:,0,0]),len(Ma2[0,:,0]),3])
            #for k in range(0,3):
            #    for j in range(0,len(Ma1[:,0,k])):
            #        fractions_IR1[j,:,k] = dstest.L1.values[pass_no-1,:]/Ma1[j,:,k]
            #        fractions_IR2[j,:,k] = dstest.L2.values[pass_no-1,:]/Ma2[j,:,k]

            # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma3[:,:,:],axis=(0,2)), np.nanmax(Ma3[:,:,:],axis=(0,2)), alpha=0.2, color='red')

            #axs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            axs[o].plot(dstest.L3.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

            # fill pxs between min and max fractions
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[0].legend(fontsize=1,ncol=2)
            #axs[3].legend(fontsize=1,ncol=2)

            axs[o].format(title=f'   {MaExps[0].day}/1 ({pass_no})',fontsize=8)

            if DAYGLOW:
                xlims_main=[5e13,1.4e15]
            if NIGHTGLOW:
                xlims_main=[3e12,1.3e14]

            #axs[o].format(xlim=xlims_main)
            axs[o].format(xlim=[-3*10**12,3*10**12])
            pxs[o].format(xlim=[-5, 5])
            #pxs[o].format(xlim=[0.4, 1.6])
            
            #axs[o].format(xscale='log')
            #axs[o].format(xscale='log')
            #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
            
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)



            

            o = o + 1
directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
if NIGHTGLOW:
    savstr = f'NG'
else:
    savstr = f'DG'
fig.savefig(f'{directory}/ALL_{savstr}_IR3_v09.png',format='png')

# %%
# %%

#### BGR CHANNELS

### APPENDIX FIGURE
if DAYGLOW:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='24cm', ncols=4, nrows=5,abc='a.',sharey=1, sharex=1,
                                includepanels=False)
    #pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Dayglow IR4')
else:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='24cm', ncols=4, nrows=4,abc='a.',sharey=1, sharex=1,
                                includepanels=False)
    #pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('Nightglow IR4')

#fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

o = 0

for results in results_dict.values(): # for every day
                
    day0 = results['days'][0]
    dstest = results['ds'][0]

    if len(results['Ma1']) > 0: # if there are any encounters
        
        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            Ma1 = np.array(results['Ma1'][i])
            Ma2 = np.array(results['Ma2'][i])
            Ma3 = np.array(results['Ma3'][i])
            Ma4 = np.array(results['Ma4'][i])
            Ma1_sza = np.array(results['Ma1_sza'][i])
            Ma2_sza = np.array(results['Ma2_sza'][i])
            Ma2_coords = np.array(results['Ma2_coords'][i])
            MaExps = np.array(results['MaExps'][i])
            MaExps = pd.to_datetime(MaExps)


            labels=['MATS L','MATS M','MATS R']
            colos=['red','orange','green']
            for k in range(0,3):
                
                #meanz = np.mean(Ma1[:, :, k],axis=0)
                #axs[o].plot(meanz[10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #                    color=colos[k],linewidth=0.75)
                #pxs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut]/meanz[:][10:endcut], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)
                
                meanz = np.mean(Ma4[:, :, k],axis=0)
                axs[o].plot(meanz, zsi[10:endcut]/1000, alpha=1,label=labels[k],
                    color=colos[k],linewidth=0.75)
                #pxs[o].plot(dstest.L3.values[pass_no-1,:][10:endcut]/meanz[:], zsi[10:endcut]/1000, alpha=1,label=labels[k],
                #            color=colos[k],linewidth=0.75)

            #for k in range(0,3):
            #    for j in range(0,len(Ma1[:,0,k])):
            #        fractions_IR1[j,:,k] = dstest.L1.values[pass_no-1,:]/Ma1[j,:,k]
            #        fractions_IR2[j,:,k] = dstest.L2.values[pass_no-1,:]/Ma2[j,:,k]

            # instead of plotting the Ma1[i,:,0] lines, shade the area between the lines
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma1[:,:,:],axis=(0,2))[10:endcut], np.nanmax(Ma1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(Ma3[:,:,:],axis=(0,2)), np.nanmax(Ma3[:,:,:],axis=(0,2)), alpha=0.2, color='red')

            #axs[o].plot(dstest.L1.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            axs[o].plot(dstest.L4.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')

            # fill pxs between min and max fractions
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR1,axis=(0,2))[10:endcut], np.nanmax(fractions_IR1[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')            
            #pxs[o].fill_betweenx(zsi[10:endcut]/1000, np.nanmin(fractions_IR2,axis=(0,2))[10:endcut], np.nanmax(fractions_IR2[:,:,:],axis=(0,2))[10:endcut], alpha=0.2, color='red')
            #axs[3].plot(dstest.L2.values[pass_no-1,:][10:endcut], dstest.altitude.values[10:endcut],color='blue',label='OSIRIS')
            #axs[0].legend(fontsize=1,ncol=2)
            #axs[3].legend(fontsize=1,ncol=2)

            axs[o].format(title=f'   {MaExps[0].day}/1 ({pass_no})',fontsize=8)

            if DAYGLOW:
                xlims_main=[5e13,1.4e15]
            if NIGHTGLOW:
                xlims_main=[3e12,1.3e14]

            #axs[o].format(xlim=xlims_main)
            axs[o].format(xlim=[-3*10**12,3*10**12])
            pxs[o].format(xlim=[-5, 5])
            #pxs[o].format(xlim=[0.4, 1.6])
            
            #axs[o].format(xscale='log')
            #axs[o].format(xscale='log')
            #axs[2].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #            xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #             xlabel='[m-2 s-1 str-1 nm-1]',xlim=[0.50, 1.5])
            
            #axs[3].format(title=f'Limb radiance ({chan_str})',ylabel='Tangent altitude [km]',
            #                xlabel='[m-2 s-1 str-1 nm-1]',xlim=xlims_main)



            

            o = o + 1
directory='/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/PAPER/profiles'
if NIGHTGLOW:
    savstr = f'NG'
else:
    savstr = f'DG'
fig.savefig(f'{directory}/ALL_{savstr}_IR3_v09.png',format='png')

# %%

# plot dark1 and dark2

#### BGR CHANNELS

if NIGHTGLOW:
    fig, axs = pplt.subplots(figwidth='16cm',figheight='8cm', ncols=2, nrows=1,abc='a.',
                                includepanels=False)
    #pxs=axs.panel('r', space=0, width='4em')

    fig.suptitle('dark1 and dark2')

#fig.format(suptitle=(day0.strftime('%Y-%m-%d') + " (Pass #"+str(pass_no)+" - ODIN time: "+str(dstest.time[pass_no-1].values)))

o = 0

for results in results_dict.values(): # for every day
                
    day0 = results['days'][0]
    dstest = results['ds'][0]

    if len(results['Ma1']) > 0: # if there are any encounters
        
        for i in range(0,len(dstest.time)): # for every pass
            pass_no=i+1

            dark1 = np.array(results['dark1'][i])
            dark2 = np.array(results['dark2'][i])
            
            for k in range(0,1):
                counts = np.linspace(0,len(dark1[:,k]))
                axs[0].scatter(dark1[:,k], alpha=1,label='dark1',
                    color='red',linewidth=0.75)
                axs[1].scatter(dark2[:,k], alpha=1,label='dark2',
                    color='blue',linewidth=0.75)
                
                
                
# %%
