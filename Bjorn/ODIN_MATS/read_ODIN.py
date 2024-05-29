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


#%%
def calc_distance(coord1, coord2):
    '''function to calculate distance between the two satellite positions'''
    distance = geodesic(coord1, coord2).kilometers
    return distance # km

def prepare_measurment(ch,ir3,ir4,subtract=True,endcut=-25):

    z1,p1=prepare_profile(ch)
    _,p3=prepare_profile(ir3)
    _,p4=prepare_profile(ir4)

    if subtract:

        p3 = p3[10:endcut]
        p4 = p4[10:endcut]

        p3=p3-p3[-4:].mean()/1.05
        p4=p4-p4[-4:].mean()/1.05

        p1=p1[10:endcut]-(p3+p4)/2

    return np.array(z1),np.array(p1)

def prepare_profile(ch):
    # This function averages some columns and
    # calculates tangent heights
    
    image = image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))

    # multiply with factor to get right values (depends on version?)
    profile = np.array(image[:, col-1:col+1].mean(axis=1)*1e12)

    # set heights
    common_heights = np.arange(60000,110250,250)
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile


#%%

# NLC
NLC=False

if NLC:
    data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/Odin/OSIRIS_NLCscans.mat')

# AIRGLOW
else:    
    data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/Odin/sm_ABand_2022_2023.mat')
# %%

# sattime lat lon
if NLC:
    lat=np.array(data['lat'])
    lon=np.array(data['lon'])
    time=np.array(data['t'])

    flat_list = []
    for xs in time:
        for x in xs:
            flat_list.append(x)
    timeO=flat_list

    flat_list = []
    for xs in lat:
        for x in xs:
            flat_list.append(x)
    lat=flat_list

    flat_list = []
    for xs in lon:
        for x in xs:
            flat_list.append(x)
    lon=flat_list

    timeO = np.array(timeO)
    lon = np.array(lon)
    lat = np.array(lat)

else:
    lat=np.array(data['sm']['lat'])
    lon=np.array(data['sm']['lon'])
    sza=data['sm']['sza']
    timeO=np.array(data['sm']['mjd'])

    for i in range(0,len(sza)):
        if type(sza[i]) is list:
            sza[i] = np.nan 

    sza = np.array(sza)

    time = np.array(timeO)

jd = timeO + 2400000.5
t = Time(jd, format='jd')
timeO = t.to_datetime()

# xarray
if NLC:
    ds_day_all = xr.Dataset(
        data_vars=dict(
            lat=(["time"], lat),
            lon=(["time"], lon)
        ),
        coords=dict(
            time=timeO,
        ),
        attrs=dict(description="ODIN"),
    )

else:
    ds = xr.Dataset(
        data_vars=dict(
            lat=(["time"], lat),
            lon=(["time"], lon),
            sza=(["time"], sza)
        ),
        coords=dict(
            time=timeO,
        ),
        attrs=dict(description="ODIN"),
    )
    ds_night_all=ds.where(ds.sza > 90)
    ds_day_all=ds.where(ds.sza < 90)

if NLC:
    ds_day_all=ds_day_all.sel(time = slice("2022-02-01", "2023-03-01")).dropna(dim='time')

else:
    ds_night_all=ds_night_all.sel(time = slice("2022-02-01", "2023-03-01")).dropna(dim='time')
    ds_day_all=ds_day_all.sel(time = slice("2022-02-01", "2023-03-01")).dropna(dim='time')


#%%
## COVERAGE MAP
fig, axs = pplt.subplots(figwidth='15cm',ncols=1, nrows=1,abc='a.',sharex=0,proj='cyl')
axs.format(coast=True,landzorder=5)
axs.scatter(ds_day_all.lon,ds_day_all.lat, s=0.3, c='red')
if not NLC:
    axs.scatter(ds_night_all.lon,ds_night_all.lat, s=0.3, c='blue')
    axs.format(title='OSIRIS 22/23 airglow coverage --- red: dayglow; blue: nightglow')
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/dayglow_nightglow_map.png',format='png')
else:
    axs.format(title='OSIRIS 22/23 NLC coverage')
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/nlc_map.png',format='png')

### NLC DAYGLOW AIRGLOW COVERAGE
fig, axs = pplt.subplots(figwidth='20cm',ncols=3, nrows=1,abc='a.',sharex=0)
axs.format(coast=True,landzorder=5)


if not NLC:
    m=axs[0].scatter(ds.time,ds.lat, c=lon, s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    fig.colorbar(m,loc='r',title='longitude')

    axs[1].scatter(ds_night_all.time,ds_night_all.lat, c=ds_night_all.lon, s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    axs[2].scatter(ds_day_all.time,ds_day_all.lat, c=ds_day_all.lon, s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    axs[0].format(title='airglow coverage')
    axs[1].format(title='nightglow coverage')
    axs[2].format(title='dayglow coverage')
    fig.format(suptitle='OSIRIS 22/23 MLT AIRGLOW COVERAGE')
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/glow_emissions.png',format='png')

else:
    axs[0].format(title='NLC coverage')
    axs[0].scatter(ds_day_all.time,ds_day_all.lat, c=ds_day_all.lon, s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/NLC.png',format='png')

# %% 
# example of a single pass
#starttime=datetime(2023,2,15,18,30)
#stoptime=datetime(2023,2,15,18,34)

starttime=datetime(2023,2,1,0,0)
stoptime=datetime(2023,3,1,0,0)


channel='IR1'
l1b_version="0.6"

run=True

if run:
# some filter settings - dayglow (ALL SZA)
    if NLC:
        tplat0, tplat1 = -90, -40
        ccdsel0, ccdsel1 = 6, 6
        filter={'TPlat': [tplat0, tplat1], 'CCDSEL': [ccdsel0,ccdsel1]}
        dftop_all=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)

    else:
        dmin,dmax = 0, 95
        tplat0, tplat1 = -90, -70
        ccdsel0, ccdsel1 = 1, 1
        filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1], 'CCDSEL': [ccdsel0,ccdsel1]}
        dftop_all=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)

    dftop_all['EXPDate'] = pd.to_datetime(dftop_all['EXPDate'])
#%%
if run:
    fig, axs = pplt.subplots(figwidth='10cm',ncols=1, nrows=1,abc='a.',sharex=0)
    axs.format(coast=True,landzorder=5)
    m=axs[0].scatter(dftop_all.EXPDate,dftop_all.TPlat, c=dftop_all.TPlon, s=0.5,alpha=1,cmap='viridis',vmax=360,vmin=0)
    axs[0].format(title='MATS coverage')
    fig.colorbar(m,loc='r',title='longitude')

    if not NLC:
        fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/MATS_coverage.png',format='png')
    else:
        fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/MATS_coverage_NLC.png',format='png')

#%% 

# read in BGR channels
ccdsel0, ccdsel1 = 2, 3
filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1], 'CCDSEL': [ccdsel0,ccdsel1]}
df_bgr =read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)

ir3 = df_bgr[df_bgr['channel'] == 'IR3'].dropna().reset_index()#[0:10]
ir4 = df_bgr[df_bgr['channel'] == 'IR4'].dropna().reset_index()#[0:10]

ir3=ir3.set_index('EXPDate')
ir4=ir4.set_index('EXPDate')

# find nearest times
#ir3_times, ir4_times = [], []
#ir3.index.get_loc(MATS_times, method='nearest')
#ir4.index.get_loc(MATS_times[0], method='nearest')
#% DEV


#%%
if run:
    dftop_all = dftop_all.set_index('EXPDate')
# RUN DAY-BY-DAY
interpolate=False
polar=True
dtime = 2 # minutes before and after Odin MEAS
mindis=500


for dayi in range(0,31):
    MATS_times, ODIN_times = [], []
    day0 = starttime + DT.timedelta(dayi)
    day1 = day0 + DT.timedelta(1)

    ds_day=ds_day_all.sel(time=slice(day0,day1))
    dftop=dftop_all.loc[str(day0)[0:-9]:str(day1)[0:-9]]

    if (len(dftop.index) > 0) and (len(ds_day.time) > 0):
        # CONVERT TO EPOCH TIME
        test_EXPdate = np.zeros(len(dftop.index.values))
        test_odint = np.zeros(len(ds_day.time))

        for i in range(0,len(dftop.index.values)):
            test_EXPdate[i]=int(dftop.index.values[i])
        for i in range(0,len(ds_day.time)):
            test_odint[i]=int(ds_day.time.values[i])

        if interpolate:
            test_ttot=np.sort(np.concatenate((test_EXPdate,test_odint)))
            fig, axs = pplt.subplots(figwidth='20cm',ncols=1, nrows=1,abc='a.',sharex=0)
            axs.plot(test_ttot)
            axs.format(title='union time')

            # BEFORE INTERPOLATE - REMOVE TOO LARGE STEPS
            j=0
            timz_epoch=np.array([])
            timz=[]
            indz=[]
            for Oi in test_odint:
                tind = np.where(test_ttot == Oi)[0][0]
                
                # check that time between Oi and previous tot time < 5 min
                if ((Oi - test_ttot[tind-1])/1e9 > 5*60):
                    timz_epoch=np.append(timz_epoch, Oi)
                    timz.append(ds_day.time[j].values)
                    indz.append(tind)
                    #test_odint = np.delete(test_odint,j)
                elif len(test_ttot) != tind+1:
                    if ((test_ttot[tind+1] - Oi)/1e9 > 5*60):
                        timz_epoch=np.append(timz_epoch, Oi)
                        timz.append(ds_day.time[j].values)
                        indz.append(tind)

                j = j + 1

            print(f'{len(ds_day.time)} ODIN times before')
            for tim in timz:
                # remove corresponding in ds_day
                ds_day = ds_day.where(ds_day.time != tim,drop=True)
            print(f'{len(ds_day.time)} ODIN times after')

            # recreate in case needed
            test_odint = np.zeros(len(ds_day.time))
            for i in range(0,len(ds_day.time)):
                test_odint[i]=int(ds_day.time.values[i])

            test_ttot=np.sort(np.concatenate((test_EXPdate,test_odint)))

        else:
            test_ttot=test_EXPdate
            test_EXPdate_str = []
            for i in range(0,len(dftop.index.values)):
                test_EXPdate_str.append(dftop.index[i])

        # INTERPOLATE
        if interpolate:
            TPlati=np.interp(test_ttot,test_EXPdate,dftop.TPlat.values)
            TPloni=np.interp(test_ttot,test_EXPdate,dftop.TPlon.values)
        else:
            TPlati = dftop.TPlat.values
            TPloni = dftop.TPlon.values

        print(test_EXPdate_str)

        # XARRAY STUFF
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
            attrs=dict(description="MATS on union time grid"),
        )

        # MAX TIME BETWEEN MATS AND ODIN MEASUREMENT

        # GET RID OF DUBLICATED TIMES (instrument anomalies?)
        MATSint = MATSint.drop_duplicates(dim = 'time')

        # ISOLATE TEN MINUTES BEFORE AND AFTER
        time0=((ds_day.time[0] - np.timedelta64(dtime, 'm')).values)
        time1=((ds_day.time[0] + np.timedelta64(dtime, 'm')).values)
        MATSi=MATSint.sel(time=slice(int(time0), int(time1)),drop=True)
        for i in range(0, len(test_odint)):
            time0=((ds_day.time[i] - np.timedelta64(dtime, 'm')).values)
            time1=((ds_day.time[i] + np.timedelta64(dtime, 'm')).values)
            MATSi=xr.concat([MATSi,MATSint.sel(time=slice(int(time0), int(time1)),drop=True)],dim='time')
        MATSi = MATSi.drop_duplicates(dim = 'time').dropna(dim='time')

        if len(MATSi.time) > 0:
            # INTERPOLATED PLOTS
            if interpolate:
                fig, axs = pplt.subplots(figwidth='20cm',ncols=1, nrows=2,abc='a.',sharex=0)
                axs[0].scatter(test_odint,ds_day.lat.values, c='blue', s=1,alpha=1)
                axs[0].scatter(MATSi.time.values,MATSi.TPlati.values, c='red', s=1,alpha=1)
                axs[1].scatter(test_odint,ds_day.lon.values, c='blue', s=1,alpha=1)
                axs[1].scatter(MATSi.time.values,MATSi.TPlatilon.values, c='red', s=1,alpha=1)

            # MINIMUM DISTANCE BETWEEN ODIN AND MATS POINTS

            #MATS_measurements
            distance_bool = np.zeros([len(ds_day.time),len(MATSi.time)])
            for i in range(0,len(ds_day.time)):
                
                #isolate odin with nearby mats (in time)
                coord_ODIN = (ds_day.lat.values[i], ds_day.lon.values[i])
                time0=((ds_day.time[i] - np.timedelta64(dtime, 'm')).values)
                time1=((ds_day.time[i] + np.timedelta64(dtime, 'm')).values)
                MATS_dt=MATSi.where(MATSi.time > int(time0))
                MATS_dt=MATS_dt.where(MATS_dt.time < int(time1))

                for j in range(0,len(MATS_dt.time)):
                    if MATS_dt.TPlati[j].notnull():
                        coord_MATS = (MATS_dt.TPlati.values[j],MATS_dt.TPlatilon.values[j])
                        dis=calc_distance(coord_ODIN, coord_MATS)
                        if dis < mindis:
                            distance_bool[i,j] = 1

            MATS_meas = MATSi.where(np.sum(distance_bool,axis=0) > 0)
            ODIN_meas = ds_day.where(np.sum(distance_bool,axis=1) > 0)

            if (len(ODIN_meas.dropna('time').lat) > 0) and (len(MATS_meas.dropna('time').time) > 0):
                # PLOTTING OF MAPS
                time_arrow=True #requires interpolation = True

                if NLC or polar:
                    fig, axs = pplt.subplots(figwidth='20cm',ncols=1, nrows=1,abc='a.',sharex=0,proj='splaea')
                    axs.format(coast=True, boundinglat=-30,landzorder=5)

                else:
                    pplt.rc.axesfacecolor = 'gray4'
                    fig, axs = pplt.subplots(figwidth='15cm',ncols=1, nrows=1,abc='a.',sharex=0,proj='cyl')
                    axs.format(coast=True,landzorder=5)


                axs.scatter(ODIN_meas.dropna(dim='time').lon,ODIN_meas.dropna(dim='time').lat, s=15, c='blue',marker='x')

                for j in range(0,len(ds_day.time)):
                    first=True
                    append_ODIN=False
                    if np.sum(distance_bool[j,:],axis=0) > 0:

                        for i in range(0,len(MATSi.time),1):

                            if MATSi.time.values[i] == int(ds_day.time.values[j]):
                                if time_arrow:
                                    axs.plot([MATSi.TPlatilon[i], ds_day.lon[j]],
                                                [MATSi.TPlati[i], ds_day.lat[j]],
                                                c='green', alpha=0.75,
                                                linestyle='-',linewidth=0.75)

                            if distance_bool[j,i] == 1:
                                
                                append_ODIN=True
                                MATS_times.append(MATSi.time_str.values[i])

                                axs.plot([MATSi.TPlatilon[i], ds_day.lon[j]],
                                            [MATSi.TPlati[i],ds_day.lat[j]],
                                            c='grey6', alpha=0.6, linestyle='--',
                                            linewidth=0.5)
                                m=axs.scatter(MATSi.TPlatilon[i], MATSi.TPlati[i],
                                            c='red',
                                            s=5, vmin=-4, vmax=4,cmap='coral')
                                
                                if first:
                                    timediff= np.around((MATSi.time.values[i] - int(ds_day.time.values[j]))/(1e9*60),
                                                        decimals=2)
                                    axs.text(s=f'      dt: {str(timediff)} min',
                                                x=MATSi.TPlatilon.values[i],
                                                y=MATSi.TPlati.values[i],
                                                transform='map',fontsize=6,
                                                c='black', weight="bold",bbox=False,
                                                bboxalpha=0.1,bboxcolor='red')
                                    first=False
                                    axs.plot([MATSi.TPlatilon[i], ds_day.lon[j]],
                                                [MATSi.TPlati[i], ds_day.lat[j]], c='black',
                                                alpha=0.75, linestyle='--',linewidth=0.75)
                        if append_ODIN:
                            ODIN_times.append(ds_day.time.values[j])

                    # plot texts for ODIN points
                    for i in range(0,len(ODIN_meas.dropna(dim='time').time)):
                        axs.text(s='      '+str(ODIN_meas.dropna(dim='time').time.values[i])[0:-10],
                                    x=ODIN_meas.dropna(dim='time').lon.values[i],
                                    y=ODIN_meas.dropna(dim='time').lat.values[i],transform='map',fontsize=6,
                                    c='green', weight="bold",bbox=False,
                                    bboxalpha=0.1,bboxcolor='green')

                fig.format(suptitle=f'ODIN / MATS: {str(day0)} -- {str(day1)}')
                axs.format(title='red: MATS - blue: Odin //')

                if not NLC:
                    fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/map_mindis{mindis}_mintime{dtime}_{str(day0)}.png',format='png')
                else:
                    fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/NLC/map_mindis{mindis}_mintime{dtime}_{str(day0)}.png',format='png')


            #  GENERATE MEANS FROM THE ISOLATED DATA

            # mean of all OSIRIS measurements done at ODIN_times (epoch)
            # RECREATE DS now with profiles included:

            # OSIRIS NEEDS TO BE PUT ON A COMMON GRID
            # 60 to 105?

            OSheights = np.arange(60,110,0.25)

            # messy OSIRIS stuff
            OSinter=np.zeros([len(data['sm']['mjd']),len(OSheights)])
            for i in range(0,len(data['sm']['mjd'])):
                flip=False # for sorting

                zt = data['sm']['zt'][i]
                flat_list = []
                for xs in zt:
                    for x in xs:
                        flat_list.append(x)
                zt=flat_list
                if flat_list[0] > flat_list[1]:
                    zt=zt[::-1]
                    flip=True

                Li = data['sm']['L1'][i]
                flat_list = []
                for xs in Li:
                    for x in xs:
                        flat_list.append(x)
                Li=flat_list
                if flip:
                    Li=Li[::-1]
                    

                OSinter[i,:]=np.interp(OSheights,zt,Li)

            dstest = xr.Dataset(
                data_vars=dict(
                    lat=(["time"], lat),
                    lon=(["time"], lon),
                    sza=(["time"], sza),
                    Li=(["time", "altitude"], OSinter),

                ),
                coords=dict(
                    time=timeO,
                    altitude=OSheights
                ),
                attrs=dict(description="ODIN"),
            )
            dstest=dstest.sel(time=ODIN_times) # xarray
            dstest.Li.mean(dim='time').plot.line()

            common_heights = np.arange(60000,110250,250)

            # compute means encounter by encounter (MATS)
            pass_no=0
            print(f'compute means encounter by encounter ({len(ODIN_times)} encounters)')
            for j in range(0,len(distance_bool[:,0])):

                MATS_times_tmp=MATSi.time_str.where(distance_bool[j,:] > 0).dropna(dim='time')
                if len(MATS_times_tmp) > 0:
                    pass_no=pass_no+1
                    # mean of all MATS measurements done at MATS_times (epoch)
                    MATS_df=dftop_all.loc[MATS_times_tmp] # pandas
                    MATS_df['EXPDate'] = MATS_df.index


                    MALi = np.zeros([len(MATS_df), len(common_heights)-10-25])
                    for i in range(0,len(MATS_df)):
                        #heights, profile = prepare_profile(MATS_df.iloc[i])

                        # nearest time
                        ir3_p=ir3.iloc[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
                        ir4_p=ir4.iloc[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')

                        # set EXPDate accordingly
                        ir3_p['EXPDate'] = ir3.index[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
                        ir4_p['EXPDate'] = ir4.index[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]

                        #ir3_p=ir3_p.reset_index()
                        #ir4_p=ir4_p.reset_index()

                        zsi,MALi[i,:] = prepare_measurment(MATS_df.iloc[i], ir3_p, ir4_p,subtract=True,endcut=-25)

                    fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=1,abc='a.',sharey=0)
                    fig.format(suptitle=f'{day0}-{day1} (Pass #{pass_no})')
                    axs[0].plot(np.mean(MALi[:,:],axis=0), zsi[10:-25]/1000, color='red',label='MATS')
                    axs[0].plot(dstest.Li.values[pass_no-1,:], dstest.altitude.values,color='blue',label='OSIRIS')
                    axs[0].legend()
                    axs[0].format(title='Limb radiance (IR1)',ylabel='Tangent altitude [km]',
                                xlabel='[m-2 s-1 str-1 nm-1]')

                    axs[0].format(xlim=[0,1e15])

                    axs[1].format(title=str(dstest.time[pass_no-1].values))
                    print(MATS_df.EXPDate)
                    axs[1].plot(MATS_df.EXPDate,MATS_df.TPsza.values, color='red')
                    
                    axs[1].axhline(dstest.sza[pass_no-1], color='blue',linestyle='--',alpha=0.5)
                    axs[1].axvline(dstest.time[pass_no-1].values, color='blue',linestyle='--',alpha=0.5)
                    #axs[1].text(str(dstest.time[pass_no-1].values),fontsize=1,loc='l')
                    axs[1].scatter(dstest.time[pass_no-1],dstest.sza[pass_no-1], s=15, marker='x', color='blue')
                    #axs[1].format(title='solar zenith angles',xlim=)
                    fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/PROFILES/{day0}_pass{pass_no}.png',format='png')

# %%
