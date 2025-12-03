#%%
#CMTM model
import xarray as xr
import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd


#%%
def CMTM_generate_old(urnal, diurnal_waves,LST, altitude, mon, variables):
    datadir = '/Users/lindamegner/MATS/MATS-retrieval/data/CMTM5541913/'
    datax = xr.Dataset()
    for variable in variables:
        variablename='_'+variable
        first=True

        for urnal_str in urnal:

            for wave in diurnal_waves:

                #if (wave == 'w1' and urnal_str == 'diurnal') or (wave == 'w2' and urnal_str == 'semidiurnal'):
                    #print(f'SKIPPED! {wave} and {urnal_str}')
                #    b=5

                #else:

                lon=np.linspace(-np.pi,np.pi,360)

                diurnal = xr.open_dataset(f'{datadir}ctmt_{urnal_str}_2002_2008.nc')
                diurnal = diurnal.sel(lev=altitude)
                diurnal = diurnal.sel(month=mon)

                phasevar=f"phase_{wave}{variablename}"
                ampvar=f"amp_{wave}{variablename}"

                try:

                    amp = diurnal[ampvar]
                    units = amp.units
                    phase = diurnal[phasevar]

                    amp = np.resize(diurnal[ampvar],(len(lon),len(diurnal[ampvar])))
                    phase = np.resize(diurnal[phasevar],(len(lon),len(diurnal[phasevar])))
                    lon = np.resize(lon, (len(diurnal[phasevar]),len(lon)))

                    if wave[0] == "w":
                        s = - int(wave[1])
                        #print(f's = {s}')
                    if wave[0] == "e":
                        s = int(wave[1])
                        #print(f's = {s}')
                    
                    if wave[0] == "s":
                        s = 0
                        #print(f's = {s}')

                    if urnal_str == 'diurnal':
                        n = 1

                    else:
                        n = 2

                    if first:
                        data=amp.T*np.cos(n*2*np.pi/24*(LST-phase.T)-(s+n)*lon)
                        first=False
                    else:
                        #print(wave)
                        data=data+amp.T*np.cos(n*2*np.pi/24*(LST-phase.T)-(s+n)*lon)


                except:
                    #print(f'No {wave} in {urnal_str}')
                    b=5

        datax[variable] = (('lat','lon'), data)
        datax[variable].attrs["units"] = units

    datax['lon'] = lon[0,:]/np.pi*180
    datax['lat'] =  diurnal.lat

    return datax


def CMTM_generate(urnal, diurnal_waves, LST, altitude, mon, variables):
    import numpy as np
    import xarray as xr

    datadir = '/Users/lindamegner/MATS/MATS-retrieval/data/CMTM5541913/'

    if isinstance(altitude, (float, int)):
        altstart = altitude - 0.5
        altstop = altitude + 0.5
    elif isinstance(altitude, list):
        altstart = altitude[0]
        altstop = altitude[-1]
    else:
        raise Exception('Altitude must be a float, int, or list of [min, max]')

    datax = xr.Dataset()

    for variable in variables:
        variablename = '_' + variable
        first = True

        for urnal_str in urnal:
            for wave in diurnal_waves:
                lon = np.linspace(-np.pi, np.pi, 360)

                try:
                    diurnal = xr.open_dataset(f'{datadir}ctmt_{urnal_str}_2002_2008.nc')
                    diurnal = diurnal.sel(lev=slice(altstart, altstop))
                    diurnal = diurnal.sel(month=mon)

                    ampvar = f"amp_{wave}{variablename}"
                    phasevar = f"phase_{wave}{variablename}"

                    amp = diurnal[ampvar]
                    phase = diurnal[phasevar]
                    units = amp.units

                    lev_len = len(diurnal.lev)
                    lat_len = len(diurnal.lat)
                    lon_len = len(lon)

                    if wave[0] == "w":
                        s = -int(wave[1])
                    elif wave[0] == "e":
                        s = int(wave[1])
                    else:
                        s = 0

                    n = 1 if urnal_str == 'diurnal' else 2

                    # Initialize data array
                    if first:
                        data = np.zeros((lev_len, lat_len, lon_len))
                        first = False

                    # Compute wave contribution for each latitude
                    for i in range(lat_len):
                        amp_slice = amp[:, i].values  # shape (lev,)
                        phase_slice = phase[:, i].values  # shape (lev,)
                        for j in range(lev_len):
                            data[j, i, :] += amp_slice[j] * np.cos(
                                n * 2 * np.pi / 24 * (LST - phase_slice[j]) - (s + n) * lon
                            )

                except Exception as e:
                    print(f"Skipped {wave} in {urnal_str}: {e}")

        # Assign to xarray DataArray
        datax[variable] = xr.DataArray(
            data,
            dims=('lev', 'lat', 'lon'),
            coords={
                'lev': diurnal.lev,
                'lat': diurnal.lat,
                'lon': lon / np.pi * 180
            },
            attrs={"units": units}
        )

    datax = datax.assign_coords(lon=lon / np.pi * 180)
    datax = datax.assign_coords(lat=diurnal.lat)

    return datax



def generate_CMTM_lookup_table(urnal, diurnal_waves, altitude, mon, variables, interval_minutes=5):
    """
    Generates a combined xarray.Dataset for every 5-minute LST interval.
    
    Returns:
    - xarray.Dataset with dimensions (LST, lat, lon)
    """
    import numpy as np
    import xarray as xr

    LST_values = np.arange(0, 24, interval_minutes / 60)  # LST in hours
    datasets = []

    for LST in LST_values:
        print(f"Processing LST = {LST:.2f} hours")
        datax = CMTM_generate(urnal, diurnal_waves, LST, altitude, mon, variables)
        if datax is not None:
            datax = datax.expand_dims(LST=[LST])
            datasets.append(datax)

    combined = xr.concat(datasets, dim="LST")
    return combined


def compute_phase(s, n, GMT, longitude, phase_offset=0):
    """
    Compute atmospheric tide phase for arrays of time and longitude.

    Parameters:
    - s: zonal wavenumber
    - n: harmonic of solar day
    - GMT Greenwich Mean Time (GMT from 0 to 24)
    - longitude longitudes in degrees
    - phase_offset: optional phase offset (default 0)

    Returns:
    - xarray.DataArray of phase values
    """

    # Earth's rotation rate in degrees/day
    Omega = 2*np.pi

    # Compute phase
    phase = n * Omega * GMT/24. + s * np.deg2rad(longitude) - phase_offset
    

    # Wrap to [0, 360)

    phase_deg =np.rad2deg(phase)%360

    return phase_deg

def compute_utc_hour_from_LST(LST, lon):
    """
    Computes the Greenwich Mean Time (GMT) from Local Time (LST) and longitude.

    Parameters:
    - LST (float): Local Sidereal Time in hours (0 to 24)
    - lon (float): Observer's longitude in degrees (positive east of Greenwich, negative west)

    Returns:
    - float: Greenwich Mean Time in hours (0 to 24)
    """
    # Convert longitude from degrees to hours (360° = 24h, so 15° = 1h)
    lon_in_hours = lon / 15.0

    # Subtract longitude in hours from LST to get GMST
    GMT = LST - lon_in_hours

    # Normalize the result to be within the 0–24 hour range
    GMT = GMT % 24

    return GMT    




def get_utc_hour(timestamp):
    """
    Extracts UTC hour from a datetime timestamp.

    Parameters:
    - timestamp (datetime64 or datetime): UTC time

    Returns:
    - float: Time in hours (e.g., 18.5 for 18:30)
    """
    # Convert to pandas Timestamp (handles numpy.datetime64 safely)
    ts = pd.to_datetime(timestamp)

    hour = ts.hour
    minute = ts.minute
    second = ts.second + ts.microsecond / 1e6

    return hour + minute / 60 + second / 3600



#%%
altitude=[80, 90]
CMTMds3d = generate_CMTM_lookup_table(
    urnal=['diurnal','semidiurnal'], #this translates to n=1 (diurnal) and n=2 (semidiurnal) below
    diurnal_waves=['w1','w2','w3','w4','e1','e2','e3','e4','s0'],
    #diurnal_waves=['w2'],
    altitude=altitude,
    mon=3,
    variables=["w", "t"]
)

print(CMTMds3d)


#%%
var='w'
CMTMds3d[var].sel(LST=17.25).mean(dim='lon').plot(x='lat', y='lev')

#%%
for lat in CMTMds3d.lat.values[::2]:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # Morning (LST=5.25)
    CMTMds3d[var].sel(LST=5.25, lat=lat).plot(x='lon', y='lev', ax=ax1)
    ax1.set_title(f'CMTM w at LST=5.25 (morning), lat={lat:.2f}')
    ax1.set_xlabel('Longitude (deg)')
    ax1.set_ylabel('Altitude (lev)')

    # Evening (LST=17.5)
    CMTMds3d[var].sel(LST=17.5, lat=lat).plot(x='lon', y='lev', ax=ax2)
    ax2.set_title(f'CMTM w at LST=17.5 (evening), lat={lat:.2f}')
    ax2.set_xlabel('Longitude (deg)')
    ax2.set_ylabel('Altitude (lev)')

    plt.tight_layout()
    plt.show()
#CMTMds3d['w'].sel(LST=17.25, lat=0.0).plot(x='lon', y='lev')


#ds3d['VER'].sel(img_time=slice("2023-03-06T8:00:43", "2023-03-06T19:30:43")).plot(x='img_time', y='alt_coord', cmap='inferno', robust=True)


#%% Old stuff: Fetches level 1b MATS data
#ds = fetch_MATS_l1b_data('IR1')
# Plot the 1000th image using Xarrays built in plotting fuction
#ds.ImageCalibrated.isel(time=1000).plot()


# Path to your NetCDF file
file_path = "/Users/lindamegner/MATS/MATS-retrieval/Lukas1DRegressionOfL1b/linear_L2_2023_03_05-15.nc"
# Load the dataset
ds3d = xr.open_dataset(file_path)
# Display basic info about the dataset
print(ds3d)

#Separate ascending and descending node
# Compute latitude difference along time
lat = ds3d['latitude']
lat_diff = lat.diff(dim='img_time')
# Pad to match original length (prepend NaN)
lat_diff_padded = xr.concat([xr.full_like(lat.isel(img_time=0), np.nan), lat_diff], dim='img_time')
# Create node labels
node_labels = xr.where(lat_diff_padded > 0, 'ascending', 'descending')
# Assign to dataset
ds3d['node'] = node_labels
#ds3d = ds3d.where(ds3d['node'] == 'ascending', drop=True)


#%%

# Comput GMT in hours
ds3d['GMT'] = xr.apply_ufunc(
    get_utc_hour,
    ds3d['img_time'],
    vectorize=True,
    output_dtypes=[float]
)


CMTMds3d['GMT']=compute_utc_hour_from_LST(CMTMds3d['LST'], CMTMds3d['lon'])




#%%
#ds3d['VER'].sel(img_time=slice("2023-03-06T8:00:43", "2023-03-06T19:30:43")).plot(x='img_time', y='alt_coord', cmap='inferno', robust=True)
ds3d['VER'].plot(x='img_time', y='alt_coord', cmap='inferno', robust=True)


#%%
ds3d.img_time.sel(img_time=slice("2023-03-12T8:00:43", "2023-03-13T09:30:43")).plot(marker='.', linestyle='None')

#%%
if False: #Change to true to plot as fucntion of longitude
    # select equator and plot as function of longitude and altitude
    ds3deq = ds3d.where((ds3d['latitude'] > 60) & (ds3d['latitude'] < 90), drop=True)
    #ds3d_eq['VER'].plot(x='img_time', y='alt_coord', cmap='inferno', robust=True)

    # Flatten and filter valid data
    lon_flat = ds3deq['longitude'].values.flatten()
    alt_flat = ds3deq['altitude'].values.flatten()
    VER_flat = ds3deq['VER'].values.flatten()

    # Filter out NaNs
    mask = np.isfinite(lon_flat) & np.isfinite(alt_flat) & np.isfinite(VER_flat)
    lon_valid = lon_flat[mask]
    alt_valid = alt_flat[mask]
    VER_valid = VER_flat[mask]

    # Create a grid using interpolation (optional but helpful for contours)
    from scipy.interpolate import griddata

    # Define grid resolution
    lon_grid = np.linspace(np.min(lon_valid), np.max(lon_valid), 300)
    alt_grid = np.linspace(np.min(alt_valid), np.max(alt_valid), 300)
    lon_mesh, alt_mesh = np.meshgrid(lon_grid, alt_grid)

    # Interpolate VER onto the grid
    VER_grid = griddata((lon_valid, alt_valid), VER_valid, (lon_mesh, alt_mesh), method='linear')

    # Plot
    plt.figure(figsize=(10, 6))
    contour = plt.contourf(lon_mesh, alt_mesh, VER_grid, levels=100, cmap='viridis')
    plt.colorbar(contour, label='VER')
    plt.xlabel('Longitude')
    plt.ylabel('Altitude (m)')
    plt.title('VER vs Longitude and Altitude (Contour)')
    plt.tight_layout()
    plt.show()




#%%

ds=ds3d.sel(alt_coord=slice(altitude[0]*1000, altitude[-1]*1000)).mean(dim='alt_coord')
CMTMds=CMTMds3d.sel(lev=slice(altitude[0], altitude[-1])).mean(dim='lev')
CMTMds['w'] = -CMTMds['w']  #invert winds to be able to directly compare with signal strenght
#dstry=ds3d.sel(alt_coord=slice(altitude[0], altitude[-1]))


#%%
for s in range(-4,-3,1):
    # Example parameters
    #s = 1# zonal wavenumber
    n = 2  # harmonic of solar day


    # Compute phase
    ds['phase'] = compute_phase(s, n, ds['GMT'], ds['longitude'])
    CMTMds['phase']=compute_phase(s, n, CMTMds['GMT'],  CMTMds['lon'])



    # # Select MATS data
    # ds_sel = ds.sel(time=slice("2023-02-11", "2023-02-13"))
    # ds_selected=ds_sel.where((ds_sel['TPlat'] >= 59) & (ds_sel['TPlat'] <= 60), drop=True)
    # ds_selected['signal'] = ds_selected['ImageCalibrated'].isel(im_row=140, im_col=22)


    # Set plot dimension mode
    plotdim = 2

    # Define phase bins
    bins = np.linspace(0, 360, 73)  # 72 bins of 5 degrees
    labels = (bins[:-1] + bins[1:]) / 2

    # Define latitude bins for MATS
    minlat=-90
    maxlat=90
    lat_bins = np.linspace(minlat, maxlat, int((maxlat-minlat)/5))
    lat_labels = (lat_bins[:-1] + lat_bins[1:]) / 2


    # --- Select data ---
    if plotdim == 1:
        mylat=60
        CMTMds_sel = CMTMds.sel(lat=mylat)
        ds_selected = ds.where((ds['latitude'] >= mylat-1) & (ds['latitude'] <= mylat+1), drop=True).copy()
        signal_raw = ds_selected['VER'] #.sel(alt_coord=altitude * 1000)
    elif plotdim==2:
        CMTMds_sel = CMTMds
        ds_selected = ds.where((ds['latitude'] >= minlat) & (ds['latitude'] <= maxlat), drop=True).copy()
        signal_raw = ds_selected['VER']
    else:
        Exception('Set plotdim to 1 or 2!')

    # --- Normalize to zero mean and unit std (same std for both) ---
    ds_selected['signal'] = (signal_raw - signal_raw.mean(skipna=True)) / signal_raw.std(skipna=True)
    CMTMds_sel['signal'] = (CMTMds_sel[var] - CMTMds_sel[var].mean(skipna=True)) / CMTMds_sel[var].std(skipna=True)
  

    # --- Create Figure ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), constrained_layout=True)

    # --- MATS Plot ---
    if plotdim == 1:
        print('Note 1d is not the demeaned with regards to latitude yet as 2d is')
        phase_mats = ds_selected['phase'].sel(alt_coord=altitude * 1000)
        mean_by_phase_mats = ds_selected['signal'].groupby_bins(phase_mats, bins=bins, labels=labels).mean()
        mean_by_phase_mats = mean_by_phase_mats.reindex({'phase_bins': labels})
        ax1.plot(labels, mean_by_phase_mats, marker='o', color='tab:blue')
        ax1.set_xlabel("Phase (degrees)")
        ax1.set_ylabel("Signal")
        ax1.set_title("MATS: Signal vs Phase, s=" +str(s)+ ' n = '+str(n))
        ax1.grid(True)



    elif plotdim == 2:
        # Flatten arrays
        phase_flat = ds_selected['phase'].values.flatten()
        signal_flat = ds_selected['signal'].values.flatten()
        lat_flat = ds_selected['latitude'].values.flatten()
        lon_flat = ds_selected['longitude'].values.flatten()

        # Create DataFrame
        df_mats = pd.DataFrame({
            'phase': phase_flat,
            'signal': signal_flat,
            'latitude': lat_flat,
            'longitude': lon_flat
        }).dropna()

        # Bin phase and latitude
        df_mats['phase_bin'] = pd.cut(df_mats['phase'], bins=bins, labels=labels)
        df_mats['lat_bin'] = pd.cut(df_mats['latitude'], bins=lat_bins, labels=lat_labels)



        # Remove mean signal per latitude bin
        df_mats['signal_demeaned'] = df_mats.groupby('lat_bin')['signal'].transform(lambda x: x - x.mean())
        # --- Normalize to zero mean and unit std (same std for both) ---
        df_mats['signal_scaled'] = (df_mats['signal_demeaned'] - df_mats['signal_demeaned'].mean(skipna=True)) / df_mats['signal_demeaned'].std(skipna=True)

        # Group and pivot the demeaned signal
        mean_by_phase_lat = df_mats.groupby(['lat_bin', 'phase_bin'])['signal_scaled'].mean().unstack()
        #mean_by_phase_lat = df_mats.groupby(['lat_bin', 'phase_bin'])['signal'].mean().unstack()
        Z_mats = mean_by_phase_lat.values
        X_mats, Y_mats = np.meshgrid(labels, lat_labels)

        # Plot
        # Compute standard deviation
        sigma = np.nanstd(Z_mats)
        # Plot with fixed color scale
        pcm1 = ax1.pcolor(X_mats, Y_mats, Z_mats, cmap='viridis', vmin=-3*sigma, vmax=3*sigma)
        fig.colorbar(pcm1, ax=ax1)

        # Axis labels and title
        ax1.set_title(f"MATS: Signal vs Phase and Latitude, s={s} n={n}")
        ax1.set_xlabel("Phase (degrees)")
        ax1.set_ylabel("Latitude")

    # --- CMTM Plot ---
    phase_cmtm = CMTMds_sel['phase']
    mean_by_phase_cmtm = CMTMds_sel['signal'].groupby_bins(phase_cmtm, bins=bins, labels=labels).mean()
    mean_by_phase_cmtm = mean_by_phase_cmtm.reindex({'phase_bins': labels})

    if plotdim == 1:
        ax2.plot(labels, mean_by_phase_cmtm, marker='o', color='tab:blue')
        ax2.set_xlabel("Phase (degrees)")
        ax2.set_ylabel("Signal")
        ax2.set_title("CMTM: Signal vs Phase")
        ax2.grid(True)

    elif plotdim == 2:

        # Convert phase_bins and lat to numeric values for axis ticks
        phase_ticks = mean_by_phase_cmtm['phase_bins'].values.astype(float)
        lat_ticks = mean_by_phase_cmtm['lat'].values.astype(float)

        # Create meshgrid for plotting
        X_cmtm, Y_cmtm = np.meshgrid(phase_ticks, lat_ticks)

        pcm2 = ax2.pcolor(X_cmtm, Y_cmtm, mean_by_phase_cmtm.values, cmap='viridis')
        fig.colorbar(pcm2, ax=ax2)
        ax2.set_title("CMTM: Signal vs Phase")
        ax2.set_xlabel("Phase (degrees)")
        ax2.set_ylabel("Latitude")
        ax2.set_xticks(phase_ticks[::8])  # Show every 8th tick for clarity
        ax2.set_yticks(lat_ticks[::4])    # Show every 4th tick for clarity


    plt.show()

# %%


urnal=['semidiurnal']  #this translates to n=1 (diurnal) and n=2 (semidiurnal) below
diurnal_waves=['w1','w2','w3','w4','e1','e2','e3','s0']
#diurnal_waves=['w2'],

mon=3
variables=[var]

LST=17.5
#altitude=[84,91]
lon=np.linspace(-np.pi,np.pi,360)

for n in range(1,3):
    if n==1: urnal=['diurnal']
    if n==2: urnal=['semidiurnal']
    for s in range(-3,4):
        if s==-3: diurnal_waves=['e3']
        if s==-2: diurnal_waves=['e2']
        if s==-1: diurnal_waves=['e1']
        if s==0: diurnal_waves=['s0']
        if s==1: diurnal_waves=['w1']
        if s==2: diurnal_waves=['w2']
        if s==3: diurnal_waves=['w3']

        try:
            datax=CMTM_generate(urnal, diurnal_waves,LST, altitude, mon, variables)
            datax = datax.mean(dim='lev')


            for i, variable in enumerate(variables):
                fig, ax = plt.subplots()
                m = ax.contourf(lon/np.pi*180, datax.lat, datax[variable], robust=98, extend='both')
                fig.colorbar(m, ax=ax, label=variable)
                ax.set_title('s= '+str(s)+ ' n= '+str(n)+' '+variables[i])
        except: print('crash for s= '+str(s)+'n=' +str(n))




# %%
altitude=[80,100]
urnal=['semidiurnal']
diurnal_waves=['e2']
datax=CMTM_generate(urnal, diurnal_waves,LST, altitude, mon, variables)
# %%
