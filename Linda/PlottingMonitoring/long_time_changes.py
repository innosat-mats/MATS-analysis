# Investigate long time changes in of the MATS detectors
#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean, add_field_with_subtracted_rolling_mean2
import numpy as np
import matplotlib.pyplot as plt

from database_generation.experimental_utils import plot_CCDimage
from mats_utils.retrieval.error_flags import remove_flagged_images
from mats_utils.retrieval.hot_pix import compute_threshold_pixels

from mats_utils.data_availability.data_availability import load_hourly_time_spans, is_data_available  


#%%
#Function to add characteristics
def add_dark_mean(CCDitem):
    try:
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=130000)
        image = CCDitem.ImageCalibrated.copy()

        # Set pixels below threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[:threshold_row, col] = np.nan
        dark_mean = np.nanmean(image)
    except Exception as e:
        print(f"Error in add_dark_mean: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        dark_mean = np.nan
    return dark_mean
def add_dark_median(CCDitem):
    try:
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=130000)
        image = CCDitem.ImageCalibrated.copy()

        # Set pixels below threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[:threshold_row, col] = np.nan
        dark_median = np.nanmedian(image)+50
    except Exception as e:
        print(f"Error in add_dark_median: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        dark_median = np.nan
    return dark_median
def add_dark_count(CCDitem, threshold=10):
    try:
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=130000)
        image = CCDitem.ImageCalibrated.copy()

        # Set pixels below threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[:threshold_row, col] = np.nan
        dark_count = np.nansum(image > threshold)
    except Exception as e:
        print(f"Error in add_dark_count: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        dark_count = np.nan
    return dark_count
def add_bright_mean(CCDitem):
    try:
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=90000)
        image = CCDitem.ImageCalibrated.copy()
        maxrow = image.shape[0]
        # Set pixels above threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[threshold_row:maxrow, col] = np.nan
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=80000)
            # Set pixels below threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[:threshold_row, col] = np.nan
        bright_mean = np.nanmean(image)
    except Exception as e:
        print(f"Error in add_bright_mean: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        bright_mean = np.nan
    return bright_mean
def add_mean(CCDitem):
    try:
        image = CCDitem.ImageCalibrated.copy()
        mean = np.nanmean(image)
    except Exception as e:
        print(f"Error in add_mean: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        mean = np.nan
    return mean





#%%
available_times_l1b = '/Users/lindamegner/MATS/MATS-retrieval/MATS-utility-functions/data/data_availability/available_times_level1b-v1.0.1.txt'
# Load time spans
spans_l1b = load_hourly_time_spans(available_times_l1b)

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'
readindata=False


first_start_time = DT.datetime(2023, 5, 16, 12, 0, 0) #start with crop F and when bad pointing starts
final_stop_time = DT.datetime(2025, 4, 10, 12, 0, 0) #end with when satellite was turned off
version='1.0.1'
level='1b'

# pick time slots 
freq_days=10
freq=str(freq_days)+'D'
starttimes = pd.date_range(start=first_start_time, end=final_stop_time, freq=freq)
deltat_minutes=96

#datafilter = {'TPlat': [60, 70]}
name ='df_all_'+version+first_start_time.strftime('%Y%m%d_%H%M%S')+'_'+final_stop_time.strftime('%Y%m%d_%H%M%S')

#%%
if readindata:
    if 'df' in locals():
        del df

    for i in range(len(starttimes)):
        start_time = starttimes[i]
        stop_time = start_time + DT.timedelta(minutes=deltat_minutes)

        # Check data availability
        availability = is_data_available(start_time, spans_l1b, nearest_data_after=True)
        if availability['available']:
            stop_time = start_time + DT.timedelta(minutes=deltat_minutes)
            try:
                df_short = read_MATS_data(start_time, stop_time, level=level, version=version)
            except Exception as e:
                print(f"Failed to read MATS data for {start_time} to {stop_time}: {e}")
        else:
            # if nearest data is within half the frequency
            if 'nearest_after' in availability and availability['nearest_after']!=None and availability['nearest_after'] - start_time <= DT.timedelta(days=int(freq_days/2)):
                start_time = availability['nearest_after'] + DT.timedelta(hours=12) #add half a day since to make it the most likely time for the data to start, since nearest_after is in integer days
                stop_time = start_time + DT.timedelta(minutes=deltat_minutes)
                try:
                    df_short = read_MATS_data(start_time, stop_time, level=level, version=version)
                except Exception as e:
                    print(f"Failed to read nearest time MATS data for {start_time} to {stop_time}: {e}")
            else:
                print(f"No data available for Start time: {start_time}, Stop time: {stop_time}")
                if 'nearest_prior' in availability:
                    print(f"Nearest prior data ends at: {availability['nearest_prior']}")
                if 'nearest_after' in availability:
                    print(f"Nearest after data starts at: {availability['nearest_after']}")
                continue


        print(f"Read {len(df_short)} records from {start_time} to {stop_time}")
        # add to main dataframe
        if 'df' not in locals():
            df = df_short
        else:
            df = pd.concat([df, df_short], ignore_index=True)

    if 'df' in locals():
        df.to_pickle(data_folder+name+'.pkl')
    else:
        print("No data was read.")

else:
    df = pd.read_pickle(data_folder+name+'.pkl')

# %%
dflong=df.copy()
#%%
df['darkmean'] = df.apply(lambda x: add_dark_mean(x), axis=1)
df['brightmean'] = df.apply(lambda x: add_bright_mean(x), axis=1)
df['mean'] = df.apply(lambda x: add_mean(x), axis=1)

#%%
df['darkmedian'] = df.apply(lambda x: add_dark_median(x), axis=1)
df['darkcount'] = df.apply(lambda x: add_dark_count(x, threshold=10), axis=1)
#%%
df.to_pickle(data_folder+name+'_withmeanvalues.pkl')
#df=pd.read_pickle(data_folder+name+'_withmeanvalues.pkl')

#%%

df_ascen=df[df['satlat']<df['TPlat']]#ascending node
print(f"Number of records for ascending node: {len(df_ascen)}")
#df=df[df['satlat']>df['TPlat']]#descending node
df_descen=df[df['satlat']>df['TPlat']]#descending node
print(f"Number of records for descending node: {len(df_descen)}")
#%%

def plot_scatter(df_channel, lat_spans, colors, region_label, channel):
    # Scatter plots
    fig_scat, axs_scat = plt.subplots(5, 1, figsize=(14, 18), sharex=True)
    ylabels = ['Dark Mean',  'Dark Median', 'Dark Count','Bright Mean', 'Overall Mean']
    for idx, (latspan, latlabel) in enumerate(lat_spans):
        df_lat = df_channel[(df_channel['satlat'] >= latspan[0]) & (df_channel['satlat'] < latspan[1])].reset_index(drop=True)
        print(f"Number of records for channel {channel} in latitude range {latlabel} ({latspan[0]}° to {latspan[1]}°): {len(df_lat)}")
        if len(df_lat) < 5:
            print(f"Not enough data for channel {channel} in latitude range {latlabel}. Skipping plot.")
            continue
        for i, field in enumerate(['darkmean','darkmedian', 'darkcount','brightmean', 'mean']):
            axs_scat[i].plot(df_lat['TMHeaderTime'], df_lat[field], label=latlabel, color=colors[idx], marker='.', linestyle='None')
    for i, ax in enumerate(axs_scat):
        ax.set_ylabel(ylabels[i])
        ax.set_title(f'{channel} {ylabels[i]} over Time by Latitude ({region_label})')
        ax.legend()
        ax.grid()
    axs_scat[-1].set_xlabel('Date')
    plt.tight_layout()
    plt.show()
    return axs_scat

def plot_binned(df_channel, lat_spans, colors, region_label, freq, freq_days, channel, binned_axes=None):
    # Binned errorbar plots
    fig_err, axs_err = plt.subplots(5, 1, figsize=(14, 18), sharex=True)
    ylabels = ['Dark Mean', 'Dark Median', 'Dark Count', 'Bright Mean', 'Overall Mean']
    for idx, (latspan, latlabel) in enumerate(lat_spans):
        df_lat = df_channel[(df_channel['satlat'] >= latspan[0]) & (df_channel['satlat'] < latspan[1])].reset_index(drop=True)
        if len(df_lat) < 5:
            continue
        df_lat = df_lat.set_index('TMHeaderTime')
        df_binned = df_lat.resample(freq).agg({
            'darkmean': ['mean', 'std', 'count'],
            'darkmedian': ['mean', 'std', 'count'],
            'darkcount': ['mean', 'std', 'count'],
            'brightmean': ['mean', 'std', 'count'],
            'mean': ['mean', 'std', 'count'],
            'satlat': ['mean', 'std', 'count']
        })
        df_binned.columns = ['_'.join(col).strip() for col in df_binned.columns.values]
        df_binned.reset_index(inplace=True)
        #print(f"Binned data for channel {channel} in {latlabel} into {len(df_binned)} time bins of {freq_days} days.")
        for i, field in enumerate(['darkmean', 'darkmedian', 'darkcount','brightmean', 'mean']):
            axs_err[i].errorbar(df_binned['TMHeaderTime'], df_binned[f'{field}_mean'], yerr=df_binned[f'{field}_std'], fmt='o-', label=latlabel, color=colors[idx])
    for i, ax in enumerate(axs_err):
        ax.set_ylabel(ylabels[i])
        ax.set_title(f'{channel} Binned {ylabels[i]} over Time by Latitude ({region_label}, {freq_days}-day bins)')
        ax.legend()
        ax.grid()
    axs_err[-1].set_xlabel('Date')
    plt.tight_layout()
    plt.show()
    if binned_axes is not None:
        for i in range(5):
            binned_axes[i] = axs_err[i]
    return axs_err

def plot_period_comparison(df_channel, lat_spans, colors, period1, period2, labels, channel, binned_axes):
    ylabels = ['Dark Mean', 'Dark Median', 'Dark Count', 'Bright Mean', 'Overall Mean']
    fields = ['darkmean', 'darkmedian', 'darkcount', 'brightmean', 'mean']
    fig, axs = plt.subplots(5, 1, figsize=(14, 18), sharex=True)
    for idx, (latspan, latlabel) in enumerate(lat_spans):
        df_lat = df_channel[(df_channel['satlat'] >= latspan[0]) & (df_channel['satlat'] < latspan[1])].reset_index(drop=True)
        df_period1 = df_lat[(df_lat['TMHeaderTime'] >= period1[0]) & (df_lat['TMHeaderTime'] < period1[1])].reset_index(drop=True)
        df_period2 = df_lat[(df_lat['TMHeaderTime'] >= period2[0]) & (df_lat['TMHeaderTime'] < period2[1])].reset_index(drop=True)
        print(f"{latlabel}: Period 1 records: {len(df_period1)}, Period 2 records: {len(df_period2)}")
        if len(df_period1) < 5 or len(df_period2) < 5:
            print(f"Not enough data for channel {channel} in {latlabel} in one of the periods. Skipping comparison plot.")
            continue
        for i, field in enumerate(fields):
            means = [df_period1[field].mean(), df_period2[field].mean()]
            stds = [df_period1[field].std(), df_period2[field].std()]
            errs = [stds[0]/np.sqrt(len(df_period1)), stds[1]/np.sqrt(len(df_period2))]
            axs[i].errorbar(labels, means, yerr=stds, fmt='o-', color=colors[idx], label=latlabel)
            # Print and annotate change
            mean_diff = means[1] - means[0]
            err_diff = np.sqrt((stds[0]**2/len(df_period1)) + (stds[1]**2/len(df_period2)))
            print(f"{latlabel} {ylabels[i]} Period 1 {channel}: {means[0]:.2f} ± {errs[0]:.2f}, Period 2: {means[1]:.2f} ± {errs[1]:.2f}")
            print(f"{latlabel} {ylabels[i]} change: {mean_diff:.2f} ± {err_diff:.2f}")
            if binned_axes is not None:
                binned_axes[i].text(0.5, 0.9-idx*0.1, f"{latlabel} {ylabels[i]} change: {mean_diff:.2f} ± {err_diff:.2f}", transform=binned_axes[i].transAxes, color=colors[idx])
    for i, ax in enumerate(axs):
        ax.set_ylabel(ylabels[i])
        ax.set_title(f'{channel} {ylabels[i]} Comparison by Latitude (NH)')
        ax.legend()
        ax.grid()
    axs[-1].set_xlabel('Period')
    plt.tight_layout()
    plt.show()

direction = 'both'  # 'ascend', 'descend' or 'both'
for channel in ['UV1', 'UV2', 'IR1', 'IR2']:
    if direction == 'ascend':
        df_channel = df_ascen[df_ascen['channel'] == channel].reset_index(drop=True)
    elif direction == 'descend':
        df_channel = df_descen[df_descen['channel'] == channel].reset_index(drop=True)
    else:
        df_channel = df[df['channel'] == channel].reset_index(drop=True)

    lat_spans_north = [
        ((70, 80), 'NH High Latitude'),
        ((50, 60), 'NH Nominal'),
        ((20, 40), 'NH Low Latitude')
    ]
    lat_spans_south = [
        ((-40, -20), 'SH Low Latitude'),
        ((-60, -50), 'SH Nominal'),
        ((-80, -70), 'SH High Latitude')
    ]
    colors_north = ['tab:blue', 'tab:red', 'tab:green']
    colors_south = ['tab:purple', 'tab:brown', 'tab:gray']

    # North + Equator
    axsNerr = [None]*5
    #plot_scatter(df_channel, lat_spans_north, colors_north, 'NH + Equator', channel)
    plot_binned(df_channel, lat_spans_north, colors_north, 'NH + Equator', freq, freq_days, channel, binned_axes=axsNerr)
    # South
    axsSerr = [None]*5
    #plot_scatter(df_channel, lat_spans_south, colors_south, 'SH', channel)
    #plot_binned(df_channel, lat_spans_south, colors_south, 'SH', freq, freq_days, channel, binned_axes=axsSerr)

    # Period comparison for NH
    period1 = (DT.datetime(2023, 12, 1, tzinfo=DT.timezone.utc), DT.datetime(2024, 1, 31, tzinfo=DT.timezone.utc))
    period2 = (DT.datetime(2024, 12, 1, tzinfo=DT.timezone.utc), DT.datetime(2025, 1, 31, tzinfo=DT.timezone.utc))
    labels = ['2023-12 to 2024-01', '2024-12 to 2025-01']
    plot_period_comparison(df_channel, lat_spans_north, colors_north, period1, period2, labels, channel, axsNerr)


# %%
for channel in ['UV1','UV2','IR1','IR2']:
    df_channel = df[df['channel'] == channel].reset_index(drop=True)
    plt.figure(figsize=(14, 6))
    plt.hist(df_channel.iloc[0].ImageCalibrated.flatten(), bins=100, color='blue', alpha=0.7)
    plt.title(f'{channel} ImageCalibrated Histogram ')
    plt.xlabel('Pixel Value')
    plt.ylabel('Frequency (log scale)')

# %%
