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
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=120000)
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
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=120000)
        image = CCDitem.ImageCalibrated.copy()

        # Set pixels below threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[:threshold_row, col] = np.nan
        dark_median = np.nanmedian(image)
    except Exception as e:
        print(f"Error in add_dark_median: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        dark_median = np.nan
    return dark_median
def add_dark_count(CCDitem, threshold=40):
    try:
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=120000)
        image = CCDitem.ImageCalibrated.copy()

        # Set pixels below threshold altitude to NaN
        for col, threshold_row in enumerate(pixel_indices):
            image[:threshold_row, col] = np.nan
        dark_count = np.nansum(image > threshold)/np.nansum(~np.isnan(image))
    except Exception as e:
        print(f"Error in add_dark_count: {e}, channel: {getattr(CCDitem, 'channel', 'Unknown')}")
        dark_count = np.nan
    return dark_count
def add_bright_mean(CCDitem):
    try:
        pixel_indices = compute_threshold_pixels(CCDitem, thresholdalt=95000)
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


def plot_scatter(df_channel, lat_spans, colors, region_label, channel):
    # Scatter plots
    plt.rcParams.update({'font.size': 16})  # Increase font size
    fig_scat, axs_scat = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
    unittext = ' Brightness in 10^12 photons/m^2/s/nm/sr'
    ylabels = [
        'Dark Mean',
        'Dark Median',
        'Dark Count'
    ]
    fields = ['darkmean', 'darkmedian', 'darkcount']
    for idx, (latspan, latlabel) in enumerate(lat_spans):
        df_lat = df_channel[(df_channel['satlat'] >= latspan[0]) & (df_channel['satlat'] < latspan[1])].reset_index(drop=True)
        print(f"Number of records for channel {channel} in latitude range {latlabel} ({latspan[0]}° to {latspan[1]}°): {len(df_lat)}")
        if len(df_lat) < 5:
            print(f"Not enough data for channel {channel} in latitude range {latlabel}. Skipping plot.")
            continue
        for i, field in enumerate(fields):
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
    plt.rcParams.update({'font.size': 16})  # Increase font size
    fig_err, axs_err = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
    unittext = ' Brightness in 10^12 photons/m^2/s/nm/sr'
    ylabels = [
        'Dark Mean',
        'Dark Median',
        'Dark Count'
    ]
    fields = ['darkmean', 'darkmedian', 'darkcount']
    for idx, (latspan, latlabel) in enumerate(lat_spans):
        df_lat = df_channel[(df_channel['satlat'] >= latspan[0]) & (df_channel['satlat'] < latspan[1])].reset_index(drop=True)
        if len(df_lat) < 5:
            continue
        df_lat = df_lat.set_index('TMHeaderTime')
        df_binned = df_lat.resample(freq).agg({
            'darkmean': ['mean', 'std', 'count'],
            'darkmedian': ['mean', 'std', 'count'],
            'darkcount': ['mean', 'std', 'count'],
            'satlat': ['mean', 'std', 'count']
        })
        df_binned.columns = ['_'.join(col).strip() for col in df_binned.columns.values]
        df_binned.reset_index(inplace=True)
        for i, field in enumerate(fields):
            axs_err[i].errorbar(
                df_binned['TMHeaderTime'],
                df_binned[f'{field}_mean'],
                yerr=df_binned[f'{field}_std'],
                fmt='o-',
                label=latlabel,
                color=colors[idx]
            )
    for i, ax in enumerate(axs_err):
        ax.set_ylabel(ylabels[i])
        ax.set_title(f'{channel} Binned {ylabels[i]} over Time by Latitude ({region_label}, {freq_days}-day bins)')
        ax.legend()
        ax.grid()
    axs_err[-1].set_xlabel('Date')
    plt.tight_layout()
    plt.show()
    if binned_axes is not None:
        for i in range(3):
            binned_axes[i] = axs_err[i]
    return axs_err


def compute_change(df_channel, lat_spans, colors, period1, period2, labels, channel, binned_axes=None, plot=False):
    unittext=' Brightness in 10^12 photons/m^2/s/nm/sr'
    ylabels = ['Dark Mean'+unittext, 'Dark Median'+unittext, 'Dark Count', 'Bright Mean'+unittext]
    fields = ['darkmean', 'darkmedian', 'darkcount', 'brightmean']
    results = []
    if plot:
        fig, axs = plt.subplots(5, 1, figsize=(14, 18), sharex=True)
    for idx, (latspan, latlabel) in enumerate(lat_spans):
        df_lat = df_channel[(df_channel['satlat'] >= latspan[0]) & (df_channel['satlat'] < latspan[1])].reset_index(drop=True)
        df_period1 = df_lat[(df_lat['TMHeaderTime'] >= period1[0]) & (df_lat['TMHeaderTime'] < period1[1])].reset_index(drop=True)
        df_period2 = df_lat[(df_lat['TMHeaderTime'] >= period2[0]) & (df_lat['TMHeaderTime'] < period2[1])].reset_index(drop=True)
        entry = {
            'latlabel': latlabel,
            'latspan': latspan,
            'period1_n': len(df_period1),
            'period2_n': len(df_period2),
            'fields': {}
        }
        if len(df_period1) < 5 or len(df_period2) < 5:
            entry['insufficient_data'] = True
            results.append(entry)
            continue
        entry['insufficient_data'] = False
        for i, field in enumerate(fields):
            mean1 = df_period1[field].mean()
            mean2 = df_period2[field].mean()
            std1 = df_period1[field].std()
            std2 = df_period2[field].std()
            n1 = len(df_period1)
            n2 = len(df_period2)
            err1 = std1 / np.sqrt(n1) if n1 > 0 else np.nan
            err2 = std2 / np.sqrt(n2) if n2 > 0 else np.nan
            mean_diff = mean2 - mean1
            err_diff = np.sqrt((std1**2/n1) + (std2**2/n2)) if n1 > 0 and n2 > 0 else np.nan
            entry['fields'][field] = {
                'mean1': mean1,
                'mean2': mean2,
                'std1': std1,
                'std2': std2,
                'err1': err1,
                'err2': err2,
                'mean_diff': mean_diff,
                'err_diff': err_diff
            }
            if plot:
                axs[i].errorbar(labels, [mean1, mean2], yerr=[std1, std2], fmt='o-', color=colors[idx], label=latlabel if i == 0 else "")
                if binned_axes is not None:
                    binned_axes[i].text(0.5, 0.9-idx*0.1, f"{latlabel} {ylabels[i]} change: {mean_diff:.2f} ± {err_diff:.2f}", transform=binned_axes[i].transAxes, color=colors[idx])
        results.append(entry)
    if plot:
        for i, ax in enumerate(axs):
            ax.set_ylabel(ylabels[i])
            ax.set_title(f'{channel} {ylabels[i]} Comparison by Latitude')
            ax.legend()
            ax.grid()
        axs[-1].set_xlabel('Period')
        plt.tight_layout()
        plt.show()
    return results

# --- New function for histogram comparison (dark region only) ---
def stack_dark_images(df,thresholdalt=130000,  compare_to_darkest=False):
    imgs = []
    for _, x in df.iterrows():
        image = x.ImageCalibrated.copy()
        if compare_to_darkest:
            # find the darkest percentile value in the image
            darkest_value = np.nanpercentile(image, 1)  # 1st percentile
            # Subtract this value from the entire image
            image = image - darkest_value
        pixel_indices = compute_threshold_pixels(x, thresholdalt=thresholdalt)
        # Mask bright part (below thresholdalt) as NaN
        mask = np.zeros_like(image, dtype=bool)
        for col, threshold_row in enumerate(pixel_indices):
            mask[:threshold_row, col] = True
        dark_part = image.copy()
        dark_part[mask] = np.nan
        imgs.append(dark_part)
    if len(imgs) == 0:
        return None
    return imgs

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
#pick out only channels IR1, IR2, UV1, UV2
df=dflong[dflong['channel'].isin(['IR1', 'IR2', 'UV1', 'UV2'])].reset_index(drop=True)
df['darkmean'] = df.apply(lambda x: add_dark_mean(x), axis=1)
print('Darkmean done')
df['brightmean'] = df.apply(lambda x: add_bright_mean(x), axis=1)
print('Brightmean done')
df['darkmedian'] = df.apply(lambda x: add_dark_median(x), axis=1)
print('Darkmedian done')
df['darkcount'] = df.apply(lambda x: add_dark_count(x, threshold=40), axis=1)
print('Darkcount done')
#%%
#df.to_pickle(data_folder+name+'_withmeanvalues.pkl')
df=pd.read_pickle(data_folder+name+'_withmeanvalues.pkl')



#%%

df_ascen=df[df['satlat']<df['TPlat']]#ascending node
print(f"Number of records for ascending node: {len(df_ascen)}")
#df=df[df['satlat']>df['TPlat']]#descending node
df_descen=df[df['satlat']>df['TPlat']]#descending node
print(f"Number of records for descending node: {len(df_descen)}")
#%%


#%%
direction = 'both'  # 'ascend', 'descend' or 'both'
for channel in ['IR1']:
    if direction == 'ascend':
        df_channel = df_ascen[df_ascen['channel'] == channel].reset_index(drop=True)
    elif direction == 'descend':
        df_channel = df_descen[df_descen['channel'] == channel].reset_index(drop=True)
    else:
        df_channel = df[df['channel'] == channel].reset_index(drop=True)
    #remove flagged images

    print(f"Number of records for channel {channel} before removing flagged images: {len(df_channel)}")
    df_channel, mask=remove_flagged_images(df_channel, [8, 9], return_mask=True) # Remove images where desmearing malfunctioned
    print(f"Number of records for channel {channel} after removing flagged images: {len(df_channel)}")
    # Define latitude spans and colors

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

    # North 
    axsNerr = [None]*5
    plot_scatter(df_channel, lat_spans_north, colors_north, 'NH', channel)
    plot_binned(df_channel, lat_spans_north, colors_north, 'NH', freq, freq_days, channel, binned_axes=axsNerr)
    # South
    axsSerr = [None]*5
    plot_scatter(df_channel, lat_spans_south, colors_south, 'SH', channel)
    plot_binned(df_channel, lat_spans_south, colors_south, 'SH', freq, freq_days, channel, binned_axes=axsSerr)

    # Period comparison for NH
    period1 = (DT.datetime(2023, 12, 1, tzinfo=DT.timezone.utc), DT.datetime(2024, 1, 31, tzinfo=DT.timezone.utc))
    period2 = (DT.datetime(2024, 12, 1, tzinfo=DT.timezone.utc), DT.datetime(2025, 1, 31, tzinfo=DT.timezone.utc))
    labels = ['2023-12 to 2024-01', '2024-12 to 2025-01']
    change_results = compute_change(df_channel, lat_spans_north, colors_north, period1, period2, labels, channel, axsNerr, plot=False)
    # To plot, set plot=True in the call above
#%%


#%%

for channel in ['IR1', 'IR2', 'UV1', 'UV2']:
    #channel = 'IR1'
    df_channel = df[(df['channel'] == channel) & (df['TPheight'] > 110000)].reset_index(drop=True)
    #Remove flagged images
    df_channel, mask=remove_flagged_images(df_channel, [8, 9], return_mask=True) # Remove images where desmearing malfunctioned

    most_common_shape = df_channel.ImageCalibrated.apply(lambda x: x.shape).mode()[0]
    df_channel = df_channel[df_channel.ImageCalibrated.apply(lambda x: x.shape) == most_common_shape].reset_index(drop=True)
    df_p1 = df_channel[(df_channel['TMHeaderTime'] >= period1[0]) & (df_channel['TMHeaderTime'] < period1[1])]
    df_p2 = df_channel[(df_channel['TMHeaderTime'] >= period2[0]) & (df_channel['TMHeaderTime'] < period2[1])]
    print(f"Number of {channel} records in Period 1: {len(df_p1)}, Period 2: {len(df_p2)}")


    # Stack and compute mean images for both periods
    images_stack_p1 = stack_dark_images(df_p1, thresholdalt=120000, compare_to_darkest=False)
    images_stack_p2 = stack_dark_images(df_p2, thresholdalt=120000, compare_to_darkest=False)
    median_img_p1 = np.nanmedian(images_stack_p1, axis=0) if images_stack_p1 is not None else None
    median_img_p2 = np.nanmedian(images_stack_p2, axis=0) if images_stack_p2 is not None else None


    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2], width_ratios=[1, 1], hspace=0.3, wspace=0.15)

    # Top row: Mean images side by side
    ax1 = fig.add_subplot(gs[0, 0])
    if channel in ['IR1', 'IR2', 'IR3', 'IR4']:
        clim = (0, 50)
    else:
        clim = (0, 200)
    if median_img_p1 is not None:
        plot_CCDimage(median_img_p1, fig=fig, axis=ax1, cmap='inferno', clim=clim, title=f'{channel} Median ImageCalibrated (Dark region)\nPeriod 1')
    else:
        ax1.set_title('No data for Period 1')
        ax1.axis('off')

    ax2 = fig.add_subplot(gs[0, 1])
    if median_img_p2 is not None:
        plot_CCDimage(median_img_p2, fig=fig, axis=ax2, cmap='inferno', clim=clim, title=f'{channel} Median ImageCalibrated (Dark region)\nPeriod 2')
    else:
        ax2.set_title('No data for Period 2')
        ax2.axis('off')

    # Bottom row: Normalized histogram comparison
    ax3 = fig.add_subplot(gs[1, :])
    if channel in ['IR1', 'IR2', 'IR3', 'IR4']:
        bins = np.linspace(0, 100, 101)
        maxxlim=100
    else:
        bins = np.linspace(0, 200, 201)
        maxxlim=200
    if (images_stack_p1 is not None) and (images_stack_p2 is not None):
        ax3.hist(np.array(images_stack_p1).flatten(), bins=bins, color='tab:blue', alpha=0.6, label='Period 1', density=True)
        ax3.hist(np.array(images_stack_p2).flatten(), bins=bins, color='tab:orange', alpha=0.6, label='Period 2', density=True)
    elif images_stack_p1 is not None:
        ax3.hist(np.array(images_stack_p1).flatten(), bins=bins, color='tab:blue', alpha=0.6, label='Period 1', density=True)
    elif images_stack_p2 is not None:
        ax3.hist(np.array(images_stack_p2).flatten(), bins=bins, color='tab:orange', alpha=0.6, label='Period 2', density=True)
    ax3.set_xlim(0, maxxlim)
    ax3.set_title(f'{channel} Histogram Comparison (Dark region)')
    ax3.set_xlabel('Pixel Value')
    ax3.set_ylabel('Normalized Frequency')
    ax3.legend()
    ax3.grid()

    plt.tight_layout()
    plt.show()

#%%


# %%
from mats_utils.geolocation.coordinates import col_heights
import numpy as np
def compute_threshold_pixels(CCDitem, thresholdalt):
    """
    Calculates what pixel (ie what row) that is at the threshold altitude in each column

    Args:
        CCDitem: row of a dataframe with CCDitems
        thresholdalt: float with the altitude threshold

    Returns:
        pixel_indices: list of pixel indices (rows) at the threshold altitude in each column
    """
    thl = col_heights(CCDitem, 0, 2) # (left) the 2 makes the function return the lowest and highest pixel
    thr = col_heights(CCDitem, CCDitem.NCOL, 2) # (right) the 2 makes the function return the lowest and highest pixel
    altrowspan = thl[1] - thl[0]  # Altitude span of the image, could equally well have used thr
    altcoldiff = thl[0] - thr[0]  # Altitude difference between the sides of the image
    pixel_indices = []
    if thl[1] < thresholdalt and thr[1] < thresholdalt: #satellite pointing is too low for all pixels
        # no pixels should be selected
        pixel_indices = np.full(CCDitem.NCOL+1, CCDitem.NROW)
    #elif thl[1] < thresholdalt or thr[1] < thresholdalt: #satellite pointing may be ok for some images but we are conservative
        # no pixels should be selected
    #    pixel_indices = np.full(CCDitem.NCOL+1, CCDitem.NROW)
    elif thl[0] > thresholdalt and thr[0] > thresholdalt:# satellite pointing is high enough for all pixels to be dark, ie selected
        #set all pixels to 0
        pixel_indices = np.full(CCDitem.NCOL+1, 0)
    else:  # satellite pointing is ok for some pixels in each column
        rows=CCDitem.NROW
        dalt_dcol=altcoldiff / (CCDitem.NCOL - 1)  # Altitude difference per column
        dalt_drow=altrowspan / (rows - 1)  # Altitude difference per row
        drow_dcol=dalt_dcol / dalt_drow  # Number of rows per column
        # Calculate the pixel (ie row ) that corresponds to the threshold altitude in the first column
        firstpixel_index = int((thresholdalt - thl[0]) / altrowspan* (rows - 1))  # Pixel row index in the first column

        for col in range(CCDitem.NCOL+1):        
            pixel_index = int(firstpixel_index + col * drow_dcol)  # Calculate the pixel index for each column
            if pixel_index < 0:
                pixel_index = 0
            elif pixel_index >= rows:
                pixel_index = rows - 1
            pixel_indices.append(pixel_index)
    return pixel_indices
# %%
