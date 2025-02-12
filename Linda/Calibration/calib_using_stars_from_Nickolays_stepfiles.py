
#%%
import pandas as pd
import glob
import re
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm


def fullframe_subpix(channel,crop, column,row):
    """
    Calculate the subpixel position of the column and row in the fullframe image

    Parameters
    ----------
    channel : str
        The channel of the CCD
    crop : str
        The crop version of the CCD
    column : int
        The column of the cropped and binned image
    row : int
        The row of the cropped and binned image

    Returns
    -------
    fullframesubcolumns : nparray
        The subpixel columns of the original image that was binned
    fullframesubrows : nparray
        The subpixel rows of the original image that was binned
    """

    from mats_utils.rawdata.cropping import set_crop_settings
    import numpy as np

    NRSKIP, NRBIN, NROW, NCSKIP, NCBIN, NCOL, NCBINFPGA=set_crop_settings(channel, crop)
    # Calculate the subpixel position of the column and row in the fullframe image
    fullframe_lastsubcolumn = NCSKIP+column*NCBIN*NCBINFPGA
    fullframe_firstsubcolumn = fullframe_lastsubcolumn-NCBIN*NCBINFPGA+1    
    fullframesubcolumns = np.arange(fullframe_firstsubcolumn, fullframe_lastsubcolumn+1)

    fullframe_lastsubrow = NRSKIP+row*NRBIN 
    fullframe_firstsubrow = fullframe_lastsubrow-NRBIN+1
    fullframesubrows = np.arange(fullframe_firstsubrow, fullframe_lastsubrow+1)

    return fullframesubcolumns, fullframesubrows

def compute_avarage_flatfield_factor(ff_columns, ff_rows, channel, lenthofstreak=1):
    """
    Compute the average flatfield factor for the given columns and rows

    Parameters
    ----------
    ff_columns : nparray
        The subpixel columns of the original image that was binned
    ff_rows : nparray
        The subpixel rows of the original image that was binned
    channel : str
        The channel of the CCD
    lenthofstreak : int
        The length of the star streak in the number of binned rows

    Returns
    -------
    flatfieldfact : float
        The average flatfield factor for the given columns and rows
    """


    flatfield = np.load('/Users/lindamegner/MATS/MATS-retrieval/MATS-instrument-data/calibration_data/flatfields/flatfield_'+channel+'_HSM.npy')
    
    #allow for computation of the full streak
    if lenthofstreak!=1:
        # first=int(ff_rows.mean()-lenthofstreak*len(ff_rows)/2)
        # if first<0:
        #     first=0
        # last=int(ff_rows.mean()+lenthofstreak*len(ff_rows)/2)
        # if last>flatfield.shape[0]:
        #     last=flatfield.shape[0]

        first=ff_rows[0]
        last=int(first+lenthofstreak*len(ff_rows))
        if last>flatfield.shape[0]:
            last=flatfield.shape[0]

    
    flatfieldfact = np.mean(flatfield[np.ix_(ff_rows, ff_columns)])

    return flatfieldfact

def fitgaussian(df, column, plot=False, plotxlim=None):
    """
    Fit a Gaussian to the histogram data

    Parameters
    ----------
    df : DataFrame
        The DataFrame containing the data
    column : str
        The column containing the data to fit

    Returns
    -------
    mu : float
        The mean of the Gaussian fit
    sigma : float
        The standard deviation of the Gaussian fit
    """

    from scipy.stats import norm, kstest

    # Fit a Gaussian to the histogram data
    (mu, sigma) = norm.fit(df[column])
    (mu, sigma) = norm.fit(df[df[column].between(mu-3*sigma, mu+3*sigma)][column])
    (mu, sigma) = norm.fit(df[df[column].between(mu-3*sigma, mu+3*sigma)][column])
    (mu, sigma) = norm.fit(df[df[column].between(mu-2*sigma, mu+2*sigma)][column])

    # Calculate the goodness of fit using the Kolmogorov-Smirnov test
    ks_stat, p_value = kstest(df[df[column].between(mu-2*sigma, mu+2*sigma)][column], 'norm', args=(mu, sigma))
    
    if plot:
        fig, ax = plt.subplots()
        n, bins, patches = ax.hist(df[df[column].between(mu-4*sigma, mu+4*sigma)][column], bins=50, color='r', alpha=0.2, label=column)

        # Plot the Gaussian fit
        x = np.linspace(bins[0], bins[-1], 500)
        y = norm.pdf(x, mu, sigma) * len(df[column]) * (bins[1] - bins[0])
        ax.plot(x, y, 'k--', linewidth=2, label=f'Gaussian fit: $\mu={mu:.2f}$, $\sigma={sigma:.2f}$')
        ax.text(0.1, 0.9, f'KS test: $p={p_value:.2f}$', transform=ax.transAxes)
        ax.text(0.1, 0.8, str(ks_stat), transform=ax.transAxes)
        ax.text(0.1, 0.7, f'$\mu={mu:.2f}$, $\sigma={sigma:.2f}$', transform=ax.transAxes)
        if plotxlim:
            ax.set_xlim(plotxlim)

        ax.set_ylabel(column)
        ax.set_title(df.iloc[0].channel+' star nr: '+str(df.iloc[0].starnumber))

    return mu, sigma, ks_stat, p_value

# Define the function to fit a line through the origin and return slope with error estimates
def fit_line_through_origin(x, y, xerr=None, yerr=None):
    if xerr is None:
        xerr = np.ones_like(x)
    if yerr is None:
        yerr = np.ones_like(y)
    
    # Calculate weights as the inverse of the errors
    weights = 1 / (xerr**2 + yerr**2)
    
    # Calculate the weighted slope
    slope = np.sum(weights * x * y) / np.sum(weights * x**2)
    
    # Calculate the residuals
    residuals = y - slope * x
    
    # Calculate the standard error of the slope
    std_err = np.sqrt(np.sum(weights * residuals**2) / (len(x) - 1)) / np.sqrt(np.sum(weights * x**2))
    
    return slope, std_err
# # Define the function to fit a line through the origin and return slope with error estimates
# def fit_line_through_origin(x, y):
#     # Calculate the slope directly
#     slope = np.sum(x * y) / np.sum(x**2)
    
#     # Calculate the residuals
#     residuals = y - slope * x
    
#     # Calculate the standard error of the slope
#     std_err = np.sqrt(np.sum(residuals**2) / (len(x) - 1)) / np.sqrt(np.sum(x**2))
    
#     return slope, std_err

#%%

nickodatabasepath='/Users/lindamegner/MATS/MATS-retrieval/data/NickolaysStars/'
# Use glob to get all the files in the directory
all_files = glob.glob(nickodatabasepath +"star_analysis_steps/star_*_step2_.txt")

# Create an empty list to store the dataframes
dfs = []

# Loop through the files and read them into dataframes
for file in all_files:
    # Extract starnumber and channel from the filename
    match = re.search(r'star_(\d+)_([A-Z0-9]+)_step2_', file)
    if match:
        starnumber = int(match.group(1))
        channel = match.group(2)
    
    # Read the file into a dataframe
    df = pd.read_csv(file, delim_whitespace=True, header=None, names=['date', 'exposure_number', 'column', 'row', 'corrFitCoeff', 'count'])
    
    # Add the starnumber and channel columns
    df['starnumber'] = starnumber
    df['channel'] = channel
    
    # Append the dataframe to the list
    dfs.append(df)

# Concatenate all the dataframes into a single dataframe

allstarsdf = pd.concat(dfs, ignore_index=True)
# Convert Date column to datetime object
allstarsdf['date'] = pd.to_datetime(allstarsdf['date'], format='%Y%m%d')

# Display the first few rows of the dataframe
print(allstarsdf.head())


# Define the path to the Excel file
file_path =nickodatabasepath+ 'all_stars_20240513lm.xlsx'

# Read the Excel file into a DataFrame
starsdf = pd.read_excel(file_path, sheet_name='Allstars')
starsdf = starsdf.dropna(subset=['flux', 'count', 'name'])
print(starsdf)

allstarsdf = allstarsdf.merge(starsdf[['starnumber', 'name', 'channel', 'flux']], on=['starnumber', 'channel'], how='left')

#select date for crop d %Crop D from Feb 09 18.00 to May 3 24.00
crop_d_start = pd.Timestamp('2023-02-09 18:00:00')
crop_d_end = pd.Timestamp('2023-05-04 00:00:00')
allstarsdf = allstarsdf[(allstarsdf['date'] >= crop_d_start) ]
allstarsdf['crop']=allstarsdf['date'].apply(lambda x: 'CROPD' if x<=crop_d_end else 'CROPF')
#& (allstarsdf['date'] <= crop_d_end)]
# select only the rows, columns within the region free from  baffle interference
#ddd = (abs(column(:,3))<3) & (abs(column(:,4))<10);
#allstarsdf = allstarsdf[(abs(allstarsdf['column']) < 3) & (abs(allstarsdf['row']) < 10)]


#%%

# Add the inside_FOV column based on the conditions provided
#BASED on Nickolays step3 files
#conditions_IR3_IR4 = (abs(allstarsdf['column']) < 3) & (abs(allstarsdf['row']) < 10) & (allstarsdf['channel'].isin(['IR3', 'IR4']))
#conditions_IR1_IR2_UV1_UV2 = (abs(allstarsdf['column']) < 150) & (abs(allstarsdf['row']) < 200) & (allstarsdf['channel'].isin(['IR1', 'IR2', 'UV1', 'UV2']))
#allstarsdf['inside_FOV'] = conditions_IR3_IR4 | conditions_IR1_IR2_UV1_UV2


allstarsdf[['ff_columns', 'ff_rows']] = allstarsdf.apply(lambda row: fullframe_subpix(row['channel'], row['crop'], row['column'], row['row']), axis=1, result_type='expand')

FirstRow = 350 
LastRow = 400 
FirstCol =524
LastCol =1523
allstarsdf['inside_FOV'] = allstarsdf.apply(lambda row: (row['ff_rows'][int(len(row['ff_rows'])/2)] > FirstRow) & (row['ff_rows'][int(len(row['ff_rows'])/2)] < LastRow) & (row['ff_columns'][int(len(row['ff_columns'])/2)] > FirstCol) & (row['ff_columns'][int(len(row['ff_columns'])/2)] < LastCol), axis=1) 
#allstarsdf['inside_FOV'].sum()
#allstarsdf = allstarsdf[allstarsdf['inside_FOV'] == True]

#%%

#read in flatfields in a dictionary
#create a dictionary with the flatfields
# flatfield = {}
# for channel in allstarsdf['channel'].unique():
#      flatfield[channel] = np.load('/Users/lindamegner/MATS/MATS-retrieval/MATS-instrument-data/calibration_data/flatfields/flatfield_'+channel+'_HSM.npy')
#      plt.imshow(flatfield[channel])
#      plt.title(channel)
#      plt.colorbar()
#      plt.show()





allstarsdf['exptime'] = allstarsdf['channel'].map({'IR1': 5, 'IR2': 3, 'IR3': 1.5, 'IR4': 1.5, 'UV1': 5, 'UV2': 5})
#allstarsdf['exptime_nickolay'] = allstarsdf['Cal']*allstarsdf['count']/allstarsdf['flux'] 

#set lengthofstreak= 9 for IR1, IR2, UV1, UV2 and 3 for IR3, IR4
# As given from Nikcolay:
# %leng=50;  % IR1
# % leng=32;  % IR2
# % leng=6;  % IR3
# leng=6;  % IR4
#This gives 10 binned rows per second, ie 20 unbinned rows per second
allstarsdf['lengthofstreak']=22*allstarsdf['exptime']/allstarsdf['ff_rows'].apply(lambda x: len(x))
#allstarsdf['lengthofstreak'] = allstarsdf['channel'].map({'IR1': 50, 'IR2': 30, 'IR3': 3, 'IR4': 3, 'UV1': 9, 'UV2': 9})

lincomp=False
if lincomp==True:
    print('**********Compensate for linearity**********')
    #Compensate for nonlinearity
    from mats_l1_processing.instrument import Instrument
    global instrument
    instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')

    def compensate_for_nonlinearity(count, lengthofstreak, channel):
        """
        Compensate for nonlinearity in the count data

        Parameters
        ----------
        count : float
            The count data to compensate
        channel : str
            The channel of the CCD

        Returns
        -------
        float
            The compensated count data
        """
        # Get the CCD item for the channel
        #CCDitem = CCDitemdict[channel]
        CCDunit=instrument.get_CCD(channel)
        
        # Compensate for nonlinearity
        count = CCDunit.non_linearity.get_true_value(count/lengthofstreak)*lengthofstreak
        
        return count
    allstarsdf['count_noncomp'] = allstarsdf['count'].copy()
    allstarsdf['count'] =allstarsdf.apply(lambda row: compensate_for_nonlinearity(row['count'], row['lengthofstreak'], row['channel']), axis=1)
else:
    print('**********Not compensating for linearity**********')
    

allstarsdf['flatfieldfact'] = allstarsdf.apply(lambda row: compute_avarage_flatfield_factor(row['ff_columns'], row['ff_rows'], row['channel'], lenthofstreak=row['lengthofstreak']), axis=1)




#%%


allstarsdf['count_corrected'] = allstarsdf['count']/allstarsdf['flatfieldfact']
allstarsdf['calcoeff']=allstarsdf['flux']/allstarsdf['count']*allstarsdf['exptime']
allstarsdf['calcoeff_corrected']=allstarsdf['flux']/allstarsdf['count_corrected']*allstarsdf['exptime']


#%%
# plot the calibration coefficients for each channel, as a function of the star number
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
allstarsdf = allstarsdf.dropna(subset=['calcoeff_corrected'])

allstarsdf['outlier'] = False
figdict = {}
axdict = {}
calcoeffmeandict= {}
calcoeffstddict= {}
calcoeffmeanerrordict= {}

plotxlimdir = {'IR1': [8, 12], 'IR2': [2.7, 3.1], 'IR3': [16, 26], 'IR4': [17, 37], 'UV1': [45, 65], 'UV2': [16, 28]}

for channel in allstarsdf['channel'].unique():
    df = allstarsdf[allstarsdf['channel'] == channel]

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    figdict[channel] = fig
    axdict[channel] = ax

    # plot with title set to channel
    ax[0,0].plot(df['starnumber'], df['calcoeff'], 'o', color='b', alpha=0.2, label='not corrected')
    ax[0,0].plot(df['starnumber'], df['calcoeff_corrected'], 'o', color='r', alpha=0.2, label='flatfield corrected')
    #set limis on y-axis to mean + 3*std
    meancalcoeff = df['calcoeff_corrected'].mean()
    stdcalcoeff = df['calcoeff_corrected'].std()
    ax[0,0].set_ylim(meancalcoeff-3*stdcalcoeff, meancalcoeff+3*stdcalcoeff)
    ax[0,0].set_title(channel)
    ax[0,0].set_xlabel('Star number')
    ax[0,0].set_ylabel('Calibration coefficient')
    ax[0,0].legend()




    # plot the histogram of the calibration coefficients
    #fig, ax = plt.subplots()
    if channel in ['IR1', 'IR2', 'IR3', 'IR3']:
        nhistbins=500
    else:
        nhistbins=30
    ax[0,1].hist(df['calcoeff'], bins=nhistbins, color='b', alpha=0.2, label='not corrected')
    n, bins, patches = ax[0,1].hist(df['calcoeff_corrected'], bins=nhistbins, color='r', alpha=0.2, label='flatfield corrected')

    # Fit a Gaussian to the histogram data
    (mu, sigma) = norm.fit(df['calcoeff_corrected'])

    # Plot the Gaussian fit
    x = np.linspace(bins[0], bins[-1], 500)
    y = norm.pdf(x, mu, sigma) * len(df['calcoeff_corrected']) * (bins[1] - bins[0])
    ax[0,1].plot(x, y, 'k--', linewidth=2, label=f'Gaussian fit: $\mu={mu:.2f}$, $\sigma={sigma:.2f}$')

    #remove data outside of 3 sigma
    #n, bins, patches = ax[0,1].hist(df[df['calcoeff_corrected'].between(mu-3*sigma, mu+3*sigma)]['calcoeff_corrected'], bins=30, color='r', alpha=0.2, label='flatfield corrected')
    # Fit a Gaussian to the histogram data
    (mu_new1, sigma_new1) = norm.fit(df[df['calcoeff_corrected'].between(mu-0.6*mu, mu+0.6*mu)]['calcoeff_corrected'])
    #Fit again
    (mu_new, sigma_new) = norm.fit(df[df['calcoeff_corrected'].between(mu_new1-2*sigma_new1, mu_new1+2*sigma_new1)]['calcoeff_corrected'])

    # Plot the Gaussian fit
    y = norm.pdf(x, mu_new, sigma_new) * len(df['calcoeff_corrected']) * (bins[1] - bins[0])
    #error of mean
    err_of_mean = sigma_new / np.sqrt(len(df[df['calcoeff_corrected'].between(mu_new1-2*sigma_new1, mu_new1+2*sigma_new1)]['calcoeff_corrected'])-1)
    err_of_mean_percent = err_of_mean/mu_new*100
    calcoeffmeandict[channel] = mu_new
    calcoeffstddict[channel] = sigma_new
    calcoeffmeanerrordict[channel] = err_of_mean
    ax[0,1].plot(x, y, 'k--', linewidth=2, label=f'Gaussian fit without outliers: $\mu={mu_new:.2f}$, $\sigma={sigma_new:.2f}$', color='g')
    ax[0,1].set_title(f'Mean: {mu_new:.2f}, Error of mean: {err_of_mean:.2f}, Error of mean in %: {err_of_mean_percent:.2f}')
    ax[0,1].set_xlim(mu_new1-3*sigma_new1, mu_new1+3*sigma_new1)
    # Add labels and legend
    ax[0,1].set_xlabel('calcoeff_corrected')
    ax[0,1].set_ylabel('Frequency')
    ax[0,1].legend()

    #mark outliers larger than 3 sigma away from the mean
    allstarsdf[allstarsdf['channel'] == channel]['outlier'] = allstarsdf[allstarsdf['channel'] == channel]['calcoeff_corrected'].between(mu_new-2*sigma, mu_new+2*sigma)


    #fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    scatter = ax[1,0].scatter(1/df['flatfieldfact'], df['calcoeff']/df['calcoeff'].mean(), c=df['starnumber'], cmap='viridis', alpha=0.3)
    ax[1,0].set_xlabel('1/flatfieldfactor')
    ax[1,0].set_ylabel('Normalised calibration coefficient (uncorrected)')
    ax[1,0].set_title(channel)
    ax[1,0].set_ylim(1-2*stdcalcoeff/meancalcoeff, 1+2*stdcalcoeff/meancalcoeff)
    plt.colorbar(scatter, ax=ax[1,0], label='Starnumber')
    plt.tight_layout()




#remove all instances with outliers
#allstarsdf = allstarsdf[allstarsdf['outlier'] == False]
#Now fit gaussian again seperatly for each star and channel


#Relative calibration
# prepare the dataframe by adding columns for average calibration coefficients for each star and channel, as well as their errors
# Initialize columns for calcoeffmean and calcoeffmean_error in starsdf
starsdf['calcoeffmean'] = None
starsdf['calcoeffstd'] = None
starsdf['calcoeffmean_error'] = None


# Loop through each row in starsdf and calculate calcoeffmean and calcoeffmean_error for each star and channel
for index, row in starsdf.iterrows():
    starnumber = row['starnumber']
    channel = row['channel']
    
    # Filter allstarsdf for matching starnumber and channel
    matching_rows = allstarsdf[(allstarsdf['starnumber'] == starnumber) & (allstarsdf['channel'] == channel)]#& (allstarsdf['outlier'] == False)]


    # Calculate mean and standard deviation for matching rows
    calcoeffmean = matching_rows['calcoeff_corrected'].mean()
    calcoeffstd = matching_rows['calcoeff_corrected'].std()
    calcoeffmean_error = calcoeffstd / np.sqrt(len(matching_rows) - 1)
    
    # Assign values to starsdf
    starsdf.at[index, 'calcoeffmean'] = calcoeffmean
    starsdf.at[index, 'calcoeffstd'] = calcoeffstd
    starsdf.at[index, 'calcoeffmean_error'] = calcoeffmean_error


    if len(matching_rows) > 0:
        
        mu, sigma, ks_stat, p_value = fitgaussian(matching_rows, 'calcoeff_corrected', plot=False, plotxlim=plotxlimdir[channel])
        # Assign values to starsdf
        starsdf.at[index, 'calcoeff_gauss_mu'] =mu
        starsdf.at[index, 'calcoeff_gauss_sigma'] =sigma
        starsdf.at[index, 'calcoeff_gauss_mu_error'] =sigma / np.sqrt(len(matching_rows) - 1)
        starsdf.at[index, 'calcoeff_gauss_ks_stat'] =ks_stat
        starsdf.at[index, 'calcoeff_gauss_p_value'] =p_value

# Display the first few rows of starsdf
print(starsdf.head())


starsdf['relcal'] = np.nan
starsdf['relcal_error'] = np.nan

#stare_stars=[4,5,6,20]
#starsdf=starsdf[starsdf['starnumber'].isin(stare_stars)]
for index in starsdf.index:
    istarnumber=starsdf.loc[index, 'starnumber']
    # for IR channels:
    if starsdf.loc[index, 'channel'] in ['IR1', 'IR2', 'IR3', 'IR4']:
        #check if istarnumber is in the list of starnumbers for channel IR2
        if istarnumber not in starsdf[starsdf['channel'] == 'IR2'].starnumber.values:
            continue

        # pick out the coefficiens for starnumber =istarnumer and channel=IR2
        refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'IR2')].calcoeffmean.values[0]
        #error_in_refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'IR2')].calcoeffmean_error.values[0]
        error_in_refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'IR2')].calcoeffstd.values[0]
        starsdf['relcal'].loc[index] = starsdf.loc[index, 'calcoeffmean'] / refcoeff
        #starsdf['relcal_error'].loc[index] = np.sqrt((starsdf.loc[index, 'calcoeffmean_error']/starsdf.loc[index, 'calcoeffmean'])**2 + (error_in_refcoeff/refcoeff)**2)
        starsdf['relcal_error'].loc[index] = np.sqrt((starsdf.loc[index, 'calcoeffstd']/starsdf.loc[index, 'calcoeffmean'])**2 + (error_in_refcoeff/refcoeff)**2)
    # for UV channels:
    if starsdf.loc[index, 'channel'] in ['UV1', 'UV2']:
        #check if istarnumber is in the list of starnumbers for channel UV1
        if istarnumber not in starsdf[starsdf['channel'] == 'UV2'].starnumber.values:
            continue

        # pick out the coefficiens for starnumber =istarnumer and channel=UV1
        refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'UV2')].calcoeff_gauss_mu.values[0]
        #error_in_refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'UV2')].calcoeffmean_error.values[0]
        error_in_refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'UV2')].calcoeff_gauss_sigma.values[0]
        starsdf['relcal'].loc[index] = starsdf.loc[index, 'calcoeffmean'] / refcoeff
        #starsdf['relcal_error'].loc[index] = np.sqrt((starsdf.loc[index, 'calcoeffmean_error']/starsdf.loc[index, 'calcoeffmean'])**2 + (error_in_refcoeff/refcoeff)**2)
        starsdf['relcal_error'].loc[index] = np.sqrt((starsdf.loc[index, 'calcoeffstd']/starsdf.loc[index, 'calcoeffmean'])**2 + (error_in_refcoeff/refcoeff)**2)


#loop through the channels in figdict and axdict and plot the relative calibration coefficients
for channel in figdict.keys():
    fig = figdict[channel]
    ax=axdict[channel]
    #plot relcal with errorbars
    ax[1,1].errorbar(starsdf[starsdf['channel'] == channel].starnumber, starsdf[starsdf['channel'] == channel].relcal, yerr=starsdf[starsdf['channel'] == channel].relcal_error, fmt='o')
    ax[1,1].set_xlabel('Star number')
    ax[1,1].set_ylabel('Relative calibration coefficient')
    mean = starsdf[starsdf['channel'] == channel].relcal.mean()
    std = starsdf[starsdf['channel'] == channel].relcal.std()
    err_of_mean = std / np.sqrt(len(starsdf[starsdf['channel'] == channel])-1)
    #inverse-variance weighting is used to calculate the mean and the variance of the mean
    var_ivw = 1 / np.sum(1 / starsdf[starsdf['channel'] == channel].relcal_error**2)
    mean_ivw = np.sum(starsdf[starsdf['channel'] == channel].relcal / starsdf[starsdf['channel'] == channel].relcal_error**2) * var_ivw
    err_of_mean_ivw = np.sqrt(var_ivw)
    ax[1,1].text(0.1, 0.9, 'Relative calibration :', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.85, f'Mean: {mean:.2f}, Std: {std:.2f}, Error of mean: {err_of_mean:.2f}', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.8, f'Std in %: {std/mean*100:.2f}, Error of mean in %: {err_of_mean/mean*100:.2f}', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.75, f'Error of mean in %: {err_of_mean/mean*100:.2f}', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.7, f'Mean_ivw: {mean_ivw:.2f}, Error of mean_ivw: {err_of_mean_ivw:.2f}', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.65, f'Error of mean ivw in %: {err_of_mean_ivw/mean_ivw*100:.2f}', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.4, 'Relative calibration from absolute values:', transform=ax[1,1].transAxes)
    relcalmean = calcoeffmeandict[channel]/calcoeffmeandict['IR2']
    relcalmeanerror = np.sqrt((calcoeffmeanerrordict[channel]/calcoeffmeandict[channel])**2 + (calcoeffmeanerrordict['IR2']/calcoeffmeandict['IR2'])**2)
    relcalmeanerrorfromstd = np.sqrt((calcoeffstddict[channel]/calcoeffmeandict[channel])**2 + (calcoeffstddict['IR2']/calcoeffmeandict['IR2'])**2)
    ax[1,1].text(0.1, 0.35, f'Mean: {relcalmean:.3f}, Error of mean: {relcalmeanerror:.3f}, Error using std: {relcalmeanerrorfromstd:.3f}', transform=ax[1,1].transAxes)
    ax[1,1].text(0.1, 0.30, f'Error of mean in %: {relcalmeanerror/relcalmean*100:.3f}, Error using std in %: {relcalmeanerrorfromstd/relcalmean*100:.3f}', transform=ax[1,1].transAxes)            
    ax[1,1].text(0.1, 0.25, f'Abs Mean: {calcoeffmeandict[channel]:.3f}, Abs Error of mean: {calcoeffmeanerrordict[channel]:.3f}, Abs Std: {calcoeffstddict[channel]:.3f}', transform=ax[1,1].transAxes, color='grey')

plt.tight_layout()

#%%
relcaldict = {'IR2': 1, 'UV2': 1}
relcalerrordict = {'IR2': 0, 'UV2': 0}
#for channel in starsdf['channel'].unique():
for channel in ['IR1', 'IR3', 'IR4', 'UV1']:
    # Filter data for IR1 and IR2 channels
    ichannel_data = starsdf[starsdf['channel'] == channel]
    if channel in ['IR1', 'IR3', 'IR4']:
        ir2_data = starsdf[starsdf['channel'] == 'IR2']

        # Merge IR1 and IR2 data on starnumber
        mergeddf = pd.merge(ichannel_data[['starnumber', 'calcoeff_gauss_mu', 'calcoeff_gauss_sigma','calcoeff_gauss_mu_error', 'calcoeff_gauss_ks_stat', 'calcoeff_gauss_p_value']],
                            ir2_data[['starnumber', 'calcoeff_gauss_mu', 'calcoeff_gauss_sigma','calcoeff_gauss_mu_error', 'calcoeff_gauss_ks_stat', 'calcoeff_gauss_p_value']],
                            on='starnumber',
                            suffixes=('_'+channel, '_IR2'))
        mergeddf = mergeddf.dropna(subset=['calcoeff_gauss_mu_IR2', 'calcoeff_gauss_mu_'+channel])
        xcoeff=mergeddf['calcoeff_gauss_mu_IR2']
        xcoefferr=mergeddf['calcoeff_gauss_mu_error_IR2']
    elif channel in ['UV1']:
        uv2_data = starsdf[starsdf['channel'] == 'UV2']

        # Merge IR1 and IR2 data on starnumber
        mergeddf = pd.merge(ichannel_data[['starnumber', 'calcoeff_gauss_mu', 'calcoeff_gauss_sigma','calcoeff_gauss_mu_error', 'calcoeff_gauss_ks_stat', 'calcoeff_gauss_p_value']],
                            uv2_data[['starnumber', 'calcoeff_gauss_mu', 'calcoeff_gauss_sigma','calcoeff_gauss_mu_error', 'calcoeff_gauss_ks_stat', 'calcoeff_gauss_p_value']],
                            on='starnumber',
                            suffixes=('_'+channel, '_UV2'))
        mergeddf = mergeddf.dropna(subset=['calcoeff_gauss_mu_UV2', 'calcoeff_gauss_mu_'+channel])
        xcoeff=mergeddf['calcoeff_gauss_mu_UV2']
        xcoefferr=mergeddf['calcoeff_gauss_mu_error_UV2']

    
    ycoeff=mergeddf['calcoeff_gauss_mu_'+channel]

    ycoefferr=mergeddf['calcoeff_gauss_mu_error_'+channel]
    # Plotting
    fig3, ax3 = plt.subplots()
    ax3.errorbar(xcoeff, ycoeff, xerr=xcoefferr, yerr=ycoefferr, fmt='o', ecolor='r', capthick=2)
    # Fit the line and get the slope and error estimates
    slope, std_err = fit_line_through_origin(xcoeff, ycoeff, xcoefferr, ycoefferr)

    x = np.linspace(xcoeff.min(), xcoeff.max(), 100)
    y = slope * x 
    ax3.plot(x, y, 'k--', label=f'Fit: $y={slope:.2f}x$, $Standard\ error={std_err:.2f}$') 
    std_err_percent = std_err / slope * 100

    ax3.legend()


    if channel in ['IR1', 'IR3', 'IR4']:
        ax3.set_xlabel('Calibration coefficient IR2 [ph cm$^{-2}$ nm$^{-1} /$ counts]')
    elif channel in ['UV1']:
        ax3.set_xlabel('Calibration coefficient UV2 [ph cm$^{-2}$ nm$^{-1} /$ counts]')
    # set label with unit [ph cm$^{-2}$ nm$^{-1} /$ counts]
    ax3.set_ylabel('Calibration coefficient ' + channel+' [ph cm$^{-2}$ nm$^{-1} /$ counts]')
    ax3.grid(True)
    plt.show()
    print(f'The slope of the line is {slope:.2f} with a standard error of {std_err:.2f} ({std_err_percent:.2f}%)')
    # Save the figure in directory ../output
    fig3.savefig(f'../output/{channel}_vs_IR2_relative_calib.png')
    relcaldict[channel] = slope
    relcalerrordict[channel] = std_err

print('*********** FINAL ABSOLUTE CALIBRATION COEFFICIENTS ***********')
print(calcoeffmeandict)
print('*********** FINAL ABSOLUTE CALIBRATION COEFFICIENT ERRORS ***********')
print(calcoeffstddict)
print('*********** FINAL RELATIVE CALIBRATION COEFFICIENTS ***********')
print(relcaldict)
print('*********** FINAL RELATIVE CALIBRATION COEFFICIENT ERRORS ***********')
print(relcalerrordict)


# %%
#Investigate flatfield correction
for channel in allstarsdf['channel'].unique():
    df = allstarsdf[allstarsdf['channel'] == channel]

    fig, ax = plt.subplots()
    scatter = ax.scatter(df['column'], df['calcoeff']/df['calcoeff'].mean(), alpha=0.3, color='b', label='not corrected')
    scatter = ax.scatter(df['column'], df['calcoeff_corrected']/df['calcoeff_corrected'].mean(), alpha=0.3, color='r', label='flatfield corrected')
    ax.set_xlabel('Column where star appears')
    ax.set_ylabel('Normalised calibration coefficient')
    ax.set_title(channel)
    ax.set_ylim(1-6*stdcalcoeff/meancalcoeff, 1+3*stdcalcoeff/meancalcoeff)
    ax.legend()
    plt.tight_layout()
# %%

#%%

# print the relative calibration coefficients with 1 standard deviation error bars in absolute values and percentages
for channel in starsdf['channel'].unique():
    print(channel)
    for index in starsdf[starsdf['channel'] == channel].index:
        if np.isnan(starsdf.loc[index, 'relcal']):
            continue
        print(starsdf.loc[index, 'name'], starsdf.loc[index, 'relcal'])
        mean = starsdf[starsdf['channel'] == channel].relcal.mean()
        std = starsdf[starsdf['channel'] == channel].relcal.std()
        err_of_mean = std / np.sqrt(len(starsdf[starsdf['channel'] == channel])-1)
        #inverse-variance weighting is used to calculate the mean and the variance of the mean
        var_ivw = 1 / np.sum(1 / starsdf[starsdf['channel'] == channel].relcal_error**2)
        mean_ivw = np.sum(starsdf[starsdf['channel'] == channel].relcal / starsdf[starsdf['channel'] == channel].relcal_error**2) * var_ivw
        err_of_mean_ivw = np.sqrt(var_ivw)

        

    print('Mean:', mean, 'Std:', std, 'Error of mean:', err_of_mean, 'Std in %:', std/mean*100, 'Error of mean in %:', err_of_mean/mean*100)
    print('Mean_ivw:', mean_ivw, 'Error of mean_ivw:', err_of_mean_ivw, 'Std in %:', std/mean_ivw*100, 'Error of mean in %:', err_of_mean_ivw/mean_ivw*100)
 


#%%

plt.plot(allstarsdf[(allstarsdf['channel'] == 'IR2') & (allstarsdf['starnumber'] == 1)]['count_corrected'], 'ro')

# %%
allstarsdf[(allstarsdf['channel'] == 'UV1')].starnumber.unique()

# %%
