
#%%
import pandas as pd
import numpy as np




# Define the path to the Excel file
file_path = 'all_stars_20240513lm.xlsx'

# Read the Excel file into a DataFrame
starsdf = pd.read_excel(file_path, sheet_name='Allstars')


#remove all rows with NaNs in the flux column or the count column
starsdf = starsdf.dropna(subset=['flux', 'count'])
starsdf = starsdf.dropna(subset=['name'])
print(starsdf)

#Add column with the exposure time according to the channel
#IR1: 5 s, IR2: 3 s, IR3: 1.5 s, IR4: 1.5 s, UV1: 5 s, UV2: 5 s
#Cal is already calculated in the Excel file, but we reverse engineer it here
starsdf['exptime'] = starsdf['channel'].map({'IR1': 5, 'IR2': 3, 'IR3': 1.5, 'IR4': 1.5, 'UV1': 5, 'UV2': 5})
starsdf['exptime2'] = starsdf['Cal']*starsdf['count']/starsdf['flux'] 
starsdf['calcoeff']=starsdf['flux']/starsdf['count']*starsdf['exptime']
starsdf['error_in_calcoeff']=starsdf['error_in_counts']/starsdf['count']*starsdf['calcoeff']


# %%
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
        refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'IR2')].calcoeff.values
        error_in_calcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'IR2')].error_in_calcoeff.values
        starsdf['relcal'].loc[index] = starsdf.loc[index, 'calcoeff'] / refcoeff
        starsdf['relcal_error'].loc[index] = np.sqrt((starsdf.loc[index, 'error_in_calcoeff']/starsdf.loc[index, 'calcoeff'])**2 + (error_in_calcoeff/refcoeff)**2)
    # for UV channels:
    if starsdf.loc[index, 'channel'] in ['UV1', 'UV2']:
        #check if istarnumber is in the list of starnumbers for channel UV1
        if istarnumber not in starsdf[starsdf['channel'] == 'UV2'].starnumber.values:
            continue

        # pick out the coefficiens for starnumber =istarnumer and channel=UV1
        refcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'UV2')].calcoeff.values
        error_in_calcoeff = starsdf[(starsdf['starnumber'] == istarnumber) & (starsdf['channel'] == 'UV2')].error_in_calcoeff.values
        starsdf['relcal'].loc[index] = starsdf.loc[index, 'calcoeff'] / refcoeff
        starsdf['relcal_error'].loc[index] = np.sqrt((starsdf.loc[index, 'error_in_calcoeff']/starsdf.loc[index, 'calcoeff'])**2 + (error_in_calcoeff/refcoeff)**2)

# %%
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
 



    
# %%
import pandas as pd
import numpy as np

# Create the dataframe from the provided data
data = {
    'starnumber': [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 16, 20, 23],
    'Pulkovo': [2491, 2326, 5340, '5459, 5460', np.nan, 1708, 1713, 2943, 472, 5267, 2990, 4730, 4763],
    'name': ['Sirius', 'Canopus', 'Arcturus', 'Rigil Kent.', 'Vega', 'Capella', 'Rigel', 'Procyon', 'Archenar', 'Hadar', 'Pollux', 'Acrux', 'Gacrux'],
    'class': ['A', 'F', 'K', 'G', 'A', 'G', 'B', 'F', 'B', 'B', 'K', 'B', 'M'],
    'flux': [18791.25, 12067.49, 15818.84, 13183.00, 5231.20, 9601.09, 4472.03, 5411.20, 3018.41, 2479.00, 4057.25, 2095.33, 8198.07],
    'count': [9271.1, 6110.5, 7839.8, 6325.5, 2488.0, 4826.6, 2258.1, 2824.1, 1449.6, 1189.2, 2066.5, 1045.0, 4062.1],
    'error_in_counts': [17.8, 12.0, 33.0, 12.0, 28.0, 50.6, 16.0, 16.0, 118.0, 12.3, 13.8, 16.0, 22.0],
    'Cal': [10.134315, 9.874388, 10.088803, 10.420520, 10.512862, 9.946018, 9.902197, 9.580397, 10.411182, 10.422973, 9.816719, 10.025502, 10.090926],
    'err_percent': [0.191994, 0.196383, 0.420929, 0.189708, 1.125402, 1.048357, 0.708560, 0.566552, 8.140177, 1.034309, 0.667796, 1.531100, 0.541592],
    'slope': [-7.5, -5.7, -9.6, -7.9, 0.9, 1.4, -3.5, -3.5, -6.5, 0.8, -0.9, -2.0, -10.9],
    'slope_err': [2.0, 1.5, 3.7, 1.4, 3.5, 8.6, 2.0, 2.0, 20.0, 1.5, 1.73, 1.90, 2.77],
    'rel_slope': [-0.000809, -0.000933, -0.001225, -0.001249, 0.000362, 0.000290, -0.001550, -0.001239, -0.004484, 0.000673, -0.000436, -0.001914, -0.002683],
    'rel_Err': [0.000216, 0.000245, 0.000472, 0.000221, 0.001407, 0.001782, 0.000886, 0.000708, 0.013797, 0.001261, 0.000837, 0.001818, 0.000682],
    'N': [34, 50, 12, 50, 6, 28, 51, 37, 3, 26, 46, 17, 20],
    'channel': ['IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1', 'IR1'],
    'meastype': [np.nan, np.nan, np.nan, 'stare', 'stare', 'stare', np.nan, np.nan, np.nan, np.nan, np.nan, 'stare', np.nan],
    'exptime': [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
    'exptime2': [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
    'calcoeff': [10.134315, 9.874388, 10.088803, 10.420520, 10.512862, 9.946018, 9.902197, 9.580397, 10.411182, 10.422973, 9.816719, 10.025502, 10.090926],
    'error_in_calcoeff': [0.019457, 0.019392, 0.042467, 0.019769, 0.118312, 0.104270, 0.070163, 0.054278, 0.847489, 0.107806, 0.065556, 0.153501, 0.054652]
}

starsdf = pd.DataFrame(data)

# Calculate the best estimate of the calibration constant
weights = starsdf['error_in_calcoeff'] ** -2
calcoeff_best = (starsdf['calcoeff'] * weights).sum() / weights.sum()

# Calculate the error in the best estimate of the calibration constant
error_in_calcoeff_best = (weights.sum()) ** -0.5

print(f"The best estimate of the calibration constant is {calcoeff_best:.6f}")
print(f"The error in the best estimate of the calibration constant is {error_in_calcoeff_best:.6f}")
# %%
