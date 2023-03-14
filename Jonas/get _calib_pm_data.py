#%%
from mats_utils.rawdata.read_data import read_MATS_PM_data
from calibrate_pm import calibrate_pm
from matplotlib import pyplot as plt
from datetime import datetime, timezone

# Import level1a pm data
start_date = datetime(2023, 2, 24, 00, 00, tzinfo=timezone.utc)
end_date = datetime(2023, 2, 24, 15, 55, tzinfo=timezone.utc)
df = read_MATS_PM_data(start_date,end_date)
df.columns
#%%
df_cal = calibrate_pm(df)
df.columns

#%% Plot calibrated data
# Plot temperatures, calibrated
plt.plot(df['PMTime'], df['pmBkg_Tpd'], label="Photodiode Bkg photometer", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['pmBkg_Tif'], label="Filter Bkg photometer", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['pmAband_Tpd'], label="Photodiode A-band photometer", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['pmAband_Tif'], label="Filter A-band photometer", marker='+', markersize=2.0, linewidth=0)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Temperature [°C]')
plt.show()

# Plot photometer signal, calibrated
plt.plot(df['PMTime'], df['pmBkg_Sig'], label="Bkg Photometer", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['pmAband_Sig'], label="A-band Photometer", marker='+', markersize=2.0, linewidth=0)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Signal [photons cm-2 str-1 s-1]')
plt.show()

# Plot photometer ratio vs satlat and satlon (higher value means higher emission altitude, 0.6 = surface; 13.5 = 10 km cloud top)
plt.figure()
plt.scatter(df['satlon'],df['satlat'],c=df['pmAband_Sig']/df['pmBkg_Sig'], marker='.',clim=[0.6,1.4])
plt.colorbar()
