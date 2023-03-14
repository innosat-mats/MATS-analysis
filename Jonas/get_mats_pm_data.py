#%%
from mats_utils.rawdata.read_data import read_MATS_PM_data
from datetime import datetime, timezone
import matplotlib.pylab as plt

start_date = datetime(2023, 2, 24, 00, 00, tzinfo=timezone.utc)
end_date = datetime(2023, 2, 24, 15, 55, tzinfo=timezone.utc)

# Get the l1a photometer data (uncalibrated data)
df = read_MATS_PM_data(start_date,end_date)

#%%
# Plot temperatures, uncalibrated
plt.plot(df['PMTime'], df['PM1A']/df['PM1ACNTR'], label="PM1_Tphotodiode", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM1B']/df['PM1BCNTR'], label="PM1_Tifilter", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM2A']/df['PM2ACNTR'], label="PM2_Tphotodiode", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM2B']/df['PM2BCNTR'], label="PM2_Tifilter", marker='+', markersize=2.0, linewidth=0)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Temperature [bits]')
plt.show()

# Plot photometer signal, uncalibrated
plt.plot(df['PMTime'], df['PM1S']/df['PM1SCNTR'], label="PM1_Signal, Bkg Phot", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM2S']/df['PM1SCNTR'], label="PM2_Signal, A-band Phot", marker='+', markersize=2.0, linewidth=0)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Signal [bits]')
plt.show()

# %%
# Plot photometer ratio vs satlat and satlon (higher value means higher emission altitude, 0.6 = surface; 13.5 = 10 km cloud top)
plt.figure()
plt.scatter(df['satlon'],df['satlat'],c=df['PM2S']/df['PM1S'], marker='.',clim=[0.6,1.4])
plt.colorbar()
# %%
