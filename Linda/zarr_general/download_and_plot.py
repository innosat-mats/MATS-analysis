#%%

#python get_zarr.py -c IR1 -b 2023 2 20 0 0 0 -e 2023 3 1 0 0 0  -f TPlon 10 30


import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

outputpath="/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output"
# Assign web address
path = 'https://bolin.su.se/data/s3/data/mats-level-1b-limb-cropd-1.0/mats-level-1b-limb-cropd-IR1.zarr'

# Get xarray Dataset using zarr engine
data = xr.open_zarr(path)

# Print data structure
print(data)

# Retrieve first time step
time = data.time[0].values

# Retrieve first calibrated image
image = data.ImageCalibrated[0,:,:].values
# %%
meantpheight=data["TPheight"].mean(dim="time")
print(meantpheight)


#%% Select the Time Slice
subdata = data.sel(time=slice('2023-02-16T15:00:00', '2023-02-16T20:00:00'))
#%%
fig, ax = plt.subplots()
img=ax.imshow(subdata['ImageCalibrated'].isel(time=0), cmap='viridis', aspect=1/4)
ax.invert_yaxis()  # Invert y-axis to match image orientation
ax.set_title(str(subdata.time.values[0]))

def animate(i):
    img.set_data(subdata['ImageCalibrated'].isel(time=i))
    ax.set_title(str(subdata.time.values[i]))
    return img, ax

ani = animation.FuncAnimation(fig, animate, frames=len(subdata.time), interval=20, blit=True, repeat=False)

#%%
print('hello')
# Save the animation as a GIF file
ani.save(f'{outputpath}/mats_animation.gif', writer='pillow', fps=10)






# %%

def play_movie(images, interval=0.1):
    import matplotlib.pyplot as plt
    import time

    fig, ax = plt.subplots()
    im = ax.imshow(images.isel(time=0), cmap='gray')
    title = ax.set_title(str(images.time.values[0]))

    plt.ion()
    plt.show()

    for i in range(len(images.time)):
        im.set_array(images.isel(time=i).values)
        title.set_text(str(images.time.values[i]))
        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(interval)

    plt.ioff()
    plt.show()

#%%
images = subdata['ImageCalibrated']
play_movie(images, interval=0.1) # Adjust interval for speed


# %%
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

# Assuming your dataset is called "ds"
# and satlat is (time) and satheight is also (time)
# Check if satheight exists, otherwise adjust the variable name
lat = data["satlat"]
height = data["satheight"]

# Make sure time is a datetime index
time_index = pd.DatetimeIndex(data.time.values)

# Group into 15-day bins
groups = time_index.to_period("15D")

# Iterate through groups
for period, idx in pd.Series(range(len(time_index)), index=groups).groupby(groups):
    subset_lat = lat.isel(time=idx.values)
    subset_height = height.isel(time=idx.values)
    
    # Skip if data is empty or all NaN
    if subset_lat.notnull().sum() == 0 or subset_height.notnull().sum() == 0:
        continue

    plt.figure(figsize=(8, 6))
    plt.scatter(subset_lat, subset_height, s=1, alpha=0.5)
    plt.xlabel("Satellite Latitude")
    plt.ylabel("Satellite Height")
    plt.title(f"satlat vs satheight ({period.start_time.date()} → {period.end_time.date()})")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# %%
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Variables
lat = data["satlat"]
height = data["satheight"]

# Ensure time is datetime
time_index = pd.to_datetime(data.time.values)

# Get start and end of dataset
start_time = time_index.min()
end_time = time_index.max()

# Create 15-day edges starting from the first timestamp
bin_edges = pd.date_range(start=start_time.normalize(), end=end_time, freq="1D")
if bin_edges[-1] < end_time:
    bin_edges = bin_edges.append(pd.DatetimeIndex([end_time]))

# Number of bins
n_groups = len(bin_edges) - 1

# Subplot layout
ncols = 3
nrows = int(np.ceil(n_groups / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows), sharex=True, sharey=True)
axes = axes.flatten()

for i in range(n_groups):
    t0, t1 = bin_edges[i], bin_edges[i+1]
    mask = (time_index >= t0) & (time_index < t1)
    
    subset_lat = lat.sel(time=mask)
    subset_height = height.sel(time=mask)
    
    ax = axes[i]
    if subset_lat.notnull().sum() > 0 and subset_height.notnull().sum() > 0:
        ax.scatter(subset_lat, subset_height, s=1, alpha=0.5)
    ax.set_title(f"{t0.date()} → {t1.date()}")
    ax.set_xlabel("satlat")
    ax.set_ylabel("satheight")
    ax.grid(True)

# Hide any unused subplots
for ax in axes[n_groups:]:
    ax.axis("off")

fig.suptitle("Satellite Latitude vs Height (non-overlapping 15-day chunks)", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()



# %%
# select the weird data, namely from the 24th of april 2023 to the 25th of april, and plot satheight as a function of time

weird_data = data.sel(time=slice("2023-04-24", "2023-04-25"))
plt.figure(figsize=(10, 6))
plt.plot(weird_data.time, weird_data["satheight"], marker='o', linestyle='-')
plt.title("Satellite Height from 24th to 25th April 2023")
plt.xlabel("Time")
plt.ylabel("Satellite Height")
plt.grid()
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# %%
