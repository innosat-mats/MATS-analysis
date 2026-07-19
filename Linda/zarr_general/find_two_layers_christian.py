#%%

#python get_zarr.py -c IR1 -b 2023 2 20 0 0 0 -e 2023 3 1 0 0 0  -f TPlon 10 30


import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

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

#%% Select the Time Slice
subdata = data.sel(time=slice('2023-02-18T20:00:00', '2023-02-18T23:00:00'))
#%%


def _smooth_profile(values, window=9):
    values = np.asarray(values, dtype=float)
    if values.size < 3:
        return values.copy()
    if window is None:
        window = max(5, (values.size // 50) | 1)
    window = max(3, int(window))
    if window % 2 == 0:
        window += 1
    window = min(window, values.size if values.size % 2 == 1 else values.size - 1)
    if window < 3:
        return values.copy()

    kernel = np.ones(window, dtype=float) / window
    padded = np.pad(values, (window // 2,), mode="reflect")
    return np.convolve(padded, kernel, mode="valid")

def _merge_close_peaks(peak_indices, peak_values, min_distance):
    if len(peak_indices) <= 1:
        return list(peak_indices)
    if len(peak_values) <= 1:
        return list(peak_values)


    order = np.argsort(peak_indices)
    peak_indices = np.asarray(peak_indices)[order]
    peak_values = np.asarray(peak_values)[order]

    merged_idx = [int(peak_indices[0])]
    merged_val = [float(peak_values[0])]

    for idx, val in zip(peak_indices[1:], peak_values[1:]):
        if idx - merged_idx[-1] <= min_distance:
            if val > merged_val[-1]:
                merged_idx[-1] = int(idx)
                merged_val[-1] = float(val)
        else:
            merged_idx.append(int(idx))
            merged_val.append(float(val))

    return merged_idx

def check_layer(column, rows=None, row_threshold=100, return_fit=False,
    smooth_window=9, prominence_factor=2.5):
    column = np.asarray(column, dtype=float)
    if column.ndim != 1:
        raise ValueError("check_layer expects a 1-D column")
    n = column.size
    if n < 3:
        if return_fit:
            return 0, lambda x_new: np.full_like(np.asarray(x_new, dtype=float), np.nan, dtype=float)
        return 0

    smooth = _smooth_profile(column, window=smooth_window)
    if smooth.size != n:
        smooth = column.copy()

    residual = column - smooth
    mad = np.median(np.abs(residual - np.median(residual)))
    noise = 1.4826 * mad
    if not np.isfinite(noise) or noise <= 0:
        noise = np.std(residual)
    if not np.isfinite(noise) or noise <= 0:
        noise = max(1.0, 0.01 * np.ptp(smooth))

    # Use an index-based cutoff as primary behavior for robustness.
    i_max_idx = min(n - 1, max(2, int(row_threshold) - 1))

    # If rows are provided and useful, allow them to expand the search window.
    if rows is not None:
        rows = np.asarray(rows, dtype=float)
        if rows.shape == column.shape:
            candidates = np.where(rows < row_threshold)[0]
            if candidates.size > 0:
                i_max_idx = max(i_max_idx, int(np.max(candidates)))

    y = smooth[: i_max_idx + 1]
    if y.size < 3:
        if return_fit:
            x = np.arange(n, dtype=float)
            return 0, lambda x_new: np.interp(np.asarray(x_new, dtype=float), x, smooth, left=smooth[0], right=smooth[-1])
        return 0

    peaks = []

    # Edge peak at top
    if y[0] > y[1]:
        peaks.append(0)

    dy = np.diff(y)
    for i in range(1, y.size - 1):
        if dy[i - 1] > 0 and dy[i] <= 0:
            peaks.append(i)

    # Edge peak at bottom side of the selected interval
    if y[-1] > y[-2]:
        peaks.append(y.size - 1)

    prom_window = max(6, y.size // 8)
    prom_floor = max(prominence_factor * noise, 0.06 * np.ptp(y))

    kept_idx = []
    kept_val = []

    for idx in sorted(set(peaks)):
        left = max(0, idx - prom_window)
        right = min(y.size, idx + prom_window + 1)

        if idx == 0:
            base = np.min(y[1:right]) if right - 1 > 0 else y[idx]
        elif idx == y.size - 1:
            base = np.min(y[left:idx]) if idx - left > 0 else y[idx]
        else:
            left_min = np.min(y[left:idx]) if idx > left else y[idx]
            right_min = np.min(y[idx + 1:right]) if idx + 1 < right else y[idx]
            base = max(left_min, right_min)

        prominence = y[idx] - base
        if prominence >= prom_floor:
            kept_idx.append(int(idx))
            kept_val.append(float(y[idx]))

    if len(kept_idx) > 1:
        min_dist = max(5, y.size // 10)
        kept_idx = _merge_close_peaks(kept_idx, kept_val, min_dist)

    layer_count = min(2, len(kept_idx))

    x = np.arange(n, dtype=float)
    def fit_function(x_new):
        x_new = np.asarray(x_new, dtype=float)
        return np.interp(x_new, x, smooth, left=smooth[0], right=smooth[-1])

    if return_fit:
        return layer_count, fit_function
    return layer_count

def check_layers(*args, **kwargs):
    return check_layer(*args, **kwargs)

#%%

midcolumn = subdata["ImageCalibrated"].isel(time=slice(850, 900)).isel(im_col=subdata.sizes["im_col"] // 2)

#plot vs time and im_row
fig, ax = plt.subplots()
img=ax.imshow(midcolumn.T, cmap='viridis', aspect=1/8)
ax.invert_yaxis()  # Invert y-axis to match image orientation
ax.set_title('Middle Column vs Time and Row')
ax.set_xlabel('Time')
ax.set_ylabel('Row')
plt.colorbar(img, ax=ax, label='Middle Column Value')
plt.show()

#%%
#plot every 10th profile
fig, ax = plt.subplots()

for i in range(0, midcolumn.sizes["time"],2):
    ax.plot(midcolumn.isel(time=i), midcolumn["im_row"], label=f'Time {i}')
ax.set_title('Middle Column Profiles at Different Times')
ax.set_xlabel('Middle Column Value')
ax.set_ylabel('Row')
ax.legend()
plt.show()
#%%

#plot certain timesteps 
fig, ax = plt.subplots()
for i in range(37,40):
    ax.plot(midcolumn.isel(time=i), midcolumn["im_row"], label=f'Time {i}')
ax.set_title('Middle Column Profile')
ax.set_xlabel('Middle Column Value')
ax.set_ylabel('Row')
ax.legend()
plt.show()


# %%
#check if there are two layers in the middle column
rows = midcolumn["im_row"].values
layer_categories = {0: [], 1: [], 2: []}

for time_index in range(30, 40):
    profile = midcolumn.isel(time=time_index).values
    layer_count = check_layers(profile, row_threshold=120, rows=rows, smooth_window=5, prominence_factor=2.0)
    if layer_count in layer_categories and len(layer_categories[layer_count]) < 1:
        layer_categories[layer_count].append(time_index)

#%%

fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

for layer_count, ax in zip([0, 1, 2], axes):
    time_indices = layer_categories[layer_count]
    if not time_indices:
        ax.set_title(f"No example found for {layer_count} maxima")
        ax.set_ylabel("Signal")
        continue

    time_index = time_indices[0]
    profile = midcolumn.isel(time=time_index).values
    fit_count, fit_function = check_layers(profile, rows=rows, row_threshold=115, return_fit=True)
    fitted_profile = fit_function(np.arange(profile.size, dtype=float))

    ax.plot(rows, profile, color="tab:blue", alpha=0.6, label="Raw profile")
    ax.plot(rows, fitted_profile, color="tab:orange", linewidth=2, label="Smoothed fit")
    ax.set_title(f"Example with {fit_count} maxima at time index {time_index}")
    ax.set_ylabel("Signal")
    ax.legend()

axes[-1].set_xlabel("Row")
plt.tight_layout()
plt.show()
# %%
