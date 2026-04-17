#%%
import xarray as xr

filepath = "/Users/lindamegner/MATS/MATS-retrieval/data/tomotestdata/2023-04-25-09-29-27_L2_tomo.nc"

ds = xr.open_dataset(filepath)
print(ds)

# %% Print all variables
for var in ds.data_vars:
    print(var, ds[var])

#%%
import matplotlib.pyplot as plt
import numpy as np

# Plot 5 altitude levels from 80 to 100 km, every 5 km
target_alts_m = np.arange(80e3, 101e3, 5e3)  # 80, 85, 90, 95, 100 km in metres
alt_indices = [int(np.argmin(np.abs(ds.alt_coord.values - a))) for a in target_alts_m]

fig, axes = plt.subplots(5, 1, figsize=(12, 20))


for ax, alt_idx in zip(axes, alt_indices):
    data = ds['T'].isel(alt_coord=alt_idx).squeeze()
    pcm = ax.pcolormesh(
        ds.alongtrack_coord / 1e3,
        ds.acrosstrack_coord / 1e3,
        data,
        shading='auto',
        cmap='viridis'
    )
    plt.colorbar(pcm, ax=ax, label='T')
    ax.set_xlabel('Along-track [km]')
    ax.set_ylabel('Across-track [km]')
    ax.set_title(f'T | alt = {ds.alt_coord.values[alt_idx]/1e3:.1f} km')

plt.tight_layout()
plt.show()

7h# %%
