#script to try out xarray

#%%
import xarray as xr
import numpy as np

# %%
temperature_data = np.random.rand(4, 3)
temperature = xr.DataArray(
    temperature_data,
    dims=["time", "location"],
    coords={"time": ["2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04"], "location": ["A", "B", "C"]},
    name="temperature"
)
print(temperature)


# %%
precipitation_data = np.random.rand(4, 3)
precipitation = xr.DataArray(
    precipitation_data,
    dims=["time", "location"],
    coords={"time": ["2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04"], "location": ["A", "B", "C"]},
    name="precipitation"
)
print(precipitation)
# %%
dataset = xr.Dataset(
    {
        "temperature": temperature,
        "precipitation": precipitation
    }
)
print(dataset)

# %%
dataset["temperature"].mean(dim="time")
# %%
