#%%
import glob
import numpy as np

folder = "/Users/lindamegner/MATS/MATS-retrieval/data/atmos_smallerfiles/"
files = sorted(glob.glob(folder + "*.dat"))

total_rows = 0
total_intensity = 0.0

for f in files:
    data = np.genfromtxt(f, comments="!", skip_header=1)
    if data.ndim == 0 or data.size == 0:
        print(f"{f.split('/')[-1]}: 0 rows, intensity = 0")
        continue
    if data.ndim == 1:
        data = data[np.newaxis, :]
    n_rows = len(data)
    intensity_sum = data[:, -1].sum()
    total_rows += n_rows
    total_intensity += intensity_sum
    intensity_per_row = intensity_sum / n_rows
    print(f"{f.split('/')[-1]}: {n_rows} rows, intensity = {intensity_sum:.4e}, intensity/row = {intensity_per_row:.4e}")

print(f"\nTotal rows:      {total_rows}")
print(f"Total intensity: {total_intensity:.4e}")

# %%
