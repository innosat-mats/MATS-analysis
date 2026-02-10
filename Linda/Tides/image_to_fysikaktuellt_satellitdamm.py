
#%%import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import requests
from io import StringIO

# Set publication-quality parameters
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
rcParams['font.size'] = 11
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 14
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10

# Fetch data from Our World in Data
print("Fetching satellite launch data...")
url = "https://ourworldindata.org/grapher/yearly-number-of-objects-launched-into-outer-space.csv"
response = requests.get(url)
data = pd.read_csv(StringIO(response.text))

print("Data columns:", data.columns.tolist())
print("First few rows:")
print(data.head(20))

# Save raw data
data.to_csv('/home/sandbox/satellite_launches_raw_data.csv', index=False)
print("\nRaw data saved to satellite_launches_raw_data.csv")


# %%
