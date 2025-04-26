#%%
import pickle
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from database_generation.experimental_utils import plot_CCDimage
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans

# Load the data
idiff = 1
channel = 'IR1'
nrimages = 1600
seed = 42
crop = 'CROPD'
with open(f'testdata/df_random_allchannels_{crop}_idiff_{idiff}_nimg_{nrimages}_seed_{seed}.pkl', 'rb') as f:
    dfchannelsdict = pickle.load(f)

dfo = dfchannelsdict[channel]

# Select which field to search for anomalies in and cut the dataset
field = 'ImageCalibratedDiff' + str(idiff)
dfo = dfo[~(((dfo['satlat'] <= 0) & (dfo['satlat'] >= -60) & ((dfo['satlon'] > 300) | (dfo['satlon'] < 30))))]

df_day = dfo[(dfo['TPsza'] <= 90)]
df_night = dfo[(dfo['TPsza'] >= 100)]

df = df_day.copy()
df.reset_index(drop=True, inplace=True)


#%%
# Assuming df['IMAGE'] contains the image data as numpy arrays
field='ImageCalibrated'
image_data = np.array(df[field].tolist())

# Flatten the images if necessary
image_data = image_data.reshape(len(image_data), -1)

# Apply K-Means clustering
kmeans = KMeans(n_clusters=5, random_state=42)
kmeans.fit(image_data)

# Get the cluster labels
labels = kmeans.labels_
# %%
