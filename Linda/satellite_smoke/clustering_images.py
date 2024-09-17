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

def cluster_images(image_data, method='kmeans'):
    if method == 'kmeans':
        return kmeans_cluster(image_data)
    elif method == 'dbscan':
        return dbscan_cluster(image_data)
    else:
        raise ValueError("Invalid clustering method")

def kmeans_cluster(image_data):       
    # Apply K-Means clustering
    kmeans = KMeans(n_clusters=5, random_state=42)
    kmeans.fit(image_data)
    # Get the cluster labels
    labels = kmeans.labels_
    return labels

def dbscan_cluster(image_data):
    determine_eps = True
    if determine_eps:
        # # Fit NearestNeighbors to determine the optimal eps value
        n_neighbors = 2
        neighbors = NearestNeighbors(n_neighbors=n_neighbors)
        neighbors_fit = neighbors.fit(image_data)
        distances, indices = neighbors_fit.kneighbors(image_data)

        # Sort the distances
        distances = np.sort(distances, axis=0)
        k_distances = distances[:, n_neighbors - 1]

        # Plot the k-distance graph to help choose eps
        plt.plot(k_distances)
        plt.xlabel('Points sorted by distance')
        plt.ylabel(f'Distance to {n_neighbors}-th nearest neighbor')
        plt.title('k-distance Graph')
        plt.show()
    # Apply DBSCAN with adjusted parameters
    eps_value =  1700 # Adjust this value based on the k-distance graph
    min_samples_value = 4  # Adjust this value based on your dataset

    # Apply DBSCAN clustering
    dbscan = DBSCAN(eps=eps_value, min_samples=min_samples_value)
    labels = dbscan.fit_predict(image_data)
    return labels


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
field = 'ImageCalibrated'
dfo = dfo[~(((dfo['satlat'] <= 0) & (dfo['satlat'] >= -60) & ((dfo['satlon'] > 300) | (dfo['satlon'] < 30))))]

df_day = dfo[(dfo['TPsza'] <= 90)]
df_night = dfo[(dfo['TPsza'] >= 100)]

df = df_day.copy()
df.reset_index(drop=True, inplace=True)
#%%

# Flatten the image data
image_data = np.array([img.flatten() for img in df[field]])

# Perform PCA to reduce dimensionality
#pca = PCA(n_components=50)
#image_data_reduced = pca.fit_transform(image_data)

method='dbscan'
labels=cluster_images(image_data, method=method)

print('labels', labels)

# %%


# Add labels to the DataFrame
df["Cluster"] = labels

# Display results

unique_labels = set(labels)
nimages = 4
fig, ax = plt.subplots( nimages,len(unique_labels)-1, figsize=(10, 10))
for label in unique_labels:
    if label == -1:
        # Noise
        continue
    cluster_images = df[df["Cluster"] == label][field]
    for i in range(min(nimages, len(cluster_images))):
        plot_CCDimage(cluster_images.iloc[i] , fig=fig, axis=ax[i, label])
        ax[i, label].set_title(f"Cluster {label}")
plt.show()

if method=='dbscan':
    # Identify outliers
    nimg = 10
    outlier_images = df[df["Cluster"] == -1][field]
    for i in range(min(nimg, len(outlier_images))):
        plot_CCDimage(outlier_images.iloc[i], title="Outlier")
        



# %%

# #old code for dbscanning


