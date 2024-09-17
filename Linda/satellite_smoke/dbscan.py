
"""
Test script for running DBSCAN (Density-Based Spatial Clustering of Applications with Noise).

This algorithm groups together points that are closely packed together, marking points that are in low-density regions as outliers.
It’s particularly useful for data with noise and varying densities.

Author: [Your Name]
Date: [Current Date]
"""
#%%
# Import necessary libraries

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.datasets import make_blobs
from matplotlib import cm
from sklearn.neighbors import NearestNeighbors
# Define your code here

# Generate sample data
X, _ = make_blobs(n_samples=300, centers=4, cluster_std=0.60, random_state=0)
#%%
# Fit NearestNeighbors
n_neighbours=2
neighbors = NearestNeighbors(n_neighbors=n_neighbours)
neighbors_fit = neighbors.fit(X)
distances, indices = neighbors_fit.kneighbors(X)

# Sort the distances (4th column)
distances = np.sort(distances, axis=0)
k_distances = distances[:, n_neighbours-1]

# Plot the k-distance graph
plt.plot(k_distances)
plt.xlabel('Points sorted by distance')
plt.ylabel('Distance to {}-th nearest neighbor'.format(n_neighbours))
plt.title('k-distance Graph')
plt.show()

#%%
# Apply DBSCAN
dbscan = DBSCAN(eps=0.3, min_samples=10)
clusters = dbscan.fit_predict(X)

#%%
# Plot the clusters
plt.scatter(X[:, 0], X[:, 1], c=clusters, cmap='plasma')
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.title('DBSCAN Clustering')


# Add a color bar
norm = plt.Normalize(vmin=min(clusters), vmax=max(clusters))
sm = cm.ScalarMappable(cmap='plasma', norm=norm)
sm.set_array([])
plt.colorbar(sm, label='Cluster Label')
print('labels', clusters)
plt.show()

# %%
