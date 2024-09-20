#%%
from keras.applications.vgg16 import VGG16, preprocess_input
from keras.preprocessing import image
from keras.models import Model
import numpy as np

# Load pre-trained VGG16 model + higher level layers
base_model = VGG16(weights='imagenet')
model = Model(inputs=base_model.input, outputs=base_model.get_layer('fc1').output)

def extract_features(img_path):
    img = image.load_img(img_path, target_size=(224, 224))
    img_data = image.img_to_array(img)
    img_data = np.expand_dims(img_data, axis=0)
    img_data = preprocess_input(img_data)
    features = model.predict(img_data)
    return features.flatten()

#%%

with open(f'testdata/df_random_allchannels_{crop}_idiff_{idiff}_nimg_{nrimages}_seed_{seed}.pkl', 'rb') as f:
    dfchannelsdict = pickle.load(f)

df = dfchannelsdict[channel]
#field = 'ImageCalibratedDiff' + str(idiff)
field= 'ImageCalibrated'
df['features'] = df[field].apply(extract_features)

from sklearn.cluster import KMeans

# Convert features to a numpy array
features = np.array(df['features'].tolist())

# Perform K-Means clustering
kmeans = KMeans(n_clusters=5, random_state=0).fit(features)

# Add cluster labels to DataFrame
df['cluster'] = kmeans.labels_

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Reduce dimensions to 2D for visualization
pca = PCA(n_components=2)
reduced_features = pca.fit_transform(features)

# Plot the clusters
plt.scatter(reduced_features[:, 0], reduced_features[:, 1], c=df['cluster'])
plt.show()
