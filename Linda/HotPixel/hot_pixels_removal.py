#%%
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
#%matplotlib widget



# ---------------------------------------------------------------------------
# Hot pixel / single-event removal
# ---------------------------------------------------------------------------
# Strategy: spatial sigma-clipping against a local median background.
#
# 1. Compute local median of each pixel from its (median_size x median_size)
#    neighbourhood using a fast C-backed median filter (no Python loops).
# 2. Estimate the typical noise level from the median absolute deviation (MAD)
#    of (image - local_median).
# 3. Flag pixels where (image - local_median) > sigma * noise as hot/cosmic.
# 4. Replace flagged pixels with the local median value.
#
# This handles both persistent hot pixels and single-event cosmic rays.
# References:
#   astropy CCD guide  https://www.astropy.org/ccd-reduction-and-photometry-guide/
#   scipy.ndimage.median_filter  https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.median_filter.html
# ---------------------------------------------------------------------------

def detect_hot_pixels(image, sigma=5.0, median_size=3):
    """Return boolean mask of hot pixels detected by local sigma clipping.

    Parameters
    ----------
    image : 2-D ndarray
    sigma : float
        Detection threshold in units of the local noise estimate.
    median_size : int (odd)
        Side length of the median-filter neighbourhood.

    Returns
    -------
    mask : bool ndarray, True where a hot pixel / cosmic-ray event was detected
    local_median : float ndarray, the median-filtered image (background estimate)
    """
    local_median = median_filter(image, size=median_size)
    residual = image - local_median
    # MAD-based noise estimate (robust, not influenced by the hot pixels)
    noise = 1.4826 * np.median(np.abs(residual - np.median(residual)))
    if noise == 0:
        noise = 1.0  # avoid division by zero for perfectly flat images
    mask = residual > sigma * noise
    return mask, local_median


def replace_hot_pixels(image, mask, local_median=None, median_size=3):
    """Replace hot pixels with the local median.

    Parameters
    ----------
    image : 2-D ndarray
    mask : bool ndarray  (True = hot pixel)
    local_median : 2-D ndarray, optional
        Pre-computed median image; computed here if not supplied.
    median_size : int (odd)

    Returns
    -------
    corrected : 2-D ndarray with the same dtype as image
    """
    if local_median is None:
        local_median = median_filter(image, size=median_size)
    corrected = image.copy().astype(float)
    corrected[mask] = local_median[mask]
    return corrected


def remove_hot_pixels(image, sigma=5.0, median_size=3):
    """Detect and replace hot pixels in a single image.

    Parameters
    ----------
    image : 2-D ndarray
    sigma : float  detection threshold
    median_size : int  neighbourhood size for median filter

    Returns
    -------
    corrected : 2-D float ndarray
    mask : bool ndarray  (True = pixel was hot)
    """
    mask, local_median = detect_hot_pixels(image, sigma=sigma, median_size=median_size)
    corrected = replace_hot_pixels(image, mask, local_median=local_median)
    return corrected, mask


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
#%%
pkl_path = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/hotpix_removal/dfl1a_hot_pixel_removal.pkl'
dfl1a = pd.read_pickle(pkl_path)
print(f"Loaded {len(dfl1a)} rows, IMAGE shape: {dfl1a.iloc[0]['IMAGE'].shape}")


#select only IR1 channel
dfl1a_IR1 = dfl1a[dfl1a.channel == 'IR1'].reset_index(drop=True)
dfl1a_IR2 = dfl1a[dfl1a.channel == 'IR2'].reset_index(drop=True)
dfl1a_IR3 = dfl1a[dfl1a.channel == 'IR3'].reset_index(drop=True)
dfl1a_IR4 = dfl1a[dfl1a.channel == 'IR4'].reset_index(drop=True)
dfl1a_UV1 = dfl1a[dfl1a.channel == 'UV1'].reset_index(drop=True)
dfl1a_UV2 = dfl1a[dfl1a.channel == 'UV2'].reset_index(drop=True)
print(f"After filtering for IR1: {len(dfl1a)} rows remain")


#save the dataframe as a pickle file
with open('../output/hotpix_removal/dfl1a_IR1_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a_IR1, f)
with open('../output/hotpix_removal/dfl1a_IR2_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a_IR2, f)
with open('../output/hotpix_removal/dfl1a_IR3_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a_IR3, f)
with open('../output/hotpix_removal/dfl1a_IR4_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a_IR4, f)
with open('../output/hotpix_removal/dfl1a_UV1_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a_UV1, f)
with open('../output/hotpix_removal/dfl1a_UV2_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a_UV2, f)   



df = dfl1a_IR2.copy()  # work on a copy of the IR1 subset for now
# ---------------------------------------------------------------------------
# Apply hot pixel removal to every IMAGE in dfl1a
# ---------------------------------------------------------------------------
#%%
SIGMA      = 4.0   # detection threshold (# of MAD-noise units above local median)
MED_SIZE   = 3     # neighbourhood size for the median filter (3×3 pixels)

df = df.copy()


for idx, row in df.iloc[::5].iterrows():
    img = row['IMAGE'].astype(float)
    #plot a histogram of the pixel values to check for outliers
    plt.figure(figsize=(6, 4))
    plt.hist(img[:].flatten(), bins=100, edgecolor='k')
    plt.xlabel('Pixel value')
    plt.ylabel('Number of pixels')
    plt.title(f'Row {idx} – pixel value distribution before hot pixel removal')
    plt.tight_layout()
    plt.show()
    img_corr, mask = remove_hot_pixels(img, sigma=SIGMA, median_size=MED_SIZE)
    #corrected_images.append(img_corr)
    #hot_pixel_masks.append(mask)
    n_hot = mask.sum()
    if n_hot > 0:
        print(f"Row {idx}: {n_hot} hot pixels detected")

#%%
df = df.copy()
df['IMAGE_hotpix_corrected'] = corrected_images
df['hot_pixel_mask']         = hot_pixel_masks


# ---------------------------------------------------------------------------
# Visualisation: compare a sample image before and after
# ---------------------------------------------------------------------------
#%%
sample_idx = 0   # change to inspect a different row
img_raw  = df.iloc[sample_idx]['IMAGE'].astype(float)
img_corr = df.iloc[sample_idx]['IMAGE_hotpix_corrected']
mask     = df.iloc[sample_idx]['hot_pixel_mask']
n_hot    = mask.sum()

fig, axes = plt.subplots(3, 1, figsize=(10, 8))
im0 = axes[0].imshow(img_raw,  aspect='auto', origin='upper')
axes[0].set_title(f'Raw IMAGE (row {sample_idx})')
plt.colorbar(im0, ax=axes[0])

im1 = axes[1].imshow(img_corr, aspect='auto', origin='upper',
                     vmin=img_raw.min(), vmax=img_raw.max())
axes[1].set_title(f'Hot-pixel corrected (sigma={SIGMA}, {n_hot} pixels fixed)')
plt.colorbar(im1, ax=axes[1])

im2 = axes[2].imshow(mask.astype(float), aspect='auto', origin='upper', cmap='Reds')
axes[2].set_title('Hot pixel mask (red = flagged)')
plt.colorbar(im2, ax=axes[2])

plt.tight_layout()
plt.show()

print(f"Sample image: {n_hot} hot pixels detected and replaced "
      f"(fraction: {n_hot / img_raw.size:.2e})")


# ---------------------------------------------------------------------------
# Plot every 100th image: raw vs corrected side by side
# ---------------------------------------------------------------------------
#%%
step = 20
indices = list(range(0, len(df), step))
n = len(indices)

fig, axes = plt.subplots(n, 2, figsize=(10, 3 * n))
if n == 1:
    axes = axes[np.newaxis, :]  # ensure 2-D even for a single row

for row_ax, i in zip(axes, indices):
    raw  = df.iloc[i]['IMAGE'].astype(float)
    corr = df.iloc[i]['IMAGE_hotpix_corrected']
    mask = df.iloc[i]['hot_pixel_mask']
    n_hot = mask.sum()
    vmin, vmax = raw.min(), raw.max()

    im0 = row_ax[0].imshow(raw,  aspect='auto', origin='upper', vmin=vmin, vmax=2000)
    row_ax[0].set_title(f'Row {i} – raw')
    plt.colorbar(im0, ax=row_ax[0])

    im1 = row_ax[1].imshow(corr, aspect='auto', origin='upper', vmin=vmin, vmax=2000)
    row_ax[1].set_title(f'Row {i} – corrected ({n_hot} hot px)')
    plt.colorbar(im1, ax=row_ax[1])


plt.suptitle(f'Every {step}th image: raw (left) vs hot-pixel corrected (right)', y=1.01)
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# Optional: inspect how hot the detected pixels were
# ---------------------------------------------------------------------------
#%%
if n_hot > 0:
    raw_vals  = img_raw[mask]
    corr_vals = img_corr[mask]
    excess    = raw_vals - corr_vals
    print(f"\nHot pixel excess statistics (raw - corrected):")
    print(f"  min  = {excess.min():.1f}")
    print(f"  mean = {excess.mean():.1f}")
    print(f"  max  = {excess.max():.1f}")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(excess, bins=50, edgecolor='k')
    ax.set_xlabel('Excess counts (raw - local median)')
    ax.set_ylabel('Number of pixels')
    ax.set_title(f'Hot pixel excess distribution – row {sample_idx}')
    plt.tight_layout()
    plt.show()

# %%
