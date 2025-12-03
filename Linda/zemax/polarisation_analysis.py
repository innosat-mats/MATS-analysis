

# %%
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter
from scipy.ndimage import gaussian_filter


base_path = "/Users/lindamegner/MATS/MATS-retrieval/data/zemaxoutput/"

# # Construct full paths safely using os.path.join
# file_h = os.path.join(base_path, "IR1RectangleHorizontallyPolarised.TXT")
# file_v = os.path.join(base_path, "IR1RectangleVerticallyPolarised.TXT")
# file_n = os.path.join(base_path, "IR1RectangleNonPolarised.TXT")
channel = "IR1"
file_h = os.path.join(base_path, f"5RectanglesXpolarised{channel}.TXT")
file_v = os.path.join(base_path, f"5RectanglesYpolarised{channel}.TXT")
file_n = os.path.join(base_path, f"5RectanglesNonPolarised{channel}.TXT")

# Load the data
horizontal = np.loadtxt(file_h, comments='#', skiprows=26, usecols=range(1, 257), encoding="utf-16")
vertical = np.loadtxt(file_v, comments='#', skiprows=26, usecols=range(1, 257), encoding="utf-16")
nonpolarised = np.loadtxt(file_n, comments='#', skiprows=26, usecols=range(1, 257), encoding="utf-16")

#Apply mean filter to smooth the difference image
sizesmooth = 5
horizontal_smoothed = gaussian_filter(horizontal, sigma=sizesmooth)
vertical_smoothed = gaussian_filter(vertical, sigma=sizesmooth)
nonpolarised_smoothed = gaussian_filter(nonpolarised, sigma=sizesmooth)




# Plot
fig, axes = plt.subplots(4, 1, figsize=(10, 15))

# Vertical
im0 = axes[0].imshow(vertical, cmap='viridis')
axes[0].set_title("Vertical Polarisation")
plt.colorbar(im0, ax=axes[0])

# Horizontal
im1 = axes[1].imshow(horizontal, cmap='viridis', clim=im0.get_clim())
axes[1].set_title("Horizontal Polarisation")
plt.colorbar(im1, ax=axes[1])

#Nonpolarised
imn = axes[2].imshow(nonpolarised, cmap='viridis', clim=im0.get_clim())
axes[2].set_title("Nonpolarised")
plt.colorbar(imn, ax=axes[2])

# Difference
im2 = axes[3].imshow(difference, cmap='bwr', vmin=-1, vmax=1)
axes[3].set_title("(Horizontal - Vertical) / Nonpolarised")
plt.colorbar(im2, ax=axes[3])

plt.tight_layout()
plt.show()

fig2, axes2 = plt.subplots(8, 1, figsize=(10, 18))
# Vertical Smoothed
i=0
im0s = axes2[i].imshow(vertical_smoothed, cmap='viridis')
axes2[i].set_title("Vertical Polarisation Smoothed")
plt.colorbar(im0s, ax=axes2[i]) 
i+=1
# Horizontal Smoothed
im1s = axes2[i].imshow(horizontal_smoothed, cmap='viridis', clim=im0s.get_clim())
axes2[i].set_title("Horizontal Polarisation Smoothed")
plt.colorbar(im1s, ax=axes2[i])
i+=1
#Nonpolarised Smoothed
imns = axes2[i].imshow(nonpolarised_smoothed, cmap='viridis', clim=im0s.get_clim())
axes2[i].set_title("Nonpolarised Smoothed")
plt.colorbar(imns, ax=axes2[i])
i+=1

# Mean of Horizontal and Vertical minus Nonpolarised
mean_hv = (horizontal_smoothed + vertical_smoothed) / 2
mean_minus_nonpolarised = mean_hv - nonpolarised_smoothed
im2m = axes2[i].imshow(mean_minus_nonpolarised, cmap='bwr')
axes2[i].set_title("Mean(Horizontal, Vertical) - Nonpolarised Smoothed")
plt.colorbar(im2m, ax=axes2[i])
i+=1
im2m = axes2[i].imshow(mean_minus_nonpolarised, cmap='bwr', clim=(-50, 50))
axes2[i].set_title("Mean(Horizontal, Vertical) - Nonpolarised Smoothed smaller scale")
plt.colorbar(im2m, ax=axes2[i])
i+=1

# Difference 
difference = horizontal_smoothed - vertical_smoothed
im2a = axes2[i].imshow(difference, cmap='bwr')
axes2[i].set_title("Difference (Smoothed Horizontal - Smoothed Vertical)")
plt.colorbar(im2a, ax=axes2[i])
i+=1
# Difference smaller scale
difference_smaller_scale = horizontal_smoothed - vertical_smoothed
im2a = axes2[i].imshow(difference_smaller_scale, cmap='bwr', vmin=-25, vmax=25)
axes2[i].set_title("Difference (Smoothed Horizontal - Smoothed Vertical) smaller scale")
plt.colorbar(im2a, ax=axes2[i])
i+=1
# Difference Smoothed
im2s = axes2[i].imshow(smoothed_difference, cmap='bwr', vmin=-1, vmax=1)
axes2[i].set_title("Smoothed (Horizontal - Vertical) / Nonpolarised")
plt.colorbar(im2s, ax=axes2[i])


plt.tight_layout()
plt.show()

# %%
