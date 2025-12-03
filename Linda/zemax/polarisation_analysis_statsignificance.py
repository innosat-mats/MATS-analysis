

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
# Read files for itest = 1..3 and store arrays in lists
file_h_list = []
file_v_list = []
file_n_list = []
horizontals = []
verticals = []
nonpolariseds = []
horizontals_smoothed = []
verticals_smoothed = []
nonpolariseds_smoothed = []


for itest in range(1, 4):
    fh = os.path.join(base_path, f"5RectanglesXpolarised{channel}_{itest}.TXT")
    fv = os.path.join(base_path, f"5RectanglesYpolarised{channel}_{itest}.TXT")
    fn = os.path.join(base_path, f"5RectanglesNonPolarised{channel}_{itest}.TXT")
    file_h_list.append(fh)
    file_v_list.append(fv)
    file_n_list.append(fn)
    if itest == 1:
        # Keep variables for the first experiment so the rest of the script (unchanged) still runs:
        file_h = file_h_list[0]
        file_v = file_v_list[0]
        file_n = file_n_list[0]

    horizontals.append(np.loadtxt(fh, comments='#', skiprows=26, usecols=range(1, 257), encoding="utf-16"))
    verticals.append(np.loadtxt(fv, comments='#', skiprows=26, usecols=range(1, 257), encoding="utf-16"))
    nonpolariseds.append(np.loadtxt(fn, comments='#', skiprows=26, usecols=range(1, 257), encoding="utf-16"))


#Compute the mean and the standard deviation across the three experiments
horizontal = np.mean(horizontals, axis=0)
vertical = np.mean(verticals, axis=0)
nonpolarised = np.mean(nonpolariseds, axis=0)
horizontal_std = np.std(horizontals, axis=0)
vertical_std = np.std(verticals, axis=0)
nonpolarised_std = np.std(nonpolariseds, axis=0)

#Apply mean filter to smooth the difference image
sizesmooth = 1

horizontal_smoothed = gaussian_filter(horizontal, sigma=sizesmooth)
vertical_smoothed = gaussian_filter(vertical, sigma=sizesmooth)
nonpolarised_smoothed = gaussian_filter(nonpolarised, sigma=sizesmooth)
horizontal_std_smoothed = gaussian_filter(horizontal_std, sigma=sizesmooth)
vertical_std_smoothed = gaussian_filter(vertical_std, sigma=sizesmooth)
nonpolarised_std_smoothed = gaussian_filter(nonpolarised_std, sigma=sizesmooth)

#

# Plot
fig, axes = plt.subplots(6, 1, figsize=(10, 15))

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

# Standard Deviations
im2 = axes[3].imshow(vertical_std, cmap='viridis')
axes[3].set_title("(Vertical Polarisation Standard Deviation)")
plt.colorbar(im2, ax=axes[3])

im3 = axes[4].imshow(horizontal_std, cmap='viridis', clim=im2.get_clim())
axes[4].set_title("(Horizontal Polarisation Standard Deviation)")
plt.colorbar(im3, ax=axes[4])

im4 = axes[5].imshow(nonpolarised_std, cmap='viridis', clim=im2.get_clim())
axes[5].set_title("(Nonpolarised Standard Deviation)")
plt.colorbar(im4, ax=axes[5])

plt.tight_layout()
plt.show()



# Make figure for real signal 
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

# Standard Deviations Smoothed
im2s = axes2[i].imshow(vertical_std_smoothed, cmap='viridis')
axes2[i].set_title("(Vertical Polarisation Standard Deviation Smoothed)")
plt.colorbar(im2s, ax=axes2[i])
i+=1
im3s = axes2[i].imshow(horizontal_std_smoothed, cmap='viridis', clim=im2s.get_clim())
axes2[i].set_title("(Horizontal Polarisation Standard Deviation Smoothed)")
plt.colorbar(im3s, ax=axes2[i])
i+=1
im4s = axes2[i].imshow(nonpolarised_std_smoothed, cmap='viridis', clim=im2s.get_clim())
axes2[i].set_title("(Nonpolarised Standard Deviation Smoothed)")
plt.colorbar(im4s, ax=axes2[i])
i+=1

# Differences between polarisations
im5s = axes2[i].imshow(vertical_smoothed - horizontal_smoothed, cmap='bwr')    
axes2[i].set_title("Difference Vertical - Horizontal Polarisation Smoothed")
plt.colorbar(im5s, ax=axes2[i])
i+=1

# Relative Difference #mask out when denominator is less than a small value (100) to avoid division by zero
mask = (vertical_smoothed + horizontal_smoothed) < 100
relative_difference = np.zeros_like(vertical_smoothed)
relative_difference[~mask] = (vertical_smoothed[~mask] - horizontal_smoothed[~mask]) / (vertical_smoothed[~mask] + horizontal_smoothed[~mask])

im6s = axes2[i].imshow(relative_difference, cmap='bwr', vmin=-0.5, vmax=0.5)
axes2[i].set_title("Relative Difference (V - H) / (V + H) Smoothed")
plt.colorbar(im6s, ax=axes2[i])
plt.tight_layout()
plt.show()


#Now make the same figure for the ghosts, the only diffenerence is that the colorscale limits are different
fig4, axes4 = plt.subplots(8, 1, figsize=(10, 18))
# Vertical Smoothed
i=0
im0s = axes4[i].imshow(vertical_smoothed, cmap='viridis', vmin=0, vmax=250)
axes4[i].set_title("Vertical Polarisation Smoothed (Ghosts)")
plt.colorbar(im0s, ax=axes4[i]) 
i+=1
# Horizontal Smoothed
im1s = axes4[i].imshow(horizontal_smoothed, cmap='viridis', clim=im0s.get_clim())
axes4[i].set_title("Horizontal Polarisation Smoothed (Ghosts)")
plt.colorbar(im1s, ax=axes4[i])
i+=1
#Nonpolarised Smoothed
imns = axes4[i].imshow(nonpolarised_smoothed, cmap='viridis', clim=im0s.get_clim())
axes4[i].set_title("Nonpolarised Smoothed (Ghosts)")
plt.colorbar(imns, ax=axes4[i])
i+=1 
# Standard Deviations Smoothed
im2s = axes4[i].imshow(vertical_std_smoothed, cmap='viridis', vmin=0, vmax=30)
axes4[i].set_title("(Vertical Polarisation Standard Deviation Smoothed) (Ghosts)")
plt.colorbar(im2s, ax=axes4[i])
i+=1
im3s = axes4[i].imshow(horizontal_std_smoothed, cmap='viridis', clim=im2s.get_clim())
axes4[i].set_title("(Horizontal Polarisation Standard Deviation Smoothed) (Ghosts)")
plt.colorbar(im3s, ax=axes4[i])
i+=1
im4s = axes4[i].imshow(nonpolarised_std_smoothed, cmap='viridis', clim=im2s.get_clim())
axes4[i].set_title("(Nonpolarised Standard Deviation Smoothed) (Ghosts)")
plt.colorbar(im4s, ax=axes4[i])
i+=1
# Differences between polarisations
im5s = axes4[i].imshow(vertical_smoothed - horizontal_smoothed, cmap='bwr', vmin=-30, vmax=30)    
axes4[i].set_title("Difference Vertical - Horizontal Polarisation Smoothed (Ghosts)")
plt.colorbar(im5s, ax=axes4[i])
i+=1                
# Relative Difference #mask out when denominator is less than a small value (10) to avoid division by zero
mask = (vertical_smoothed + horizontal_smoothed) < 10
relative_difference = np.zeros_like(vertical_smoothed)
relative_difference[~mask] = (vertical_smoothed[~mask] - horizontal_smoothed[~mask]) / (vertical_smoothed[~mask] + horizontal_smoothed[~mask])          
im6s = axes4[i].imshow(relative_difference, cmap='bwr', vmin=-0.3, vmax=0.3)
axes4[i].set_title("Relative Difference (V - H) / (V + H) Smoothed (Ghosts)")
plt.colorbar(im6s, ax=axes4[i])
plt.tight_layout()
plt.show()


#Make figure that shows the difference between vertical and horizontal polarisation  (smoothedfor the three experiments
fig4, axes4 = plt.subplots(3, 1, figsize=(10, 12))
for itest in range(3):
    #smooth
    horizontalsmooth= gaussian_filter(horizontals[itest], sigma=sizesmooth)
    verticalssmooth= gaussian_filter(verticals[itest], sigma=sizesmooth)
    diff = verticalssmooth - horizontalsmooth
    im = axes4[itest].imshow(diff, cmap='bwr', vmin=-30, vmax=30)
    axes4[itest].set_title(f"Difference Vertical - Horizontal Polarisation Experiment {itest+1}")
    plt.colorbar(im, ax=axes4[itest])
plt.tight_layout()
plt.show()

# %%
