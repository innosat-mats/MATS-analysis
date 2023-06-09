#%%
import numpy as np
from pywt import wavedecn, waverecn, wavedec
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from matplotlib import pyplot as plt

#%%
starttime = DT.datetime(2023, 3, 1, 14, 0)   
endtime = DT.datetime(2023, 3, 1, 16, 0)
ccd_data_cal = read_MATS_data(starttime, endtime,filter=None,level='1b',version='0.4')

#%%
data = ccd_data_cal.ImageCalibrated[ccd_data_cal.CCDSEL==1]
num_samples=1024
subset = np.stack(data.values,axis=2)[:,:,0:num_samples]
levels = int(np.log2(num_samples))
subbands = wavedec(subset, "haar", mode='constant', level=levels)
# %%

#mean image
plt.imshow(subbands[0],origin='lower',aspect=0.1)
plt.colorbar()
plt.title('Lowest temporal frequency (DC)')

#%%
#gradient image
plt.imshow(subbands[1],origin='lower',aspect=0.1)
plt.colorbar()
plt.title('Lowest temporal frequency (level 0)')

# %%
#high frequency image
plt.title('Highest temporal frequency (level 11)')
plt.imshow(subbands[10][:,:,18],origin='lower',aspect=0.1)
plt.colorbar()

# %%
plt.title('Highest temporal frequency (level 11)')
plt.imshow(subbands[10][:,:,17],origin='lower',aspect=0.1)
plt.colorbar()

# %%
plt.title('STD of Highest temporal frequency (level 11)')
plt.plot(np.std(subbands[10][:,:,:],axis=(0,1)))

# %%
plt.title('Highest temporal frequency (level 11) image 198')
plt.imshow(subbands[10][:,:,198],origin='lower',aspect=0.1)
plt.colorbar()
# %%
plt.title('Highest temporal frequency (level 11) image 198')
plt.imshow(subbands[10][:,:,110],origin='lower',aspect=0.1)
plt.colorbar()
# %%
