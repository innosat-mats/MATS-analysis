
#%%
import numpy as np
from matplotlib import pyplot as plt

#%%
tabledir = "/home/ochri/Projects/MATS/MATS-L1-processing/calibration_data/linearity/tables/20220906/"
filename_1 = 'highresIR_4.npy'
filename_2 = 'highresIR_3.npy'
filename_3 = 'channel_2_nrbin_6_ncbinfpga_1_ncbinccd_215.npy'
filename_4 = 'channel_3_nrbin_6_ncbinfpga_1_ncbinccd_215.npy'

old_table_3 = np.load(tabledir + filename_2)
new_table_3 = np.load(tabledir + filename_4)
old_table_4 = np.load(tabledir + filename_1)
new_table_4 = np.load(tabledir + filename_3)


# %%
plt.plot(old_table_3[2,:],old_table_3[0,:])
plt.plot(new_table_3[2,:],new_table_3[0,:],'--')
plt.plot(old_table_4[2,:],old_table_4[0,:])
plt.plot(new_table_4[2,:],new_table_4[0,:],'--')
plt.xlabel('Measured value')
plt.ylabel('Corrected value')
plt.legend({"IR3 NRBIN=200","IR3 NRBIN=215","IR4 NRBIN=200","IR4 NRBIN=215"})

# %%
