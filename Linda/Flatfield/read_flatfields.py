#%%
#read in numpy arrays

import numpy as np
import matplotlib.pyplot as plt
directory='/Users/lindamegner/MATS/MATS-retrieval/MATS-instrument-data/calibration_data/flatfields/'

channels=['IR1','IR2','IR3','IR4','UV1','UV2']#,'NADIR' ]
for channel in channels:
       
    file='flatfield_'+channel+'_HSM.npy'
    flatfield=np.load(directory+file)
    plt.imshow(flatfield)
    plt.colorbar()
    plt.title('Flatfield '+channel)
    plt.show()

    fig,ax=plt.subplots(1,2)
    #plot middle column
    ax[0].plot(flatfield[:,1000], '.')
    ax[0].set_title('Middle column of flatfield '+channel)

    #plot middle row
    ax[1].plot(flatfield[250,-40:], '.')
    ax[1].set_title('Middle row of flatfield '+channel)
    plt.show()


#%%
# %%
