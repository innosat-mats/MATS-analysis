#calculate jacobian for all measurements
import numpy as np
import pandas as pd
from mats_l2_processing.forward_model import calc_jacobian
import pickle
from datetime import datetime
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l2_processing.inverse_model import do_inversion

starttime=datetime(2023,3,31,21,0)
stoptime=datetime(2023,3,31,22,35)
dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.4")

df = dftop[dftop['channel'] == 'IR2'].dropna().reset_index(drop=True)

offsets = np.arange(0,len(df)-100,50)
#select part of orbit
for offset in offsets:
    num_profiles = 100 #use 50 profiles for inversion
    df_batch = df.loc[offset:offset+num_profiles-1]
    df_batch = df.reset_index(drop=True)
    columns = np.arange(0,df["NCOL"][0],2)
    rows = np.arange(0,df["NROW"][0]-10,1)

    y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local = calc_jacobian(df_batch,columns,rows)

    y = y.reshape(-1)
    y = np.matrix(y)

    x_hat = do_inversion(ks,y)

    filename = str(offset) + ".pkl"
    with open(filename, "wb") as file:
        pickle.dump((x_hat, y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local), file)


# %% reshape x_hat

#def center_grid(grid):
#    return (grid[:-1]+grid[1:])/2

#altitude_grid = center_grid(altitude_grid_edges)
#alongtrack_grid = center_grid(alongtrack_grid_edges)
#acrosstrack_grid = center_grid(acrosstrack_grid_edges)
#x_hat_reshape1 = np.array(x_hat).reshape(len(altitude_grid),len(acrosstrack_grid),len(alongtrack_grid))




## correct for elliptical earth
# ret_ecef=ecef_to_local.inv().apply(np.array([altitude_grid[0]*np.cos(lats),np.zeros(len(lats)),altitude_grid[0]*np.sin(lats)]).T) #Should this be of [0]? or mid?
# ret_lats=np.rad2deg(cart2sph(ret_ecef)[:,2])


## Plot crossections

# x_hat_reshape1 = np.array(x_hat).reshape(len(altitude_grid)
# ,len(acrosstrack_grid)
# ,len(alongtrack_grid))

# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=0))
# plt.show()
# # %%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=1))
# plt.show()

# #%%
# plt.pcolor(np.sum(x_hat_reshape1[:,:,:],axis=2))
# plt.show()
# # %%