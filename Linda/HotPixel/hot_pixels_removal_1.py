


#%%

from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
import datetime as DT
import pickle

#%%


starttime=DT.datetime(2023,3,22,9,15,0)
endtime=DT.datetime(2023,3,22,9,25,0)


#filter={'CCDSEL': [1, 1], 'NRBIN': [1, 1], 'NCBINCCDColumns': [1, 1], 'NCOL':[2047,2048], 'NROW':[511,512]} 
#dfl1b = read_MATS_data(starttime, endtime,filter,level='1b',version='0.4')
dfl1a = read_MATS_data(starttime, endtime,level='1a',version='1.0',dev=False)

#save the dataframe as a pickle file
with open('../output/hotpix_removal/dfl1a_hot_pixel_removal.pkl', 'wb') as f:
    pickle.dump(dfl1a, f)





# %%
