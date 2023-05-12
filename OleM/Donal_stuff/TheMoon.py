#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.rawdata.calibration import calibrate_dataframe
import datetime as DT
from mats_l1_processing.instrument import Instrument
#%%
starttime = DT.datetime(2023, 5, 11, 11, 10)   
endtime = DT.datetime(2023, 5, 11, 11, 40)
#starttime = datetime(2023,2, 13, 0, 30, 0)
#endtime = datetime(2023, 2, 13, 12, 45, 0)
ccd_data = read_MATS_data(starttime, endtime,filter=None,level='1a',version='0.5')

#%%
calibration_file ="/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/calibration_data.toml"    
instrument = Instrument(calibration_file)

# #%%
l1b_data = calibrate_dataframe(ccd_data,instrument)
# %%
#ccd_data_cal = read_MATS_data(starttime, endtime,filter=None,level='1b',version='0.4')

#%%

print('tmp')