#%%
from mats_utils.geolocation import satellite as sat
import datetime as DT
sat.get_tle_dateDB(DT.datetime.now(),tledb='/home/olemar/Projects/Universitetet/MATS/MATS-utility-functions/data/matsTLE.db')
# %%
