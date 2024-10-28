
# File to compare flatfields of different temperatures

#%%

from database_generation import flatfield as flatfield
from database_generation.experimental_utils import plot_CCDimage
import toml
from mats_l1_processing.instrument import CCD
import numpy as np
from matplotlib import pyplot as plt

#from mats_l1_processing.instrument import Instrument
def read_flatfield_wo_baffle(calibration_file, channel, sigmode='HSM', reporterror=False,specifyflatfiledfolder=None):

    CCDunit=CCD(channel,calibration_file)
    calibration_data=toml.load(calibration_file)
    if specifyflatfiledfolder is not None:
        flatfield_wo_baffle, flatfield_wo_baffle_err = flatfield.read_flatfield(CCDunit, sigmode, calibration_data["flatfield"][specifyflatfiledfolder], reporterror=True)
    else:
        flatfield_wo_baffle, flatfield_wo_baffle_err = flatfield.read_flatfield(CCDunit, sigmode, calibration_data["flatfield"]["flatfieldfolder_cold_unprocessed"], reporterror=True)
    if reporterror:
        return flatfield_wo_baffle, flatfield_wo_baffle_err
    else:
        return flatfield_wo_baffle


def readmyflatfield(directory, protocol, channel, expt):
    from database_generation.experimental_utils import readprotocol
    from mats_l1_processing.items_units_functions import (
        read_files_in_protocol_as_ItemsUnits,
    )
    read_from = "rac"
    df_protocol = readprotocol(directory + protocol)
    # df_only2 = df_protocol[(df_protocol.index-2) % 3 != 0]

    # The below reads all images in protocol - very inefficient. Should be only one file read in LM200810
    CCDItemsUnits = read_files_in_protocol_as_ItemsUnits(
        df_protocol, directory, 3, read_from
    )
    if channel == "NADIR":  # Hack since we dont have any nadir flat fields yet.
        raise Warning('No flatfields measurements of the NADIR channel')
        flatfield = np.ones((511, 2048))
    else:
        CCDItemsUnitsChannel = list(
            filter(lambda x: (x.imageItem["channel"] == channel), CCDItemsUnits)
        )
        CCDItemsUnitsSelect = list(
            filter(lambda x: (x.imageItem["TEXPMS"] == expt), CCDItemsUnitsChannel)
        )


    flatfield = CCDItemsUnitsSelect[0].dark
    print('Warning: dark is used instead of flatfield')
    expt=CCDItemsUnitsSelect[0].imageItem["TEXPMS"]
    print('expt:', expt)
    return flatfield

def calcmeanforcentrarea(flatfield):
    #Calculate mean for central area
    #central area
    FirstRow = 350 
    LastRow = 400 
    FirstCol =524
    LastCol =1523
    flatfieldmean=np.mean(flatfield[FirstRow:LastRow,FirstCol:LastCol])
    return flatfieldmean

#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'




	
flatfieldfolder_cold_unprocessed = "/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/calibration_data/20200330_flatfields_0C/"
flatfieldfolder_8C_unprocessed = "/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/calibration_data/20200429_flatfields_8C/"
flatfieldfolder_roomt_unprocessed = "/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/calibration_data/20200506_flatfields_roomtemp/"


#%%
#instrument=Instrument(calibration_file)
temperaturearray=np.array([0,20])
expt=12000

channel='IR1'
flatfield0C=readmyflatfield(flatfieldfolder_cold_unprocessed, 'flatfields_200330_SigMod1_LMprotocol.txt', channel, expt)
plot_CCDimage(flatfield0C, title='fullframe 0C '+channel+' '+str(expt/1000)+'s')
#%%
flatfield8C=readmyflatfield(flatfieldfolder_8C_unprocessed, 'protocolcleaned.txt', channel, expt)
plot_CCDimage(flatfield8C, title='fullframe 8C '+channel+' '+str(expt)+'s')
#%%
flatfieldroomt=readmyflatfield(flatfieldfolder_roomt_unprocessed, 'protocol.txt', channel,expt)
plot_CCDimage(flatfieldroomt, title='fullframe roomt '+channel+' '+str(expt/1000)+'s')


#%%

bias=293 #Bias is 293 for every readout
lsb_to_electron=33
print('flatfield0C mean in electrons per s:', (flatfield0C.mean()-bias)/(expt/1000)*lsb_to_electron)
print('flatfield8C mean in electrons per s:', (flatfield8C.mean()-bias)/(expt/1000)*lsb_to_electron)
print('flatfieldroomt mean in electrons per s:', (flatfieldroomt.mean()-bias)/(expt/1000)*lsb_to_electron)


# dark current increase with temperature
print('dark current increase from 0 to 20 C:', (flatfieldroomt.mean()-flatfield0C.mean())/(expt/1000)*lsb_to_electron)

