
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


    flatfield = CCDItemsUnitsSelect[0].subpic
    expt=CCDItemsUnitsSelect[0].imageItem["TEXPMS"]
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
temperaturearray=np.array([0,20])
#create dictionary rates for different channels
rates={}
exptdict={'IR1':[3, 4, 5, 6, 7, 12, 17], 'IR2':[3, 4, 5, 6], 'IR3':[3, 4, 5, 6, 7, 12, 22], 'IR4':[3, 4, 5, 6, 7, 12], 'UV1':[4, 5, 6, 7, 12], 'UV2':[3, 4, 5, 6, 7, 12]}

for channel in ['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']:
    fig, ax=plt.subplots(1,1)
    ratesarray=np.array([])


    for expt in exptdict[channel]:
        exptms=1000*expt
        print('channel:', channel, ' expt:', expt)
        flatfield0C=readmyflatfield(flatfieldfolder_cold_unprocessed, 'flatfields_200330_SigMod1_LMprotocol.txt', channel, exptms)
        flatfield0Cmean=calcmeanforcentrarea(flatfield0C)
        flatfieldroomt=readmyflatfield(flatfieldfolder_roomt_unprocessed, 'protocol.txt', channel,exptms)
        flatfieldroomtmean=calcmeanforcentrarea(flatfieldroomt)
        meanmean=(flatfield0Cmean+flatfieldroomtmean)/2.
        flatfieldmean=np.array([flatfield0Cmean/meanmean, flatfieldroomtmean/meanmean])

        rate=(flatfieldroomtmean/meanmean-flatfield0Cmean/meanmean)*100/(20-0)
        ratesarray=np.append(ratesarray, rate)
        #prtint rate rounded to 2 decimals
        ax.plot(temperaturearray, flatfieldmean, label='Int time '+str(expt)+'s, '+ 'percent per degree '+str(round(rate,2)))
        #keep this plot
        
    ax.set_xlabel('Temperature [C]')
    ax.set_ylabel('Normalised Flatfield intensity')
    ax.legend()
    ax.set_title(channel)
    rates[channel]=ratesarray

    #flatfield8C=readmyflatfield(flatfieldfolder_8C_unprocessed, 'flatfields_200429_SigMod1_LMprotocol.txt', channel, expt)




# %%
import pickle
print(rates)
pickle.dump(rates, open('rates.pkl', 'wb'))

# %%
exptdict={'IR1':[3, 4, 5, 6, 7, 12, 17], 'IR2':[3, 4, 5, 6], 'IR3':[3, 4, 5, 6, 7, 12, 22], 'IR4':[3, 4, 5, 6, 7, 12], 'UV1':[3, 4, 5, 6, 7, 12], 'UV2':[3, 4, 5, 6, 7, 12]}
fig, ax=plt.subplots(1,1)
ax.plot(exptdict['IR1'], rates['IR1'], label='IR1')
ax.plot(exptdict['IR2'], rates['IR2'], label='IR2')
ax.plot(exptdict['IR3'], rates['IR3'], label='IR3')
ax.plot(exptdict['IR4'], rates['IR4'], label='IR4')
ax.plot(exptdict['UV1'][1:], rates['UV1'][1:], label='UV1')
ax.plot(exptdict['UV2'], rates['UV2'], label='UV2')


ax.set_xlabel('Integration time [s]')
ax.set_ylabel('Rate [%/C]')
ax.legend()

# %%
