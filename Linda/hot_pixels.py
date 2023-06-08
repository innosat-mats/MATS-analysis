#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.experimental_utils import plot_CCDimage
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibration_functions import bin_image_with_BC
from database_generation.flatfield import read_pic_and_picd
from mats_utils.rawdata.calibration import  calibrate_dataframe
import datetime as DT
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lindas_own_functions import rename_CCDitem_entries

#%%
def printcorr(img1, img2):
    corrconst=np.corrcoef(img1.ravel(),img2.ravel())
    print(corrconst)
    return corrconst

def plotcorr(img1, img2, texpms):
    fig, ax=plt.subplots(2, 1)
    corrconst=np.corrcoef(img1.ravel(),img2.ravel())
    sp=plot_CCDimage(img1, fig, ax[0], title=channel, clim=[600,610])
    sp=plot_CCDimage(img2, fig, ax[1], title=str(corrconst[0,1])+' expt: '+str(texpms/1000), clim=[600,610])
    plt.tight_layout()
    return corrconst
#%%

# #November dark images
# starttime=DT.datetime(2022,11,30,13,30,0)
# endtime=DT.datetime(2022,11,30,17,30,0)

 # #binned images with hot pixels 
#starttime=DT.datetime(2023,3,30,5,20,0)
#endtime=DT.datetime(2023,3,30,5,25,0)

# #all march
starttime=DT.datetime(2023,1,1,0,0,0)
endtime=DT.datetime(2023,4,1,0,0,0)

filter={'CCDSEL': [1, 1], 'NRBIN': [1, 1], 'NCBINCCDColumns': [1, 1], 'NCOL':[2047,2048], 'NROW':[511,512]} 
#dfl1b = read_MATS_data(starttime, endtime,filter,level='1b',version='0.4')
dfl1a = read_MATS_data(starttime, endtime,filter,level='1a',version='0.5')


#%%
#df=pd.read_pickle('./output/dfl1a')
df=copy.deepcopy(dfl1a)
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
#calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/scripts/calibration_data_linda.toml'
instrument=Instrument(calibration_file)

df=calibrate_dataframe(df, instrument)



rename_CCDitem_entries(df)

print(df.columns.tolist())


#%%
CCDunit=instrument.get_CCD(df.iloc[0]["channel"])


img_tempdep=10**CCDunit.log_a_img_avr_HSM
img_const=10**CCDunit.log_b_img_avr_HSM

startrow=0
imgb_tempdep=bin_image_with_BC(df.iloc[0], img_tempdep)
imgb_const=bin_image_with_BC(df.iloc[0], img_const)
fig, ax=plt.subplots(3,1)
imgname='image_lsb'
testimage=df.iloc[1][imgname][startrow:,:]
#plot_CCDimage(img_tempdep, fig, ax[0], title='temp_dep')
plot_CCDimage(imgb_tempdep[startrow:,:], fig, ax[0], title='temp_dep binned')#, clim=[1.27,1.29])
#plot_CCDimage(img_const, fig, ax[2], title='constant term')
plot_CCDimage(testimage, fig, ax[1], title=imgname, clim=[280,330])
plot_CCDimage(imgb_const[startrow:,:], fig, ax[2], title='constant term binned',clim=[2,7])

printcorr(imgb_tempdep[startrow:,:],testimage)
printcorr(imgb_const[startrow:,:],testimage)
# corrtempdep=np.corrcoef(imgb_tempdep[startrow:,:].ravel(),testimage.ravel())
# print(corrtempdep)
# corrconst=np.corrcoef(imgb_const[startrow:,:].ravel(),testimage.ravel())
# print(corrconst)
#%%


namelist=["image_bias_sub","image_desmeared",
    "image_dark_sub","image_calib_nonflipped","ImageCalibrated"]

fig, ax=plt.subplots(5,1)
for i, name in enumerate(namelist):
    df.iloc[0][name].shape
    plot_CCDimage(df.iloc[0][name], fig, ax[i], title=name)

# fig, ax=plt.subplots(2,1)
# diffimage=df.iloc[0]["image_calib_nonflipped"]-np.fliplr(df.iloc[0]["ImageCalibrated"])
# plot_CCDimage(diffimage, fig, ax[0], title='diffimage')


#%%
iref=13
imagename='image_lsb'
for i in range(len(df)):
    fig, ax=plt.subplots(1)
    refimage=df.iloc[iref][imagename]
    image=df.iloc[i][imagename]
    diffimage=image-refimage
    #image=df.iloc[i][imagename]-df.iloc[iref][imagename]/df.iloc[iref].TEXPMS*df.iloc[i].TEXPMS
    plot_CCDimage(diffimage, fig, ax, title='texpms '+str(df.iloc[i].TEXPMS))
    ax.text(500, 200, 'mean value img= '+str(image[200:300,200:300].mean()))
    ax.text(500, 150, 'mean value ref= '+str(refimage[200:300,200:300].mean()))
    ax.text(500, 100, 'mean value diffimage= '+str(diffimage[200:300,200:300].mean()))




# %%

#read in dark and bright picture from lab (flatfields from laboratory)
from PIL import Image
#directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/Flatfield20210421/PayloadImages/'
channel='IR1'
directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/Flatfield20210421/PayloadImages/'
dfile0='1303052462743530240_1'
picd=np.float64(Image.open(directory + dfile0 + ".pnm"))
pic, picd=read_pic_and_picd(channel,directory)    
fig, ax=plt.subplots(2,1)
sp=plot_CCDimage(pic, fig, ax[0], title=channel)
sp_dark=plot_CCDimage(picd, fig, ax[1], title=channel)
clim_dark=sp_dark.get_clim()



#%%
startrow=0
endrow=512
startcol=0
endcol=2047
for index, row in dfl1a.iterrows():
    texpms=row['TEXPMS']
    plotcorr(picd[startrow:endrow,startcol:endcol],row['IMAGE'][startrow:endrow,startcol:endcol], texpms)
    print(row['temperature_HTR'])

# %%
