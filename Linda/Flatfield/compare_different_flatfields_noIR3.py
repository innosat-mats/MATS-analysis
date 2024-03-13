
# Script to compare flatfield images taken horizontally and vertically

#%%
from mats_l1_processing.items_units_functions import (
        read_files_in_protocol_as_ItemsUnits,
    )
from database_generation.experimental_utils import readprotocol
from database_generation.flatfield import scale_field, read_flatfield_wo_baffle, read_flatfield_w_baffle

import numpy as np
from PIL import Image
from database_generation.experimental_utils import plot_CCDimage
from mats_l1_processing.read_in_functions import read_CCDitem_rac_or_imgview
from matplotlib import pyplot as plt

#%%
directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/FlatfieldsIR124UV_210908/'
racdir='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/RacFiles_210906-210910/'
protocol='protocol.txt'

read_from = "rac"
df_fullprotocol = readprotocol(directory + protocol)


#df_protocol = df_fullprotocol[df_fullprotocol['channel'] == 'IR1']
df= df_fullprotocol[:40] #take only the images without shutter
channels=['IR1','IR2','IR4','UV1','UV2']
directions=['horizontal', 'vertical']
# Create a dictionaries to store the flatfields
images_vertical = {}
images_horizontal = {}

for direction in directions:
    if direction == 'horizontal':
        df_B = df[:20][df.DarkBright == "B"]
        df_SB = df[:20][df.DarkBright == "SB"]
        df_D = df[:20][df.DarkBright == "D"]
        df_SD = df[:20][df.DarkBright == "DS"] #Order of letters are swapped in the protocol file
    else:
        df_B = df[20:][df.DarkBright == "B"]
        df_SB = df[20:][df.DarkBright == "SB"]
        df_D = df[20:][df.DarkBright == "D"]
        df_SD = df[20:][df.DarkBright == "DS"]


    for i, channel in enumerate(channels):
        imageitm_B = read_CCDitem_rac_or_imgview(racdir, df_B.PicID.iloc[i], read_from)
        imageitm_SB = read_CCDitem_rac_or_imgview(racdir, df_SB.PicID.iloc[i], read_from)
        imageitm_D = read_CCDitem_rac_or_imgview(racdir, df_D.PicID.iloc[i], read_from)
        imageitm_SD = read_CCDitem_rac_or_imgview(racdir, df_SD.PicID.iloc[i], read_from)

        image= imageitm_B['IMAGE']-imageitm_SB['IMAGE']-(imageitm_D['IMAGE']-imageitm_SD['IMAGE'])

        # plot_CCDimage(imageitm_B['IMAGE'], title='Bright '+channel)
        # plot_CCDimage(imageitm_SB['IMAGE'], title='Short Bright '+channel)
        # plot_CCDimage(imageitm_D['IMAGE'], title='Dark '+channel)
        # plot_CCDimage(imageitm_SD['IMAGE'], title='Short Dark '+channel)
        # plot_CCDimage(image, title='Resulting image '+channel)
        
        if direction == 'horizontal':
            images_horizontal[channel] = image
        else:
            images_vertical[channel] = image

#%%
# Read in flatfields without baffle from warmenr temperatures
directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200506_flatfields_roomtemp/'
protocol='protocol.txt'
df_fullprotocol = readprotocol(directory + protocol)

#create a dict with shorter protocols for each channel
df_protocols = {}

channel='IR1' 
df_protocols[channel] = df_fullprotocol[df_fullprotocol.PicID.str.endswith('1')]
channel='IR2'
df_protocols[channel] = df_fullprotocol[df_fullprotocol.PicID.str.endswith('4')]
channel='IR3' 
df_protocols[channel] = df_fullprotocol[df_fullprotocol.PicID.str.endswith('3')]
channel='IR4'
df_protocols[channel] = df_fullprotocol[df_fullprotocol.PicID.str.endswith('2')]
channel='UV1' 
df_protocols[channel] = df_fullprotocol[df_fullprotocol.PicID.str.endswith('5')]
channel='UV2'
df_protocols[channel] = df_fullprotocol[df_fullprotocol.PicID.str.endswith('6')]

flats_room={}
for channel in channels:
    imageitm_B = read_CCDitem_rac_or_imgview(directory, df_protocols[channel].PicID.iloc[0], read_from)
    imageitm_D1 = read_CCDitem_rac_or_imgview(directory, df_protocols[channel].PicID.iloc[1], read_from)
    imageitm_D2 = read_CCDitem_rac_or_imgview(directory, df_protocols[channel].PicID.iloc[2], read_from)
    image= imageitm_B['IMAGE']-0.5*(imageitm_D1['IMAGE']+imageitm_D2['IMAGE'])
    flats_room[channel] = image




#%%

#Set calibration file to read in the flatfields that are currently used in the retrieval
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
flatfields_wo_baffle = {}
flatfields_w_baffle = {}
for channel in channels:
    flatfields_wo_baffle[channel], flatfield_wo_baffle_err=read_flatfield_wo_baffle(calibration_file, channel, sigmode='HSM', reporterror=True)
    flatfields_w_baffle[channel]=read_flatfield_w_baffle(calibration_file, channel)

channels=['IR1','IR2','IR4','UV1','UV2']

for channel in channels:
    fig, ax=plt.subplots(18,1, figsize=(10,40))

    #flip the image if it is IR1, IR3,UV1 or UV2
    if channel in ['IR1','IR3','UV1','UV2']:
        image_vertical= np.fliplr(images_vertical[channel])
        image_horizontal= np.fliplr(images_horizontal[channel])
        flatfield_wo_baffle = np.fliplr(flatfields_wo_baffle[channel])
        flatfield_w_baffle = np.fliplr(flatfields_w_baffle[channel])
        flat_room = np.fliplr(flats_room[channel])
    else:
        image_vertical= images_vertical[channel]
        image_horizontal= images_horizontal[channel]
        flatfield_wo_baffle = flatfields_wo_baffle[channel]
        flatfield_w_baffle = flatfields_w_baffle[channel]
        flat_room = flats_room[channel]


    vertical_image_scaled = scale_field(image_vertical)
    horizontal_image_scaled = scale_field(image_horizontal)
    flatfield_wo_baffle_scaled = scale_field(flatfield_wo_baffle)
    flatfield_w_baffle_scaled = scale_field(flatfield_w_baffle)
    flat_room_scaled= scale_field(flat_room)






    plot_CCDimage(horizontal_image_scaled,fig=fig, axis=ax[2], title='Horizontal')
    ax[3].plot(np.mean(horizontal_image_scaled[200:300,300:1750],0), label='Horizontal mean')
    plot_CCDimage(vertical_image_scaled,fig=fig, axis=ax[0], title='Vertical')
    ax[1].plot(np.mean(vertical_image_scaled[200:300,300:1750],0), label='Vertical mean')
    diffimage = horizontal_image_scaled-vertical_image_scaled
    plot_CCDimage(diffimage,fig=fig, axis=ax[4], title='Horizontal-Vertical')
    ax[5].plot(np.mean(diffimage[200:300,300:1750],0), label='Horizontal-Vertical mean')

    plot_CCDimage(flatfield_wo_baffle_scaled,fig=fig, axis=ax[6], title='Flatfield without baffle')
    ax[7].plot(np.mean(flatfield_wo_baffle_scaled[200:300,300:1750],0), label='Flatfield without baffle mean')
    diffimage = flatfield_wo_baffle_scaled-vertical_image_scaled
    plot_CCDimage(diffimage,fig=fig, axis=ax[8], title='Flatfield without baffle - Vertical')
    ax[9].plot(np.mean(diffimage[200:300,300:1750],0), label='Flatfield without baffle - Vertical mean')
    plot_CCDimage(flatfield_w_baffle_scaled,fig=fig, axis=ax[10], title='Flatfield with baffle')
    ax[11].plot(np.mean(flatfield_w_baffle_scaled[200:300,300:1750],0), label='Flatfield with baffle mean')
    diffimage = flatfield_w_baffle_scaled-vertical_image_scaled
    plot_CCDimage(diffimage,fig=fig, axis=ax[12], title='Flatfield with baffle - Vertical')
    ax[13].plot(np.mean(diffimage[200:300,300:1750],0), label='Flatfield with baffle - Vertical mean')

    plot_CCDimage(flat_room_scaled,fig=fig, axis=ax[14], title='Flatfield roomtemp')
    ax[15].plot(np.mean(flat_room_scaled[200:300,300:1750],0), label='Flatfield roomtemp mean')
    diffimage = flat_room_scaled-flatfield_wo_baffle_scaled
    plot_CCDimage(diffimage,fig=fig, axis=ax[16], title='Flatfield roomtemp - Flatfield without baffle')
    ax[17].plot(np.mean(diffimage[200:300,300:1750],0), label='Flatfield roomtemp - Flatfield without baffle mean')

    plt.suptitle(channel)

    plt.show()


    #Save the plot
    fig.savefig('../output/compare_different_flatfields_'+channel+'_excl_ir3.png')

# %%
