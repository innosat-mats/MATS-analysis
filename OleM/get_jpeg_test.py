#%%
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
from mats_utils.imagetools.imagetools import bin_image

import pandas as pd
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from matplotlib import pyplot as plt
import struct

import numpy as np
import cv2
from PIL import Image
import tifffile as tiff
import imageio
#%% 
calibration_file ="/home/olemar/Projects/Universitetet/MATS/MATS-L1-processing/scripts/calibration_data.toml"    
instrument = Instrument(calibration_file)

#%% Select on explicit time
start_time = DT.datetime(2022,11 , 23, 10, 0)
stop_time = DT.datetime(2022, 11, 23, 14, 0)

#%%
df = read_MATS_data(start_time,stop_time,version='0.3',level='0/CCD',dev=False)

# %%
n = 12
image = df.loc[n].IMAGE
plt.imshow(image)
plt.title("CCDSEL: " + str(df.loc[n].CCDSEL))
# %%
# simulate binned image
binned_image = bin_image(image,2,40)-(2*40*(df.loc[n].TBLNK-128))+128
# %%
plt.imshow(binned_image)
plt.colorbar()
plt.title("CCDSEL: " + str(df.loc[n].CCDSEL))

# %%
# dig error
error = (binned_image/16).astype(int)*16-binned_image
std_error = np.std(error)
plt.imshow((binned_image/16).astype(int)*16-binned_image)
plt.colorbar()
plt.title(f"digitization error std: {std_error:.2f}")

# %% make into 12 bit image
array = binned_image.astype(np.uint16)
array_12 = (binned_image/16).astype(np.uint16)
# array_12 =  (array << 4).astype('uint16') #no values above 12 bits

#%%
import numpy as np

#This reader is the only one that works for output from djpeg
def read_pgm(filename):
    with open(filename, 'rb') as file:
        # Read the header (the first four lines)
        magic_number = file.readline().decode().strip()
        if magic_number != 'P5':
            raise ValueError(f"Invalid magic number: {magic_number}. Only P5 format is supported.")
        
        dimensions = file.readline().decode().strip().split()
        width, height = int(dimensions[0]), int(dimensions[1])
        
        max_val = int(file.readline().decode().strip())
        if max_val <= 0 or max_val > 65535:
            raise ValueError(f"Invalid maximum pixel value: {max_val}. Only values up to 65535 are supported for 16-bit PGM.")
        
        # Read the image data
        img_data = np.frombuffer(file.read(), dtype=np.uint16, count=width*height)
        img_data = img_data.reshape((height, width))
        
    return img_data

# def write_12bit_pgm(filename, width, height, pixels):
#     assert max(pixels) <= 4095, "Pixel values should be within the range 0-4095"
#     assert min(pixels) >= 0, "Pixel values should be within the range 0-4095"
    
#     with open(filename, "wb") as f:
#         f.write(b"P5\n")
#         f.write(f"{width} {height}\n".encode())
#         f.write(b"4095\n")
        
#         for pixel in pixels:
#             # Convert 12-bit value to two bytes (big-endian)
#             byte1 = (pixel >> 8) & 0xFF
#             byte2 = pixel & 0xFF
#             f.write(bytes([byte1, byte2]))

def save_16bit_pgm(file_path, image_data):
    # Ensure the input image_data is of type uint16
    if image_data.dtype != np.uint16:
        raise ValueError("Input image_data must be of dtype uint16.")

    # Calculate image width and height
    height, width = image_data.shape

    # Create the PGM header
    header = f'P5\n{width} {height}\n65535\n'

    # Write the header to the file
    with open(file_path, 'wb') as file:
        file.write(header.encode('ascii'))

        # Write the data to the file in little-endian format using struct.pack
        for row in image_data:
            packed_row = struct.pack('<' + 'H' * width, *row)
            file.write(packed_row)


# def read_16bit_pgm_with_pil(filename):
#     # Read the 16-bit PGM image using PIL (Pillow)
#     img_pil = Image.open(filename)

#     img_data = np.array(img_pil, dtype=np.uint16)

#     return img_data

# def read_16bit_png_with_cv2(filename):
#     # Read the 16-bit PNG image using imageio
#     img = imageio.imread(filename)

#     # Convert the image data to a 16-bit NumPy array (dtype=np.uint16)
#     img_data = img.astype(np.uint16)

#     return img_data

#%%
# Save 12 bit image to pgm:

save_16bit_pgm('uncompressed12.pgm',array)

    # %%
array_pgm = read_pgm('uncompressed12.pgm')

np.all(array==array_pgm) #Check that pgm data is equalt to array

#%%
array_12 = array_12.astype(np.int32)
array_pgm = array_pgm.astype(np.int32)
array_100 = read_pgm('decompressed12_100.pgm')
#array_95 = read_pgm('decompressed_95.pgm').astype(np.int32)
array_90 = read_pgm('decompressed12_90.pgm').astype(np.int32)   
#array_85 = read_pgm('decompressed_85.pgm').astype(np.int32)
array_80 = read_pgm('decompressed12_80.pgm').astype(np.int32)

# %%
plt.imshow(array_12)
plt.colorbar()

#%%
def get_diff(array1,array2):
    diff = array1.astype(np.int64)-array2.astype(np.int64)
    fig, axes = plt.subplots(1, 3, figsize=(8, 4))
    axes[0].imshow(array1)
    axes[0].set_title("uncompressed")
    
    axes[1].imshow(array2)
    axes[1].set_title("compressed")
    
    axes[2].imshow(diff)
    axes[2].set_title("diff std: " + str(np.std(diff)).format(".2f"))
    return diff

#%%
get_diff(array_12,array_100)
# %%
get_diff(array_12,array_90)

# %%
get_diff(array_12,array_80)
# %%
#Test other compressions

import numpy as np
import tifffile

def write_12bit_tiff(filename, data):
    # Ensure the data is in the range of 0-4095 (12-bit range)
    data = np.clip(data, 0, 4095)

    # Convert the data to uint16 before writing to the TIFF file
    data = np.uint16(data)

    # Write the TIFF image using tifffile
    tifffile.imwrite(filename, data, dtype='uint16', photometric='minisblack', description='12 bit TIFF')
    # Example usage:
# Assuming 'your_data' is a 2D NumPy array containing your 12-bit pixel data.
write_12bit_tiff('output_12bit.tiff', array_12)


image = Image.fromarray(array_12)
image.save("test.png")
image.save("test.tiff")
# %%
