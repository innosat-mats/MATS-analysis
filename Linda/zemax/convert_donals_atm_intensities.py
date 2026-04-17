#%%
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
def unit_vector_from_angles(yangle_deg, xangle_deg):
    """
    Generate a 3D unit vector given yangle and xangle.
    
    Parameters
    ----------
    yangle_deg : float
        Angle (degrees) between vector and XZ plane (elevation).
    xangle_deg : float
        Angle (degrees) between vector and YZ plane (azimuth).
    
    Returns
    -------
    (vx, vy, vz) : tuple of floats
        Components of the unit vector.      
        
    """
    # Convert to radians
    yangle = np.radians(yangle_deg)
    xangle = np.radians(xangle_deg)
    
    # Compute components<
    vx = np.cos(yangle) * np.sin(xangle)
    vy = np.sin(yangle)
    vz = np.cos(yangle) * np.cos(xangle)

    # adjust vz so that vz**2+vy**2+vx**2=1, to avoid numerical issues in zemax 
    # round to vx and vy to 7 decimal places
    vx = np.round(vx, 7)
    vy = np.round(vy, 7)
    vz = np.sqrt(1. - vx**2 - vy**2)


    return vx, vy, vz



def adjust_for_curvature(angles, intensities, xangle_deg):
    """
    Adjust intensity for Earth curvature based on xangle.
    
    Parameters
    ----------
    intensity : float
        Original intensity value.
    xangle_deg : float
        X angle in degrees.
    
    Returns
    -------
    adjusted_intensity : float
        Intensity adjusted for Earth curvature.
        
    """
    Earth_radius_km=6371
    sat_height_km=600

    tangent_point_distance_km =np.sqrt((Earth_radius_km + sat_height_km)**2 - Earth_radius_km**2)

    #calculate how many degrees correspond to 1 km altitude at the tangent point
    degrees_per_km = np.degrees(np.arctan(1 / tangent_point_distance_km))

    #Calculate the km that the horizon should be lowered at the given xangle (0 km straight ahead)
    horizon_lowering_km = Earth_radius_km * (1 - np.cos(np.radians(xangle_deg)))

    #Convert horizon lowering to degrees
    horizon_lowering_deg = horizon_lowering_km * degrees_per_km

    # Interpolate shifted intensities
    shifted_angles = angles - horizon_lowering_deg
    interp_func = interp1d(angles, intensities,
                           kind='quadratic', fill_value="extrapolate")
    adjusted_intensities = interp_func(shifted_angles)


    return adjusted_intensities
def mayhit(pos, angle_deg, x_or_y):
    """
    Determine if a ray with given x position and x angle would hit anywhere close to M2 mirror.
    Parameters
    ----------
    xpos : float
        X position of the ray (mm)
    xangle_deg : float
        X angle of the ray (degrees)
    
    Returns
    -------
    bool
        True if the ray would hit the M1 mirror or close to it, False if it would miss the M1 m irror by a large margin.
    """
    if x_or_y == "x":
        M1_radius_mm=64.5/2.+5 #radius of M1 mirror plus some margin from 1704_PRE.4 MATS M1 CAH04
    elif x_or_y == "y":
        M1_radius_mm=39./2.+5 #radius of M1 mirror plus some margin from 1704_PRE.4 MATS M1 CAH04
    else:
        raise ValueError("x_or_y must be 'x' or 'y'")    


    M1_distance_mm=800 #assumes source at distance 800 from M1 

    
    if abs(pos+M1_distance_mm * np.tan(np.radians(angle_deg))) < M1_radius_mm:
        return True
    else:
        return False

     


def generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, Earthcurvature=False, no_hit_cut=True, debug=False, returntotint=False):
    """
    Generate .dat file with atmospheric intensity data.
    
    Parameters
    ----------
    output_dat : str
        Path to output .dat file
    new_angles : array
        Array of angles (degrees)
    new_intensities : array
        Array of intensities
    ypos_range : range
        Y position range
    xpos_range : range
        X position range
    xangle_range : array-like
        X angle range (degrees)

    """
    zpos = 0
    line_count = 0
    # create data frame for xpos, ypos, xangle, yangle, intensity
    df=pd.DataFrame(columns=["xpos", "ypos", "zpos", "xangle_deg", "yangle_deg", "intensity"])

    with open(output_dat, 'w') as f:
        f.write("0 80\n")
        f.write("    ! xpos ypos zpos xcomp ycomp zcomp intensity\n")

        totintensity=0
        hitcut_totintensity=0
        for xpos in xpos_range:
            for ypos in ypos_range:
                for xangle_deg in xangle_range:
                    if Earthcurvature:
                        intensities_adj = adjust_for_curvature(new_angles, new_intensities, xangle_deg) #shift intensities downwars for higher xangles to account for Earth curvature
                    else:
                        intensities_adj = new_intensities
                    
                    for yangle_deg, intensity in zip(new_angles, intensities_adj):
                        intensity_Tera=intensity*1e-12
                        totintensity+=intensity_Tera
                        
                        if no_hit_cut or (mayhit(xpos, xangle_deg, x_or_y="x") and mayhit(ypos, yangle_deg, x_or_y="y")): #only include rays that would hit M1 mirror or come close to it
                            hitcut_totintensity+=intensity_Tera
                            # Save in df
                            if debug: 
                                df = pd.concat([df, pd.DataFrame([{"xpos": xpos, "ypos": ypos, "zpos": zpos, "xangle_deg": xangle_deg, "yangle_deg": yangle_deg, "intensity": intensity_Tera}])], ignore_index=True)

                            [xv, yv, zv] = unit_vector_from_angles(yangle_deg, xangle_deg)
                            f.write(f"{xpos} {ypos} {zpos} {xv:.8f} {yv:.8f} {zv:.8f} {intensity_Tera:.6e}\n")
                            line_count += 1
                            


        
        f.seek(0)
        f.write(f"{line_count} 4\n")

    if debug:
        plot_positions_and_angles_df(df)

    print(f"✅ .dat file '{output_dat}' generated successfully.")
    print(f"📊 The output file contains {line_count} data lines (excluding header).")

    if returntotint:
        return line_count, totintensity, hitcut_totintensity
    else:
        return line_count

def plot_positions_and_angles_df(df):
        
    print(df.head())
    # plot position and angles
    plt.figure(figsize=(10, 5))
    plt.subplot(2, 1, 1)
    plt.scatter(df["xpos"], df["ypos"], c=df["xangle_deg"], cmap='viridis')
    plt.xlabel("X Position (mm)")
    plt.ylabel("Y Position (mm)")
    plt.title("X Position vs Y Position")
    plt.colorbar(label="X Angle (degrees)")

    plt.subplot(2, 1, 2)
    plt.scatter(df["xpos"], df["ypos"], c=df["yangle_deg"], cmap='viridis')
    plt.xlabel("X Position (mm)")
    plt.ylabel("Y Position (mm)")
    plt.title("Y Position vs X Position")
    plt.colorbar(label="Y Angle (degrees)")
    
    plt.tight_layout()
    plt.show()

    # plot histogram of angles and positions
    plt.figure(figsize=(10, 5))
    plt.subplot(2, 2, 1)
    plt.hist(df["xangle_deg"], bins=100)
    plt.xlabel("X Angle (degrees)")
    plt.ylabel("Frequency")
    plt.title("Histogram of X Angles")

    plt.subplot(2, 2, 2)
    plt.hist(df["yangle_deg"], bins=100)
    plt.xlabel("Y Angle (degrees)")
    plt.ylabel("Frequency")
    plt.title("Histogram of Y Angles")

    plt.subplot(2, 2, 3)
    plt.hist(df["xpos"], bins=100)
    plt.xlabel("X Position (mm)")
    plt.ylabel("Frequency")
    plt.title("Histogram of X Positions")

    plt.subplot(2, 2, 4)
    plt.hist(df["ypos"], bins=100)
    plt.xlabel("Y Position (mm)")
    plt.ylabel("Frequency")
    plt.title("Histogram of Y Positions")
    plt.tight_layout()
    plt.show()

    #plot 2d histogram of  positions and angles
    plt.figure(figsize=(8,6))
    plt.hist2d(df["xpos"], df["ypos"], bins=100)
    plt.xlabel("X Position (mm)")
    plt.ylabel("Y Position (mm)")
    plt.title("2D Histogram of X and Y Positions")
    plt.colorbar(label="Frequency")
    plt.show()

    plt.figure(figsize=(8,6))
    plt.hist2d(df["xangle_deg"], df["yangle_deg"], bins=100)
    plt.xlabel("X Angle (degrees)")
    plt.ylabel("Y Angle (degrees)")
    plt.title("2D Histogram of X and Y Angles")
    plt.colorbar(label="Frequency")
    plt.show()



    return



def extrapolate_intensities(angles, intensities, new_angle_min, new_angle_max, plot=False, diff_angles=0.003, const_intensity=False):
    #diffangles of 0.003 is the step size of Donals original data

    # Define new angle grid: use integer indices to avoid floating-point accumulation
    # that can cause np.arange to accidentally include the endpoint (new_angle_max),
    # which would create duplicate/overlapping boundary values across files.
    n_steps = int(round((new_angle_max - new_angle_min) / diff_angles))
    new_angles = new_angle_min + np.arange(n_steps) * diff_angles

    # Logarithmic interpolation with extrapolation
    log_intensities = np.log(intensities)
    interp_func = interp1d(angles, log_intensities,
                           kind='linear', fill_value="extrapolate")
    interp_logI = interp_func(new_angles)
    new_intensities = np.exp(interp_logI)
    if const_intensity:
        new_intensities = np.full_like(new_intensities, 1.)

    if plot:
        plt.figure(figsize=(8,6))
        plt.plot(intensities, angles, 'bo', label="Original Data")
        plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated")
        plt.xlabel("Intensity")
        plt.ylabel("Angle (degrees)")
        plt.title("Intensity vs Angle of Incidence")
        plt.xscale('log')
        plt.legend()
        plt.grid(True)
        plt.show()
    return -new_angles, new_intensities # Return negative angles for Zemax convention (up is negative angle)   

def file_content(dat_file):
    """Print the min/max of xpos, ypos, xangle, and yangle found in a .dat file.

    The file is expected to have a two-token header line followed by an optional
    comment line, then data rows of the form:
        xpos ypos zpos xcomp ycomp zcomp intensity

    Direction cosines (xcomp, ycomp, zcomp) are converted back to angles:
        xangle = degrees(arctan(xcomp / zcomp))
        yangle = degrees(arctan(ycomp / zcomp))
    """
    xpos_list, ypos_list, xangle_list, yangle_list = [], [], [], []

    with open(dat_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            parts = line.split()
            if len(parts) == 2:  # header line "N M"
                continue
            if len(parts) < 7:
                continue
            xp, yp = float(parts[0]), float(parts[1])
            xc, yc, zc = float(parts[3]), float(parts[4]), float(parts[5])
            if zc == 0:
                continue
            xangle_list.append(np.degrees(np.arctan2(xc, zc)))
            yangle_list.append(np.degrees(np.arctan2(yc, zc)))
            xpos_list.append(xp)
            ypos_list.append(yp)

    print(f"File: {dat_file}")
    print(f"  xpos   : min={min(xpos_list):.3f}  max={max(xpos_list):.3f}")
    print(f"  ypos   : min={min(ypos_list):.3f}  max={max(ypos_list):.3f}")
    print(f"  xangle : min={min(xangle_list):.4f}°  max={max(xangle_list):.4f}°")
    print(f"  yangle : min={min(yangle_list):.4f}°  max={max(yangle_list):.4f}°")


def print_yangle_range(dat_file):
    """Print the minimum and maximum yangle found in a .dat file.

    Reads direction cosines (xcomp, ycomp, zcomp) from each data row and
    converts ycomp/zcomp back to yangle = degrees(arctan2(ycomp, zcomp)).
    """
    yangle_list = []

    with open(dat_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            parts = line.split()
            if len(parts) == 2:  # header line "N M"
                continue
            if len(parts) < 7:
                continue
            yc, zc = float(parts[4]), float(parts[5])
            if zc == 0:
                continue
            yangle_list.append(np.degrees(np.arctan2(yc, zc)))

    print(f"File: {dat_file}")
    print(f"  yangle : min={min(yangle_list):.4f}°  max={max(yangle_list):.4f}°")





#%%



# # -----------------------------
# # CONFIGURATION
# # -----------------------------
# input_csv = "/Users/lindamegner/MATS/MATS-retrieval/data/source_intensity_vs_angle_L3.csv"
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/data/zemax_sourcefile_atmosphere_to_-1deg.dat"            # Output DAT file
# ypos_range = range(-40, 41,2)          # Y positions from -40 to 40
# xpos_value = 0                       # X position always zero


# # -----------------------------
# # STEP 1: READ CSV DATA
# # -----------------------------
# angles = []
# intensities = []

# with open(input_csv, 'r') as f:
#     reader = csv.reader(f)
#     for row in reader:
#         angle_deg = float(row[0])
#         intensity = float(row[1])
#         angles.append(angle_deg)
#         intensities.append(intensity)

# angles = np.array(angles)
# intensities = np.array(intensities)



# # -----------------------------
# # STEP 3: EXTRAPOLATE TO +2 DEG
# # -----------------------------
# diff_angles = -np.diff(angles).mean()
# # Define new angle grid from min(filtered) to +2 degrees
# new_angles =np.arange(-2.,-1.,diff_angles)

# # Logarithmic interpolation with extrapolation
# log_intensities = np.log(intensities)
# interp_func = interp1d(angles, log_intensities,
#                        kind='linear', fill_value="extrapolate")
# interp_logI = interp_func(new_angles)
# new_intensities = np.exp(interp_logI)

# # -----------------------------
# # STEP 4: PLOT DATA
# # -----------------------------
# plt.figure(figsize=(8,6))
# plt.plot(intensities, angles, 'bo', label="Original Data")
# plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated (angles > -1.2° to +2°)")
# plt.xlabel("Intensity")
# plt.ylabel("Angle (degrees)")
# plt.title("Intensity vs Angle of Incidence")
# plt.legend()
# plt.grid(True)
# plt.show()

# # %%
# # -----------------------------
# # STEP 3: GENERATE .DAT FILE
# # -----------------------------
# line_count = 0  # counter for number of lines written
# with open(output_dat, 'w') as f:
#     # Write header
#     f.write("xpos ypos zpos  sinang cosang intensity\n")
    
#     # Loop over all angles and intensities
#     for angle_deg, intensity in zip(new_angles, new_intensities):
#         angle_rad = np.radians(angle_deg)  # Convert to radians for sin/cos
#         sinang = np.sin(angle_rad)
#         cosang = np.cos(angle_rad)
        
#         # Repeat for all xpos values
#         for ypos in ypos_range:
#             xpos = 0
#             zpos = 0
#             unsure = 0
#             f.write(f"{xpos} {ypos} {zpos} {unsure} {sinang:.6f} {cosang:.6f} {intensity:.6e}\n")
#             line_count +=1
# print(f"✅ Plot displayed and .dat file '{output_dat}' generated successfully.")
# print(f"📊 The output file contains {line_count} data lines (excluding header).")
# %%

#This generates a 3D unit vector from two angles

import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
import math
import os
# -----------------------------
# CONFIGURATION
# -----------------------------
input_csv = "/Users/lindamegner/MATS/MATS-retrieval/data/source_intensity_vs_angle_L3.csv" # Input CSV file

#yangle_resolution = 0.003 # Step size for extrapolation (degrees)

# -----------------------------
# STEP 1: READ CSV DATA
# -----------------------------
angles = []
intensities = []

with open(input_csv, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        angle_deg = float(row[0])
        intensity = float(row[1])
        angles.append(angle_deg)
        intensities.append(intensity)

angles = np.array(angles)
intensities = np.array(intensities)


#%%
# Generate fields with constant intensity for the angles from which straylight comes from
constant_intensity=True
margin=1
deltapos=6
ypos_range = range(-45-margin, 45+1+margin,deltapos)        
xpos_range = range(-95-margin, 95+1+margin,deltapos)
zpos = 0
delataxangle=0.1
xangle_start=-4
xangle_stop=4.+delataxangle/2.
xangle_range = np.arange(xangle_start, xangle_stop, delataxangle)  # X angles from xangle_start to xangle_stop degrees
deltayangle=0.1

yangle_start=-1.0
yangle_stop=1.0

nohit_cut=True
new_angles, new_intensities = extrapolate_intensities(angles, intensities, yangle_start, yangle_stop, plot=True, diff_angles=deltayangle, const_intensity=constant_intensity)
if nohit_cut:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosconst2Dr2Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.3f}dav{deltayangle:.3f}_pos{deltapos}.dat"
else:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosconst2Dr2Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.3f}dav{deltayangle:.3f}_pos{deltapos}_M1cut.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, no_hit_cut=nohit_cut)
#plot new angles and intensities
plt.figure(figsize=(8,6))
plt.plot(intensities, angles, 'bo', label="Original Data")
plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated to -1deg")
plt.xlabel("Intensity")
plt.ylabel("Angle (degrees)")
plt.title("Intensity vs Angle of Incidence")
plt.legend()
plt.grid(True)
plt.show()

#%%

# Generate fields with constant intensity for ally positions and xangles but for one yangle only
# in order to study how the ghost varies
constant_intensity=True
margin=1
deltapos=5
ypos_range = range(-45-margin, 45+1+margin,deltapos)        
xpos_range = range(-95-margin, 95+1+margin,deltapos)
zpos = 0
delataxangle=0.02
xangle_start=-6
xangle_stop=6.+delataxangle/2.
xangle_range = np.arange(xangle_start, xangle_stop, delataxangle)  # X angles from xangle_start to xangle_stop degrees

deltayangle=0.02
oneangle=-0.3 #oneangle=-0.3
yangle_start=oneangle
yangle_stop=oneangle+deltayangle/2



nohit_cut=True
new_angles, new_intensities = extrapolate_intensities(angles, intensities, yangle_start, yangle_stop, plot=True, diff_angles=deltayangle, const_intensity=constant_intensity)
print(f"  yangle : min={min(new_angles):.4f}°  max={max(new_angles):.4f}°")
print(f"  intensity : min={min(new_intensities):.4e}  max={max(new_intensities):.4e}")


if nohit_cut:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmcons2Dr1Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.3f}dav{deltayangle:.3f}_pos{deltapos}.dat"
else:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmcons2Dr1Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.3f}dav{deltayangle:.3f}_pos{deltapos}_M1cut.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, no_hit_cut=nohit_cut)
#plot new angles and intensities
plt.figure(figsize=(8,6))
plt.plot(intensities, angles, 'bo', label="Original Data")
plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated to -1deg")
plt.xlabel("Intensity")
plt.ylabel("Angle (degrees)")
plt.title("Intensity vs Angle of Incidence")
plt.legend()
plt.grid(True)
plt.show()

file_content(output_dat) #check content of one of the generated files

#%%
# Generate fiels with reduces angular spread, only the angles from which straylight comes from
# zemaxcut (H176|H175|G0)&X_AYZG(3,-5.9)&AXZL(3,1.)&AXZG(3,-1.), so from -1 to +1 in the horizontal direction and up to 0.1 degree 
# below the line of sight in the vertical direction,

margin=1
deltapos=4
ypos_range = range(-45-margin, 45+1+margin,deltapos)        
xpos_range = range(-95-margin, 95+1+margin,deltapos)
zpos = 0
delataxangle=0.075
xangle_start=-5.99
xangle_stop=6.+delataxangle/2.
xangle_range = np.arange(xangle_start, xangle_stop, delataxangle)  # X angles from xangle_start to xangle_stop degrees
deltayangle=0.05

yangle_start=-0.3
yangle_stop=-0.26

nohit_cut=False
new_angles, new_intensities = extrapolate_intensities(angles, intensities, yangle_start, yangle_stop, plot=True, diff_angles=deltayangle)
if nohit_cut:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos2Dr2Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.3f}dav{deltayangle:.3f}_pos{deltapos}.dat"
else:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos2Dr2Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.3f}dav{deltayangle:.3f}_pos{deltapos}_M1cut.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, no_hit_cut=nohit_cut)
#plot new angles and intensities
plt.figure(figsize=(8,6))
plt.plot(intensities, angles, 'bo', label="Original Data")
plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated to -1deg")
plt.xlabel("Intensity")
plt.ylabel("Angle (degrees)")
plt.title("Intensity vs Angle of Incidence")
plt.legend()
plt.grid(True)
plt.show()
#%%
## Generate limited yangles at that time since there is a 1000000 ray maximum in Zemax 
intensityoutput_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/intensity_per_file_output.dat"


margin=1
deltapos=2
ypos_range = range(-45-margin, 45+1+margin,deltapos)
xpos_range = range(-95-margin, 95+1+margin,deltapos)
zpos = 0
delataxangle=0.06
xangle_start=-6.
xangle_stop=6.+delataxangle
xangle_range = np.arange(xangle_start, xangle_stop, delataxangle)  # X angles from xangle_start to xangle_stop degrees
nohit_cut=False

deltayangle=0.06
increment=0.3
filenumber=0
line_count_list=[]
all_yangles = []  # collect all generated y-angles across files
all_yintensities = []  # collect all generated y-intensities across files
yangle_edges = np.round(np.arange(-6., 1. + increment, increment), 10)

fout_intensity = open(intensityoutput_dat, 'w')
fout_intensity.write("filenumber lines totintensity hitcut_totintensity totalintensity_per_ray hitcut_totintensity_per_ray relscalingfactor yangle_start yangle_stop\n")

for i in range(len(yangle_edges) - 1):
    yangle_start = yangle_edges[i]
    yangle_stop = yangle_edges[i + 1]
    filenumber+=1
    new_angles, new_intensities = extrapolate_intensities(angles, intensities, yangle_start, yangle_stop, plot=True, diff_angles=deltayangle)
    
    print(f"Minimum angle in file {filenumber}: {min(new_angles):.3f} degrees, maximum angle: {max(new_angles):.3f} degrees")
    all_yangles.extend(new_angles)
    all_yintensities.extend(new_intensities)
    if nohit_cut:
        output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos_{filenumber:02d}_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.2f}dav{deltayangle:.2f}_pos{deltapos}.dat"
    else:
        output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos_{filenumber:02d}_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.2f}dav{deltayangle:.2f}_pos{deltapos}_M1cut.dat"
    line_count, totintensity, hitcut_totintensity = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, no_hit_cut=nohit_cut, returntotint=True)    
    
    
    print(f"Generated {line_count} rays for yangle range {yangle_start:.2f} to {yangle_stop:.2f} degrees.")
    line_count_list.append(line_count)    
    print(f"First file yangle range: {yangle_start:.2f} to {yangle_stop:.2f} degrees, min(new_angles)={min(new_angles):.3f}, max(new_angles)={max(new_angles):.3f}")
    print(f"Last file yangle range: {yangle_start:.2f} to {yangle_stop:.2f} degrees, min(new_angles)={min(new_angles):.3f}, max(new_angles)={max(new_angles):.3f}")
    print(f"File {filenumber}: lines={line_count}, total intensity={totintensity:.2e}, hitcut total intensity={hitcut_totintensity:.2e}")
    totintensity_per_row= totintensity/line_count if line_count > 0 else 0
    hitcut_totintensity_per_row= hitcut_totintensity/line_count if line_count > 0 else 0
    relscalingfactor = totintensity_per_row / hitcut_totintensity_per_row if hitcut_totintensity_per_row > 0 else 0
    print(f"File {filenumber}: total intensity per ray={totintensity_per_row:.2e}, hitcut total intensity per ray={hitcut_totintensity_per_row:.2e}")
    fout_intensity.write(f"{filenumber} {line_count} {totintensity:.6e} {hitcut_totintensity:.6e} {totintensity_per_row:.6e} {hitcut_totintensity_per_row:.6e} {relscalingfactor:.6e} {yangle_start:.4f} {yangle_stop:.4f}\n")

    if yangle_start >=-0.9: # split into 3 sub-files by simply dividing the arrays into thirds
        n_sub = 3
        splits = np.array_split(np.arange(len(new_angles)), n_sub)
        totline_count = 0
        for suffixcounter, idx in enumerate(splits, start=1):
            new_angles_j = new_angles[idx]
            new_intensities_j = new_intensities[idx]
            if nohit_cut:
                output_dat_j = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos_{filenumber:02d}_{suffixcounter:02d}_yang{new_angles_j[0]:.2f}-{new_angles_j[-1]:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.2f}dav{deltayangle:.2f}_pos{deltapos}.dat"
            else:
                output_dat_j = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos_{filenumber:02d}_{suffixcounter:02d}_yang{new_angles_j[0]:.2f}-{new_angles_j[-1]:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.2f}dav{deltayangle:.2f}_pos{deltapos}_M1cut.dat"
            line_count_j, totintensity_j, hitcut_totintensity_j = generate_dat_file(output_dat_j, new_angles_j, new_intensities_j, ypos_range, xpos_range, xangle_range, no_hit_cut=nohit_cut, returntotint=True)
            totline_count += line_count_j
            print(f"Sub-file {suffixcounter}: {len(new_angles_j)} y-angles, {line_count_j} rays (yangle {new_angles_j[0]:.3f} to {new_angles_j[-1]:.3f})")
            print(f"Sub-file {suffixcounter}: lines={line_count_j}, total intensity={totintensity_j:.2e}, hitcut total intensity={hitcut_totintensity_j:.2e}")
            totintensity_per_row_j = totintensity_j / line_count_j if line_count_j > 0 else 0
            hitcut_totintensity_per_row_j = hitcut_totintensity_j / line_count_j if line_count_j > 0 else 0
            relscalingfactor_j = totintensity_per_row_j / hitcut_totintensity_per_row_j if hitcut_totintensity_per_row_j > 0 else 0
            fout_intensity.write(f"{filenumber}_{suffixcounter} {line_count_j} {totintensity_j:.6e} {hitcut_totintensity_j:.6e} {totintensity_per_row_j:.6e} {hitcut_totintensity_per_row_j:.6e} {relscalingfactor_j:.6e} {new_angles_j[0]:.4f} {new_angles_j[-1]:.4f}\n")
        print(f"Total sub-file rays {totline_count} should equal main file rays {line_count}.")

fout_intensity.close()
print(f"Intensity summary written to: {intensityoutput_dat}")
print(f"number of rays generated per file: {line_count_list}")




# Plot all generated y-angles to check for gaps or jumps
all_yangles = np.array(all_yangles)
all_yintensities = np.array(all_yintensities)
fig, axes = plt.subplots(3, 1, figsize=(10, 10))

axes[0].plot(all_yangles, 'o-', markersize=2)
axes[0].set_xlabel("Index")
axes[0].set_ylabel("Y Angle (degrees)")
axes[0].set_title("All generated y-angles (should be equidistant)")
axes[0].grid(True)

axes[1].plot(np.diff(all_yangles), 'o-', markersize=2)
axes[1].set_xlabel("Index")
axes[1].set_ylabel("Δ Y Angle (degrees)")
axes[1].set_title(f"Step size between consecutive y-angles (should be constant ≈ {deltayangle})")
axes[1].grid(True)

axes[2].plot(all_yintensities, all_yangles, 'o-', markersize=2)
axes[2].set_xlabel("Intensity")
axes[2].set_ylabel("Y Angle (degrees)")
axes[2].set_title("Intensity vs Y Angle across all files")
axes[2].set_xscale('log')
axes[2].grid(True)

plt.tight_layout()
plt.show()

# %%
# Generate file with only the relevant angles that will hit the edge of M2 as determined by the output 
#/Users/lindamegner/MATS/MATS-retrieval/zemax_claude/zemaxfiles/lllatmos2Dr2Da_yang-3.40-0.00xang-6.0-6.1deg_dah0.10dav0.10_pos10_M1cut_w1.ZRD
    #Generate Zemax ray files for a grid of source positions


def get_phi_max(src_x, src_y):
    """Return the maximum phi (degrees) for a source at (src_x, src_y).

    Uses the empirical linear fit:  phi_max = 1.1 + 0.052 * r
    where r = sqrt(src_x² + src_y²) is the radial distance from the optical axis.
    """
    r = math.sqrt(src_x**2 + src_y**2)
    return 1.1 + 0.052 * r


def generate_zemax_rays(src_x, src_y, src_z, phi_max,
                        Dphi=1.0, dphi=0.1):
    """Generate a Zemax source ray file for a single source position.

    Produces rays with phi spanning [phi_max - Dphi, phi_max] in steps of dphi.
    Each ray lies in the plane containing the optical axis (z) and the source
    position vector, directed inward (converging toward the optical axis).

    Direction cosines from phi and radial inward unit vector (-src_x/r, -src_y/r):
        l = -sin(phi) * src_x / r
        m = -sin(phi) * src_y / r
        n =  cos(phi)

    Parameters
    ----------
    src_x, src_y, src_z : float   Source position (mm).
    phi_max             : float   Upper bound of phi range (degrees).
    Dphi                : float   Width of phi range (degrees).  Default 1.0.
    dphi                : float   Step size in phi (degrees).    Default 0.1.
   

    Returns
    -------
    rows : list of tuples  (xpos, ypos, zpos, xcomp, ycomp, zcomp, intensity)
    """
    r = math.sqrt(src_x**2 + src_y**2)
    if r == 0:
        raise ValueError("Source is on the optical axis (r=0); direction is undefined.")

    phistart = max(0.0, phi_max - Dphi)
    phis = np.arange(phistart, phi_max + 1e-9, dphi)

    rows = []
    for phi in phis:
        phi_r = math.radians(phi)
        l = -math.sin(phi_r) * src_x / r
        m = -math.sin(phi_r) * src_y / r
        n =  math.cos(phi_r)
        if m < 0:
          zangle_deg =-math.degrees(math.acos(n)) 
        else:
          zangle_deg = math.degrees(math.acos(n)) 

        intensity=find_intensity(angles, intensities, zangle_deg)
        rows.append((src_x, src_y, src_z, l, m, n, intensity, zangle_deg)) 


    return rows


def generate_zemax_ray_grid(xposrange, yposrange, deltapos, dphi=0.01, src_z=0.0,
                           output_path=None, plotphi=False):
    """Generate a Zemax source ray file for a 2-D grid of source positions.

    For each position the phi range is determined by get_phi_max() and a
    Dphi that depends on the radial distance from the optical axis.
    All rays from every position are written to a single .dat file.

    Parameters
    ----------
    xposrange   : (float, float)  [xmin, xmax] limits for x-positions (mm).
    yposrange   : (float, float)  [ymin, ymax] limits for y-positions (mm).
    deltapos    : int             Grid step size (mm).
    dphi        : float           Angular step in phi (degrees).
    src_z       : float           Source z-coordinate (mm).
    output_path : str             Path for the output .dat file.  Required to save.
    plotphi     : bool            If True, plot mean phi per grid position (default True).

    Returns
    -------
    all_rows : list of tuples  (xpos, ypos, zpos, l, m, n, intensity)
    """
    xpos_range = range(xposrange[0], xposrange[1] + 1, deltapos)
    ypos_range = range(yposrange[0], yposrange[1] + 1, deltapos)
    xy_positions = [(x, y) for y in ypos_range for x in xpos_range]

    all_rows = []
    for src_x, src_y in xy_positions:
        r = math.sqrt(src_x**2 + src_y**2)
        if r == 0:
            print(f"Skipping optical-axis position ({src_x}, {src_y})")
            continue
        phi_max = get_phi_max(src_x, src_y)
        if r > 50:   Dphi = 0.5
        elif r > 20: Dphi = 1.1
        else:        Dphi = 2.0
        all_rows.extend(generate_zemax_rays(
            src_x, src_y, src_z=src_z, phi_max=phi_max,
            Dphi=Dphi, dphi=dphi))

    if output_path:
        with open(output_path, "w") as f:
            f.write(f"{len(all_rows)} 4\n")
            f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
            for xp, yp, zp, lc, mc, nc, inten in all_rows:
                f.write(f"{xp} {yp} {zp} {lc:.8f} {mc:.8f} {nc:.8f} {inten:.6e}\n")
        print(f"Zemax ray file saved to: {output_path}  ({len(all_rows)} rays)")

    if plotphi:
        xs   = np.array([row[0] for row in all_rows])
        ys   = np.array([row[1] for row in all_rows])
        ns   = np.array([row[5] for row in all_rows])   # zcomp = n
        phis = np.degrees(np.arccos(np.clip(ns, -1.0, 1.0)))

        df_phi = pd.DataFrame({"x": xs, "y": ys, "phi": phis})
        mean_phi = df_phi.groupby(["x", "y"])["phi"].mean().reset_index()
        pivot = mean_phi.pivot(index="y", columns="x", values="phi")

        fig, ax = plt.subplots(figsize=(10, 5))
        im = ax.pcolormesh(pivot.columns.values, pivot.index.values,
                           pivot.values, cmap="viridis")
        cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
        cbar.set_label("Mean φ  (degrees)")
        ax.set_xlabel("Source x  (mm)")
        ax.set_ylabel("Source y  (mm)")
        ax.set_title("Mean φ per source position")
        ax.set_aspect("equal")
        plt.tight_layout()

        if output_path:
            plot_path = output_path.replace(".dat", "_mean_phi.png")
            plt.savefig(plot_path, dpi=150)
            print(f"Plot saved to: {plot_path}")
        else:
            plt.show()

    return all_rows


def find_intensity(angles, intensities, zangle_deg):
    """Find the intensity corresponding to a given zangle by interpolation.

    Parameters
    ----------
    angles : array-like  Angles (degrees) corresponding to the intensities.
    intensities : array-like  Intensities corresponding to the angles.
    zangle_deg : float  The vertical angle (degrees) for which to find the intensity.

    Returns
    -------
    intensity : float  The interpolated intensity at zangle_deg.
    """
    interp_func = interp1d(angles, intensities, kind='linear', fill_value=(intensities[0], intensities[-1]), bounds_error=False)
    return interp_func(zangle_deg)



base = os.path.dirname(__file__)

#base = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/"

deltapos = 5
dphi = 0.1
output_path = os.path.join(base, f"zemax_rays_phicut_dpos{deltapos:.2f}_dphi{dphi:.3f}.dat")
generate_zemax_ray_grid(xposrange=[-80, 80], yposrange=[-32, 32], deltapos=deltapos,
                            dphi=dphi, output_path=output_path, plotphi=False)    
        

