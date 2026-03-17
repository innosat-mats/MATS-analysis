#%%
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
        M1_radius_mm=39./2.+5 #radius of M2 mirror plus some margin from 1704_PRE.4 MATS M1 CAH04
    else:
        raise ValueError("x_or_y must be 'x' or 'y'")    


    M1_distance_mm=800 #assumes source at distance 800 from M1 

    
    if abs(pos+M1_distance_mm * np.tan(np.radians(angle_deg))) < M1_radius_mm:
        return True
    else:
        return False

     


def generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, Earthcurvature=False, no_hit_cut=True, debug=False):
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
        
        for xpos in xpos_range:
            for ypos in ypos_range:
                for xangle_deg in xangle_range:
                    if Earthcurvature:
                        intensities_adj = adjust_for_curvature(new_angles, new_intensities, xangle_deg) #shift intensities downwars for higher xangles to account for Earth curvature
                    else:
                        intensities_adj = new_intensities
                    
                    for yangle_deg, intensity in zip(new_angles, intensities_adj):
                        intensity_Tera=intensity*1e-12
                        
                        if no_hit_cut or (mayhit(xpos, xangle_deg, x_or_y="x") and mayhit(ypos, yangle_deg, x_or_y="y")): #only include rays that would hit M1 mirror or come close to it
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



def extrapolate_intensities(angles, intensities, new_angle_min, new_angle_max, plot=False, diff_angles=0.003):
    #diffangles of 0.003 is the step size of Donals original data

    # Define new angle grid from min(filtered) to +2 degrees
    new_angles =np.arange(new_angle_min,new_angle_max,diff_angles)

    # Logarithmic interpolation with extrapolation
    log_intensities = np.log(intensities)
    interp_func = interp1d(angles, log_intensities,
                           kind='linear', fill_value="extrapolate")
    interp_logI = interp_func(new_angles)
    new_intensities = np.exp(interp_logI)

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




# -----------------------------
# STEP 3: GENERATE .DAT FILES
# -----------------------------



#%%
# # Call the function

# margin=0
# ypos_range = range(-45-margin, 45+1+margin,7)        
# xpos_range = range(-95-margin, 95+1+margin,7)
# zpos = 0
# xangle_range = np.arange(-3, 3.1, 0.07)  # X angles from -3 to 3 degrees
# deltayangle=0.07
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -3., -1., plot=True, diff_angles=deltayangle)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos2Dr2Da_to_-1deg_dah0p07dav0p07_pos7.dat"          
# line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, no_hit_cut=True)
# #plot new angles and intensities
# plt.figure(figsize=(8,6))
# plt.plot(intensities, angles, 'bo', label="Original Data")
# plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated to -1deg")
# plt.xlabel("Intensity")
# plt.ylabel("Angle (degrees)")
# plt.title("Intensity vs Angle of Incidence")
# plt.legend()
# plt.grid(True)
# plt.show()

# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -3., 0.66, plot=True, diff_angles=deltayangle)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos2Dr2Da_to_0p66deg_dah0p08dav0p08_pos7.dat"  #dah=delta angle in horizontal direction, dav=delta angle in vertical direction
# line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, no_hit_cut=True)
# ratio=line_count_toplus1/line_count_tominus1
# print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

# #plot_positions_and_angles(output_dat)



#%%
# Generate fiels with reduces angular spread, only the angles from which straylight comes from
# zemaxcut (H176|H175|G0)&X_AYZG(3,-5.9)&AXZL(3,1.)&AXZG(3,-1.), so from -1 to +1 in the horizontal direction and up to 0.1 degree 
# below the line of sight in the vertical direction,

margin=1
deltapos=10
ypos_range = range(-45-margin, 45+1+margin,deltapos)        
xpos_range = range(-95-margin, 95+1+margin,deltapos)
zpos = 0
delataxangle=0.1
xangle_start=-6.
xangle_stop=6.+delataxangle
xangle_range = np.arange(xangle_start, xangle_stop, delataxangle)  # X angles from xangle_start to xangle_stop degrees
deltayangle=0.1

yangle_start=-3.4
yangle_stop=-0.1+deltayangle

nohit_cut=False
new_angles, new_intensities = extrapolate_intensities(angles, intensities, yangle_start, yangle_stop, plot=True, diff_angles=deltayangle)
if nohit_cut:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos2Dr2Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.2f}dav{deltayangle:.2f}_pos{deltapos}.dat"
else:
    output_dat = f"/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmos2Dr2Da_yang{yangle_start:.2f}-{yangle_stop:.2f}xang{xangle_start}-{xangle_stop}deg_dah{delataxangle:.2f}dav{deltayangle:.2f}_pos{deltapos}_M1cut.dat"
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





# #%%
# #Generate files with 2d grid of positions and angles but reduced resolution in angles
# ypos_range = range(-40, 50,20)        
# xpos_range = range(-6*40, 6*50,6*20)
# zpos = 0
# xangle_range = np.arange(-3, 3, 0.15)  # X angles from -3 to 3 degrees
# ycutoffangle = -0.6


# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2, ycutoffangle, plot=True, diff_angles=yangle_resolution)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere2Droom2Dang_to_-0p6deg_reducedresol.dat"
# line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2, 1, plot=True, diff_angles=yangle_resolution)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere2Droom2Dang_to_1deg_reducedresol.dat"  
# line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# ratio=line_count_toplus1/line_count_tominus1
# print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")   

# #%%
# #Generate files with only one ypos
# ypos_range = [0]        
# xpos_range = range(-6*40, 6*50,6*20)
# zpos = 0
# xangle_range = np.arange(-3, 3, 0.1)  # X angles from -3 to 3 degrees
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., -1., plot=True)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dang_to_-1deg.dat"
# line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dang_to_1deg.dat"  
# line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# ratio=line_count_toplus1/line_count_tominus1
# print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

# #%%
# #Generate files with only one xpos but several ypos, no xangle variation
# ypos_range = range(-40, 50,1)       
# xpos_range = [0]
# zpos = 0
# xangle_range = [0]
# ycutoffangle = -0.6
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., ycutoffangle, plot=True, diff_angles=0.08)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom1Dang_to_-0p6deg_dang0p08.dat"
# line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True, diff_angles=0.08)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom1Dang_to_1deg_dang_0p08.dat"  
# line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# ratio=line_count_toplus1/line_count_tominus1
# print(f"Lines to +1deg: {line_count_toplus1}, Lines to {ycutoffangle} deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

# #%%
# #Generate files with only one xpos but several ypos, a few steps in xangle variation
# ypos_range = range(-40, 50,1)       
# xpos_range = [0]
# zpos = 0
# xangle_range = np.arange(-2, 3, 2)  # X angles from -2 to 2 degrees in steps of 1
# ycutoffangle = -0.6
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., ycutoffangle, plot=True)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dangy3_to_-0p6deg.dat"
# line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dangy3_to_1deg.dat"  
# line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# ratio=line_count_toplus1/line_count_tominus1
# print(f"Lines to +1deg: {line_count_toplus1}, Lines to {ycutoffangle} deg: {line_count_tominus1}, Ratio: {ratio:.2f}")




# #%%
# #Generate files with only one ypos and one xpos
# ypos_range = [0]        
# xpos_range = [0]
# zpos = 0
# xangle_range = np.arange(-3, 3, 0.1)  # X angles from -3 to 3 degrees
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., -1., plot=True)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere0Droom2Dang_to_-1deg.dat"
# line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere0Droom2Dang_to_1deg.dat"  
# line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
# ratio=line_count_toplus1/line_count_tominus1
# print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

# #%%



# # %%
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/single_ray_fan_0p01deg.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     xpos = 0
#     ypos = 0
#     zpos = 0
#     intensity = 1.
#     line_count = 0
#     for yangle in np.arange(-2, 2, 0.01):
#         xangle=yangle  
#         [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
        
#         f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#         line_count += 1
#     f.write(f"{line_count} 4\n")

# # %%
# # ray fan with 0.1 degree steps in x and y directions and denser steps along lower diagonal
# # for several xpos and ypos positions
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest3.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     intensity = 1.
#     line_count = 0
#     for xpos in range(-40*6, 40*6+1, 20):
#         for ypos in range(-40, 41, 20):
#             for xangle in np.arange(-4, 4, 0.1):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
#             for yangle in np.arange(-1.5, 1.5, 0.1):
#                 xangle=0 
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
#             for yangle in np.arange(-1.5, 0, 0.1):
#                 xangle=yangle  
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
#             for yangle in np.arange(0, 1.5, 0.05): 
#                 xangle=yangle  
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")
    
# # %%
# # ray fan with 0.1 degree steps in y directions and 0,2 degree steps in x direction 
# # for several xpos and ypos positions
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest4.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     intensity = 1.
#     line_count = 0
#     for xpos in range(-40*6, 40*6+1, 20):
#         for ypos in range(-40, 41, 20):
#             for xangle in np.arange(-4, 4, 0.2):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
#             for yangle in np.arange(-1.5, 1.5, 0.1):
#                 xangle=0 
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
#             for yangle in np.arange(-1.5, 0, 0.1):
#                 xangle=yangle  
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
#             for yangle in np.arange(0, 1.5, 0.1): 
#                 xangle=yangle  
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")
    
# # %%
# # ray fan with  1 degree steps in x direction 
# # for several xpos and ypos positions
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest5.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     intensity = 1.
#     line_count = 0
#     for xpos in range(-40*6, 40*6+1, 20):
#         for ypos in range(-20, 21, 20):
#             for xangle in np.arange(-4, 4, 1):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")


# # %%
# # ray fan with  1 degree steps in x direction 
# # for ONE xpos and ypos positions
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest7.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     ypos=0

#     intensity = 1.
#     line_count = 0
#     for xpos in range(-20*6, 20*6+1, 20):
#             for xangle in np.arange(-3, 1, 1):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")

# # ray fan with  1 degree steps in x direction 
# # for ONE xpos and ypos positions
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest8.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     xpos=0  
#     ypos=0
#     intensity = 1.
#     line_count = 0
#     for hej in range(-20*6, 20*6+1, 20):
#             for xangle in np.arange(-3, 1, 1):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")



# # ray fan with  0.01 degree steps in x direction 
# # for ONE xpos and ypos positions
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest9.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     xpos=0  
#     ypos=0
#     intensity = 1.
#     line_count = 0
#     for hej in range(-20*6, 20*6+1, 20):
#             for xangle in np.arange(-3, 1, 0.01):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")
# # %%
# # ray fan with  0.1 degree steps in x direction 
# # for several xpos
# output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest10.dat" 
# with open(output_dat, 'w') as f:
#     f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
#     zpos = 0
#     xpos=0  
#     ypos=0
#     xpos_range = range(-6*40, 6*40+1,10)
#     intensity = 1.
#     line_count = 0
#     for xpos in xpos_range:
#             for xangle in np.arange(-3, 3, 0.1):
#                 yangle=0
#                 [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
#                 f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
#                 line_count += 1
            
#     f.write(f"{line_count} 4\n")
#     print(f"📊 The output file contains {line_count} data lines (excluding header)")
# %%
# #unit vector from angles function test
# yangle=0.05
# xangle=3.
# [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)



#0.07415180 0.05930640 0.99548192
xv=0.07415180
yv=0.05930640
zv=0.99548192

#print(f"Unit vector for yangle={yangle} and xangle={xangle}: ({xv:.6f}, {yv:.6f}, {zv:.6f})")
sum_of_squares = xv**2 + yv**2 + zv**2
print(f"Sum of squares: {sum_of_squares:.6f} (should be close to 1 for a unit vector)")
    # %%

# %%
