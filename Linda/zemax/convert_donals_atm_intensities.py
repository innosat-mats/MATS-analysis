#%%
import numpy as np
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
    
    # Compute components
    vx = np.cos(yangle) * np.sin(xangle)
    vy = np.sin(yangle)
    vz = np.cos(yangle) * np.cos(xangle)
    
    return vx, vy, vz


def generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range, Earthcurvature=False):
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
    
    with open(output_dat, 'w') as f:
        f.write("0 80\n")
        f.write("    ! xpos ypos zpos xcomp ycomp zcomp intensity\n")
        
        for xpos in xpos_range:
            for ypos in ypos_range:
                for yangle_deg, intensity in zip(new_angles, new_intensities):
                    intensity_Tera=intensity*1e-12
                    for xangle_deg in xangle_range:
                        [xv, yv, zv] = unit_vector_from_angles(yangle_deg, xangle_deg)
                        f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity_Tera:.6e}\n")
                        line_count += 1
        
        f.seek(0)
        f.write(f"{line_count} 4\n")
    
    print(f"✅ .dat file '{output_dat}' generated successfully.")
    print(f"📊 The output file contains {line_count} data lines (excluding header).")

    return line_count


def extrapolate_intensities(angles, intensities, new_angle_min, new_angle_max, plot=False,red_resol_fact=1.):
    diff_angles = -np.diff(angles).mean()*red_resol_fact # Determine and Adjust resolution if needed
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

# -----------------------------
# CONFIGURATION
# -----------------------------
input_csv = "/Users/lindamegner/MATS/MATS-retrieval/data/source_intensity_vs_angle_L3.csv"
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/data/zemax_sourcefile_atmosphere_to_-1deg.dat"            # Output DAT file
ypos_range = range(-40, 41,2)          # Y positions from -40 to 40
xpos_value = 0                       # X position always zero


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
# STEP 3: EXTRAPOLATE TO +2 DEG
# -----------------------------
diff_angles = -np.diff(angles).mean()
# Define new angle grid from min(filtered) to +2 degrees
new_angles =np.arange(-2.,-1.,diff_angles)

# Logarithmic interpolation with extrapolation
log_intensities = np.log(intensities)
interp_func = interp1d(angles, log_intensities,
                       kind='linear', fill_value="extrapolate")
interp_logI = interp_func(new_angles)
new_intensities = np.exp(interp_logI)

# -----------------------------
# STEP 4: PLOT DATA
# -----------------------------
plt.figure(figsize=(8,6))
plt.plot(intensities, angles, 'bo', label="Original Data")
plt.plot(new_intensities, new_angles, 'r-', label="Extrapolated (angles > -1.2° to +2°)")
plt.xlabel("Intensity")
plt.ylabel("Angle (degrees)")
plt.title("Intensity vs Angle of Incidence")
plt.legend()
plt.grid(True)
plt.show()

# %%
# -----------------------------
# STEP 3: GENERATE .DAT FILE
# -----------------------------
line_count = 0  # counter for number of lines written
with open(output_dat, 'w') as f:
    # Write header
    f.write("xpos ypos zpos  sinang cosang intensity\n")
    
    # Loop over all angles and intensities
    for angle_deg, intensity in zip(new_angles, new_intensities):
        angle_rad = np.radians(angle_deg)  # Convert to radians for sin/cos
        sinang = np.sin(angle_rad)
        cosang = np.cos(angle_rad)
        
        # Repeat for all xpos values
        for ypos in ypos_range:
            xpos = 0
            zpos = 0
            unsure = 0
            f.write(f"{xpos} {ypos} {zpos} {unsure} {sinang:.6f} {cosang:.6f} {intensity:.6e}\n")
            line_count +=1
print(f"✅ Plot displayed and .dat file '{output_dat}' generated successfully.")
print(f"📊 The output file contains {line_count} data lines (excluding header).")
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
# Call the function
ypos_range = range(-40, 50,20)        
xpos_range = range(-6*40, 6*50,6*20)
zpos = 0
xangle_range = np.arange(-3, 3, 0.36)  # X angles from -3 to 3 degrees
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2, -1, plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere2Droom2Dang_to_-1deg.dat"          
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
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

new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2, 1, plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere2Droom2Dang_to_1deg.dat"  
line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
ratio=line_count_toplus1/line_count_tominus1
print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

#%%
#Generate files with 2d grid of positions and angles but reduced resolution in angles
ypos_range = range(-40, 50,20)        
xpos_range = range(-6*40, 6*50,6*20)
zpos = 0
xangle_range = np.arange(-3, 3, 0.15)  # X angles from -3 to 3 degrees
ycutoffangle = -0.6
reduced_resoltion_factor=3.
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2, ycutoffangle, plot=True,red_resol_fact=reduced_resoltion_factor)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere2Droom2Dang_to_-0p6deg_reducedresol.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2, 1, plot=True,red_resol_fact=reduced_resoltion_factor)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere2Droom2Dang_to_1deg_reducedresol.dat"  
line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
ratio=line_count_toplus1/line_count_tominus1
print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")   

#%%
#Generate files with only one ypos
ypos_range = [0]        
xpos_range = range(-6*40, 6*50,6*20)
zpos = 0
xangle_range = np.arange(-3, 3, 0.1)  # X angles from -3 to 3 degrees
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., -1., plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dang_to_-1deg.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dang_to_1deg.dat"  
line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
ratio=line_count_toplus1/line_count_tominus1
print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

#%%
#Generate files with only one xpos but several ypos, no xangle variation
ypos_range = range(-40, 50,1)       
xpos_range = [0]
zpos = 0
xangle_range = [0]
ycutoffangle = -0.6
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., ycutoffangle, plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom1Dang_to_-0p6deg.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom1Dang_to_1deg.dat"  
line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
ratio=line_count_toplus1/line_count_tominus1
print(f"Lines to +1deg: {line_count_toplus1}, Lines to {ycutoffangle} deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

#%%
#Generate files with only one xpos but several ypos, a few steps in xangle variation
ypos_range = range(-40, 50,1)       
xpos_range = [0]
zpos = 0
xangle_range = np.arange(-2, 3, 2)  # X angles from -2 to 2 degrees in steps of 1
ycutoffangle = -0.6
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., ycutoffangle, plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dangy3_to_-0p6deg.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere1Droom2Dangy3_to_1deg.dat"  
line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
ratio=line_count_toplus1/line_count_tominus1
print(f"Lines to +1deg: {line_count_toplus1}, Lines to {ycutoffangle} deg: {line_count_tominus1}, Ratio: {ratio:.2f}")




#%%
#Generate files with only one ypos and one xpos
ypos_range = [0]        
xpos_range = [0]
zpos = 0
xangle_range = np.arange(-3, 3, 0.1)  # X angles from -3 to 3 degrees
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., -1., plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere0Droom2Dang_to_-1deg.dat"
line_count_tominus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
new_angles, new_intensities = extrapolate_intensities(angles, intensities, -2., 1., plot=True)
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/atmosphere0Droom2Dang_to_1deg.dat"  
line_count_toplus1 = generate_dat_file(output_dat, new_angles, new_intensities, ypos_range, xpos_range, xangle_range)
ratio=line_count_toplus1/line_count_tominus1
print(f"Lines to +1deg: {line_count_toplus1}, Lines to -1deg: {line_count_tominus1}, Ratio: {ratio:.2f}")

#%%



# %%
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/single_ray_fan_0p01deg.dat" 
with open(output_dat, 'w') as f:
    f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
    xpos = 0
    ypos = 0
    zpos = 0
    intensity = 1.
    line_count = 0
    for yangle in np.arange(-2, 2, 0.01):
        xangle=yangle  
        [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
        
        f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
        line_count += 1
    f.write(f"{line_count} 4\n")

# %%
# ray fan with 0.1 degree steps in x and y directions and denser steps along lower diagonal
# for several xpos and ypos positions
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest3.dat" 
with open(output_dat, 'w') as f:
    f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
    zpos = 0
    intensity = 1.
    line_count = 0
    for xpos in range(-40*6, 40*6+1, 20):
        for ypos in range(-40, 41, 20):
            for xangle in np.arange(-4, 4, 0.1):
                yangle=0
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            for yangle in np.arange(-1.5, 1.5, 0.1):
                xangle=0 
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            for yangle in np.arange(-1.5, 0, 0.1):
                xangle=yangle  
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            for yangle in np.arange(0, 1.5, 0.05): 
                xangle=yangle  
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            
    f.write(f"{line_count} 4\n")
    print(f"📊 The output file contains {line_count} data lines (excluding header)")
    
# %%
# ray fan with 0.1 degree steps in y directions and 0,2 degree steps in x direction 
# for several xpos and ypos positions
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest4.dat" 
with open(output_dat, 'w') as f:
    f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
    zpos = 0
    intensity = 1.
    line_count = 0
    for xpos in range(-40*6, 40*6+1, 20):
        for ypos in range(-40, 41, 20):
            for xangle in np.arange(-4, 4, 0.2):
                yangle=0
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            for yangle in np.arange(-1.5, 1.5, 0.1):
                xangle=0 
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            for yangle in np.arange(-1.5, 0, 0.1):
                xangle=yangle  
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            for yangle in np.arange(0, 1.5, 0.1): 
                xangle=yangle  
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            
    f.write(f"{line_count} 4\n")
    print(f"📊 The output file contains {line_count} data lines (excluding header)")
    
# %%
# ray fan with  1 degree steps in x direction 
# for several xpos and ypos positions
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest5.dat" 
with open(output_dat, 'w') as f:
    f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
    zpos = 0
    intensity = 1.
    line_count = 0
    for xpos in range(-40*6, 40*6+1, 20):
        for ypos in range(-20, 21, 20):
            for xangle in np.arange(-4, 4, 1):
                yangle=0
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            
    f.write(f"{line_count} 4\n")
    print(f"📊 The output file contains {line_count} data lines (excluding header)")


# %%
# ray fan with  1 degree steps in x direction 
# for ONE xpos and ypos positions
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest7.dat" 
with open(output_dat, 'w') as f:
    f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
    zpos = 0
    ypos=0

    intensity = 1.
    line_count = 0
    for xpos in range(-20*6, 20*6+1, 20):
            for xangle in np.arange(-3, 1, 1):
                yangle=0
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            
    f.write(f"{line_count} 4\n")
    print(f"📊 The output file contains {line_count} data lines (excluding header)")

# ray fan with  1 degree steps in x direction 
# for ONE xpos and ypos positions
output_dat = "/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/raytest8.dat" 
with open(output_dat, 'w') as f:
    f.write("! xpos ypos zpos xcomp ycomp zcomp intensity\n")
    zpos = 0
    xpos=0  
    ypos=0
    intensity = 1.
    line_count = 0
    for hej in range(-20*6, 20*6+1, 20):
            for xangle in np.arange(-3, 1, 1):
                yangle=0
                [xv, yv, zv] = unit_vector_from_angles(yangle, xangle)
                f.write(f"{xpos} {ypos} {zpos} {xv:.6f} {yv:.6f} {zv:.6f} {intensity:.6e}\n")
                line_count += 1
            
    f.write(f"{line_count} 4\n")
    print(f"📊 The output file contains {line_count} data lines (excluding header)")

# %%
