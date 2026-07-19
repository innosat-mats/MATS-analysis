# Script to calculate what angles of light that can reach onto the nadir cameras baffle surface
# Note that this is not entirrly correct since the nadir camera is pointed forward, I cannot convince myself if that actually makse a difference though.


#%%

import numpy as np

def tan_deg(x):
    return np.tan(np.radians(x))
def sin_deg(x):
    return np.sin(np.radians(x))
def cos_deg(x):
    return np.cos(np.radians(x))


phi = np.degrees(np.arctan(8.1/12.1))  # Baffle cut at this angle (degrees)

print('Baffle cut angle:', phi, 'degrees')
print('Tan of baffle cut angle:', tan_deg(phi))


#%%

beta= 40 # Azimuth angle of incoming light, zero is in limb forward pointing directiion

gamma=np.degrees(np.arctan(sin_deg(beta)*tan_deg(phi))) # Angle of incoming light that can reach the baffle surface


print('Azimuth angle of incoming light, beta:', beta, 'degrees')
print('gamma :', gamma, 'degrees below horizon')




# %%
