
#%%
import numpy as np

#coefficients based udef in version 0.6 of the l1b calibration which is based on 
# the MISU calibration (as of Jörgs mail 220929 )
coeff_misu = np.array([0.02546, 0.07961, 0.01399, 0.01055, 0.00328, 0.01633, 1.])
#coefficients based on OHB calibration (as of Jörgs mail 220929 )
coeff_ohb=np.array([0.02544, 0.08144, 0.01330, 0.00973, 0.00289, 0.01479, 1.])
#coefficients based instrument model (as of Jörgs mail 220929 ) UVchannels qhould be 3/2 times more senistive
coeff_im=np.array([0.03304, 0.09935, 0.01331, 0.01088, 0.01100*3/2, 0.01164*3/2, 1.])
#coeffiecients based on relative calibration (as of Jörgs mail 220929 
IR2=0.0796
UV2=0.0163
coeff_relcal=np.array([0.294*IR2, IR2, 0.160*IR2, 0.115*IR2, 0.41*UV2, UV2, 1.])
#print('coefficients based on relative calibration', coeff_relcal)


coeff=coeff_im
# G is Gumbels calibration in photons nm-1 m-2 count-1 sr-1
G=1/coeff*1e12
#N: Nickolays star calibration in photons nm-1 m-2 count-1 

omega=2.676e-9 #solid angle in sr

N=omega*G 

coeff_in_ph_per_nm_per_cm2_per_count=N*1e-4

print('coefficients in photons nm-1 cm-2 count-1', coeff_in_ph_per_nm_per_cm2_per_count)



# %%
