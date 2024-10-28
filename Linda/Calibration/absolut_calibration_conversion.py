
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
IR2_G=0.0796
UV2_G=0.0163
coeff_relcal=np.array([0.2945*IR2_G, IR2_G, 0.1675*IR2_G, 0.1163*IR2_G, 0.408*UV2_G, UV2_G, 1.])
#print('coefficients based on relative calibration', coeff_relcal)


coeff=coeff_misu
# G is Gumbels calibration 10^-12 counts /(photons nm-1 m-2 s-1 sr-1)
G=1/coeff*1e12
#N: Nickolays star calibration in photons nm-1 m-2 count-1 

omega=2.676e-9 #solid angle in sr

N=omega*G 

coeff_in_ph_per_nm_per_cm2_per_count=N*1e-4

print('coefficients in photons nm-1 cm-2 count-1', coeff_in_ph_per_nm_per_cm2_per_count)

#gumbel=1/ivchenko*1e12*omega*1e-4
star_coeff_in_Nunit=np.array([10.1, 2.99, 21.2, 26.9, 60.9, 23.6, 1.])
print('star calibration in N unit, ie photons nm-1 cm-2 count-1', star_coeff_in_Nunit)
star_coeff_in_Gunit=coeff_star=1/star_coeff_in_Nunit*1e12*1e-4*omega
print('star calibration in Gumbel unit, i.e. 10^-12 counts /(photons nm-1 m-2 sr-1)', star_coeff_in_Gunit)

IR2_G=star_coeff_in_Gunit[1]
UV2_G=star_coeff_in_Gunit[5]
star_and_relcal_in_Gunit=np.array([0.2945*IR2_G, IR2_G, 0.1675*IR2_G, 0.1163*IR2_G, 0.408*UV2_G, UV2_G, 1.])
print('star and relative calibration in Gumbel unit, i.e.  10^-12 counts /(photons nm-1 m-2 sr-1)', star_and_relcal_in_Gunit)
star_and_relcal_in_Nunit=1/star_and_relcal_in_Gunit*1e12*1e-4*omega
print('star and relative calibration in N unit,ie photons nm-1 cm-2 count-1', star_and_relcal_in_Nunit)




# %%
#convert counts to photons
counts=1
TEXPMS_IR2=3000
TEXPMS=TEXPMS_IR2
counts_per_sec=counts/(TEXPMS/1000)
photons_per_sec=counts_per_sec/IR2_G
print('photons per sec for 1 count', photons_per_sec)

# %%
