
#%%
import numpy as np
#coefficientts based on absolute calibration from stars
coeff_abs={'IR4': 26.55111452979997, 'IR1': 10.04057984612763, 'IR3': 21.50032560933404, 'UV1': 54.25400766434724, 'IR2': 2.9694476370436185, 'UV2': 21.6737088994594}
#coefficients based on standard deviation of fitted gaussian  of absolute calibration from stars
coeff_abs_err={'IR4': 2.177442359677955, 'IR1': 0.476611028773807, 'IR3': 1.4290472632142368, 'UV1': 4.288968632317932, 'IR2': 0.09164900776728863, 'UV2': 1.6258329746529314}
#coefficients based on relative calibration from stars
coeff_relcal={'IR2': 1, 'UV2': 1, 'IR1': 3.3474652207795854, 'IR3': 7.117104994509601, 'IR4': 8.988084772972352, 'UV1': 2.495592564313261}
#coefficients based on error of fitted line
coeff_relcal_err={'IR2': 0, 'UV2': 0, 'IR1': 0.01148764948071401, 'IR3': 0.06090434382828135, 'IR4': 0.06010805587114811, 'UV1': 0.01567013605385334}
coeff_combined={}
coeff_combined['IR2']=coeff_abs['IR2']
coeff_combined['UV2']=coeff_abs['UV2']
coeff_combined['IR1']=coeff_abs['IR2']*coeff_relcal['IR1']
coeff_combined['IR3']=coeff_abs['IR2']*coeff_relcal['IR3']
coeff_combined['IR4']=coeff_abs['IR2']*coeff_relcal['IR4']
coeff_combined['UV1']=coeff_abs['UV2']*coeff_relcal['UV1']
coeff_combined_err={}
coeff_combined_err['IR2']=coeff_abs_err['IR2']
coeff_combined_err['UV2']=coeff_abs_err['UV2']


coeff_combined_err['IR1']=(coeff_abs_err['IR2']/coeff_abs['IR2']+coeff_relcal_err['IR1']/coeff_relcal['IR1'])*coeff_combined['IR1']
coeff_combined_err['IR3']=(coeff_abs_err['IR2']/coeff_abs['IR2']+coeff_relcal_err['IR3']/coeff_relcal['IR3'])*coeff_combined['IR3']
coeff_combined_err['IR4']=(coeff_abs_err['IR2']/coeff_abs['IR2']+coeff_relcal_err['IR4']/coeff_relcal['IR4'])*coeff_combined['IR4']
coeff_combined_err['UV1']=(coeff_abs_err['UV2']/coeff_abs['UV2']+coeff_relcal_err['UV1']/coeff_relcal['UV1'])*coeff_combined['UV1']

coeff_combined_relative_error={}
coeff_combined_relative_error['IR1']=coeff_combined_err['IR1']/coeff_combined['IR1']
coeff_combined_relative_error['IR2']=coeff_combined_err['IR2']/coeff_combined['IR2']
coeff_combined_relative_error['IR3']=coeff_combined_err['IR3']/coeff_combined['IR3']
coeff_combined_relative_error['IR4']=coeff_combined_err['IR4']/coeff_combined['IR4']
coeff_combined_relative_error['UV1']=coeff_combined_err['UV1']/coeff_combined['UV1']
coeff_combined_relative_error['UV2']=coeff_combined_err['UV2']/coeff_combined['UV2']




print('************* FINAL RESULT ****************')
print('combined coefficients', coeff_combined)
print('combined coefficients errors', coeff_combined_err)
print('combined coefficients relative errors', coeff_combined_relative_error)
omega=2.676e-9 #solid angle in sr
for key in ['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']:
    print('channel', key)
    print('calibration factor (that is devided with) in code:', 1/(coeff_combined[key])*1e12*1e-4*omega)
    # The calibration factor is given in *10


#%%
print('************* BELOW OLD STUFF ****************')

# #coefficients based udef in version 0.6 of the l1b calibration which is based on 
# # the MISU calibration (as of Jörgs mail 220929 )
coeff_misu = np.array([0.02546, 0.07961, 0.01399, 0.01055, 0.00328, 0.01633, 1.])
# #coefficients based on OHB calibration (as of Jörgs mail 220929 )
coeff_ohb=np.array([0.02544, 0.08144, 0.01330, 0.00973, 0.00289, 0.01479, 1.])
# #coefficients based instrument model (as of Jörgs mail 220929 ) UVchannels qhould be 3/2 times more senistive
# coeff_im=np.array([0.03304, 0.09935, 0.01331, 0.01088, 0.01100*3/2, 0.01164*3/2, 1.])
# #coeffiecients based on relative calibration (as of Jörgs mail 220929 
# IR2_G=0.0796
# UV2_G=0.0163
# coeff_relcal=np.array([0.2945*IR2_G, IR2_G, 0.1675*IR2_G, 0.1163*IR2_G, 0.408*UV2_G, UV2_G, 1.])
# #print('coefficients based on relative calibration', coeff_relcal)



# # coeff is in [G ]Gumbels calibration 10^-12 counts /(photons nm-1 m-2 s-1 sr-1)
a=1/(coeff_ohb*1e-12)
# a in photons nm-1 m-2 s-1 sr-1 per count and is the same as factor a in the paper
print('calibration factor a in  nm-1 m-2 s-1 sr-1 per count:', a)

# #N: Nickolays star calibration unit is photons nm-1 cm-2 count-1 

Nickolayunit=omega*a*1e-4
print('MISU calibration in Nickolays unit unit is photons nm-1 cm-2 count-1', Nickolayunit)
#%%

# coeff_in_ph_per_nm_per_cm2_per_count=N*1e-4

# print('coefficients in photons nm-1 cm-2 count-1', coeff_in_ph_per_nm_per_cm2_per_count)

#gumbel=1/ivchenko*1e12*omega*1e-4
star_coeff_in_Nunit=np.array([10.1, 2.99, 21.2, 26.9, 60.9, 23.6, 1.])
print('star calibration in N unit, ie photons nm-1 cm-2 count-1', star_coeff_in_Nunit)
star_coeff_in_Gunit=(1/star_coeff_in_Nunit)*1e12*1e-4*omega
print('star calibration in Gumbel unit, i.e. 10^-12 counts /(photons nm-1 m-2 sr-1)', star_coeff_in_Gunit)

IR2_G=star_coeff_in_Gunit[1]
UV2_G=star_coeff_in_Gunit[5]
star_and_relcal_in_Gunit=np.array([0.2945*IR2_G, IR2_G, 0.1675*IR2_G, 0.1163*IR2_G, 0.408*UV2_G, UV2_G, 1.])
print('star and relative calibration in Gumbel unit, i.e.  10^-12 counts /(photons nm-1 m-2 sr-1)', star_and_relcal_in_Gunit)
star_and_relcal_in_Nunit=(1/star_and_relcal_in_Gunit)*1e12*1e-4*omega
print('star and relative calibration in N unit,ie photons nm-1 cm-2 count-1', star_and_relcal_in_Nunit)




# # %%
# #convert counts to photons
# counts=1
# TEXPMS_IR2=3000
# TEXPMS=TEXPMS_IR2
# counts_per_sec=counts/(TEXPMS/1000)
# photons_per_sec=counts_per_sec/IR2_G
# print('photons per sec for 1 count', photons_per_sec)

# # %%
# G=7.96e-6 # counts / (\unit{photon\,m^{-2} sr^{-1} nm^{-1}}) 
# print(G, 'counts / (photon m-2 sr-1 nm-1)')
# print(1/G, 'photons m-2 sr-1 nm-1 per count')

# print(1/(G*omega*1e-4), 'photons cm-2 nm-1 per count')

# %%
