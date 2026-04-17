#%%
#Throughput according to zemax
ZIR1=0.0907
ZIR2=0.1762
ZIR3=0.0598
ZIR4=0.0674

TotZ=ZIR1+ZIR2+ZIR3+ZIR4

#Throughput according to Jörg's calculations
JIR1=0.17
JIR2=0.22
JIR3=0.075
JIR4=0.060

TotJ=JIR1+JIR2+JIR3+JIR4

relative_throughput_ZIR1=ZIR1/TotZ
relative_throughput_JIR1=JIR1/TotJ
relative_throughput_ZIR2=ZIR2/TotZ
relative_throughput_JIR2=JIR2/TotJ
relative_throughput_ZIR3=ZIR3/TotZ
relative_throughput_JIR3=JIR3/TotJ
relative_throughput_ZIR4=ZIR4/TotZ
relative_throughput_JIR4=JIR4/TotJ

print(f"Throughput according to Zemax: IR1={ZIR1:.4f}, IR2={ZIR2:.4f}, IR3={ZIR3:.4f}, IIR={ZIR4:.4f}")
print(f"Throughput according to Jörg : IR1={JIR1:.4f}, IR2={JIR2:.4f}, IR3={JIR3:.4f}, IR4={JIR4:.4f}")

print(f"Relative throughput according to Zemax: IR1={relative_throughput_ZIR1:.3f}, IR2={relative_throughput_ZIR2:.3f}, IR3={relative_throughput_ZIR3:.3f}, IR4={relative_throughput_ZIR4:.3f}")
print(f"Relative throughput according to Jörg : IR1={relative_throughput_JIR1:.3f}, IR2={relative_throughput_JIR2:.3f}, IR3={relative_throughput_JIR3:.3f}, IR4={relative_throughput_JIR4:.3f}")

print(f"Relative throughput Zemax / Jörg : IR1={relative_throughput_ZIR1/relative_throughput_JIR1:.3f}, IR2={relative_throughput_ZIR2/relative_throughput_JIR2:.3f}, IR3={relative_throughput_ZIR3/relative_throughput_JIR3:.3f}, IR4={relative_throughput_ZIR4/relative_throughput_JIR4:.3f}") 

#Calibration coefficients from Jörgs instrument model
Jcal_IR1=7.95
Jcal_IR2=2.66
Jcal_IR3=19.2
Jcal_IR4=24.0

TotJcal=Jcal_IR1+Jcal_IR2+Jcal_IR3+Jcal_IR4

relative_Jcal_IR1=Jcal_IR1/TotJcal
relative_Jcal_IR2=Jcal_IR2/TotJcal
relative_Jcal_IR3=Jcal_IR3/TotJcal
relative_Jcal_IR4=Jcal_IR4/TotJcal


#Calibration coefficients from star measurements
Scal_IR1=10.0
Scal_IR2=2.97
Scal_IR3=21.5
Scal_IR4=26.6

TotScal=Scal_IR1+Scal_IR2+Scal_IR3+Scal_IR4

relative_Scal_IR1=Scal_IR1/TotScal
relative_Scal_IR2=Scal_IR2/TotScal
relative_Scal_IR3=Scal_IR3/TotScal
relative_Scal_IR4=Scal_IR4/TotScal

print(f"Calibration coefficients according to stars: IR1={Scal_IR1:.2f}, IR2={Scal_IR2:.2f}, IR3={Scal_IR3:.2f}, IR4={Scal_IR4:.2f}")
print(f"Calibration coefficients according to Jörg : IR1={Jcal_IR1:.2f}, IR2={Jcal_IR2:.2f}, IR3={Jcal_IR3:.2f}, IR4={Jcal_IR4:.2f}")
print(f"Relative calibration coefficients according to stars: IR1={relative_Scal_IR1:.3f}, IR2={relative_Scal_IR2:.3f}, IR3={relative_Scal_IR3:.3f}, IR4={relative_Scal_IR4:.3f}")
print(f"Relative calibration coefficients according to Jörg : IR1={relative_Jcal_IR1:.3f}, IR2={relative_Jcal_IR2:.3f}, IR3={relative_Jcal_IR3:.3f}, IR4={relative_Jcal_IR4:.3f}")
print(f"Relative calibration coefficients Jörg / stars : IR1={relative_Jcal_IR1/relative_Scal_IR1:.3f}, IR2={relative_Jcal_IR2/relative_Scal_IR2:.3f}, IR3={relative_Jcal_IR3/relative_Scal_IR3:.3f}, IR4={relative_Jcal_IR4/relative_Scal_IR4:.3f}")   



# %%