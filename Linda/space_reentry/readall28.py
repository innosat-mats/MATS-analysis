
# %%

"""
In CARMA the smoke is in represented in number density per cm3 in the variable
pc(nx,ny,nz,nbin,nelement) (particle concentration) whereas in CHEM2D both mixing ratio XC and number
density X is occur. However XC is the one that is most trustworthy (according
to Dave) and is therefore the one that is used. The bins XC(:,:,NCON-NDUST+1:NCON)
i.e. generally with NCON=79 and NDUST=28 XC(:,:,52:79) are used for smoke.

Beware that XC contains other stuff too! 
XC(nlat,nalt,1)= potential temp 
XC(nlat,nalt,2)= zonal mom
XC(nlat,nalt,3-7)= blank
XC(nlat,nalt,8-51)= gases/chemicals
XC(nlat,nalt,52-79)= smoke of various size bins



The following things are, and should remain set to the same values in the two
models.

Density of smoke: 
	in Carma in setupaer.f rhoelem(1)=2. g/cm3
	in Chem2d in setsmoke.f RHODUST=2000. kg/cm3

Smoke radius bins:
	in Carma in setupaer.f rmin(1)=2.2e-7 cm
		    rmrat=2. volume rate between neighbouring bins
	         in aerad.h NBIN=28

"""
#%%
import numpy as np
import matplotlib.pyplot as plt
ndays=366
outdata = np.zeros((37, 44, 29, ndays))
days = np.zeros(ndays)
years = np.zeros(ndays)
z = np.zeros(100)
p = np.zeros(100)
lat = np.zeros(100)
outxc = np.zeros((37, 44, 29, ndays))
outnd = np.zeros((37, 44, 29, ndays))
dens = np.zeros(100)

datafields = 100
count = np.zeros(datafields, dtype=int)
counter = 0
counterr = 0
#%%


filename='/Users/lindamegner/StuffFromOldComputers/redpc/CC2D/chem2d/dat/two.000'
fid = open(filename, "r")

for i in range(100000):
    fieldr = np.fromfile(fid, dtype=np.float32, count=1)

    if fieldr.size == 0:
        fid.close()
        if days[0] == days[1]: #th first value has been read in twice
            days = days[1:]
            outxc = outxc[:, :, :, 1:]
            outnd = outnd[:, :, :, 1:]
        break
    fieldr = fieldr - 1
    field=int(fieldr)
    day = int(np.fromfile(fid, dtype=np.float32, count=1))
    year = int(np.fromfile(fid, dtype=np.float32, count=1))
    nolat = int(np.fromfile(fid, dtype=np.float32, count=1))
    noz = int(np.fromfile(fid, dtype=np.float32, count=1))




    tmplat = np.fromfile(fid, dtype=np.float32, count=nolat)
    tmpz = np.fromfile(fid, dtype=np.float32, count=noz)
    data = np.fromfile(fid, dtype=np.float32, count=nolat * noz)
    data = data.reshape((noz, nolat))

    # if field == 8:
    #     counter += 1
    #     temp[:, :, counter] = data
    #     z_short = tmpz
    #     lat_short = tmplat
    #     daystemp[counter] = day

    if (field >= 251 and field <= 278) or (field >= 281 and field <= 308):
        fieldnr = int(field - 250)
        count[fieldnr] += 1
        days[count[fieldnr]] = day
        years[count[fieldnr]] = year
        z = tmpz
        lat = tmplat
        p = 1000 * np.exp(-tmpz / 7)
        if fieldnr <= 28:
            outnd[:, :, fieldnr, count[fieldnr]] = data.T
        else:
            outxc[:, :, fieldnr - 30, count[fieldnr]] = data.T


# %%
# Outnd is the smoke particle number density in cm-3 for diffent radii.
# It is a 4D array with dimensions (lat, z, smokeradius, day)
# Outxc holds other quantities, such as the h2o mixing ratio

#This should be the smoke radius size used - check this [ in m] 
r =1.e-2* np.array([2.000E-08, 2.520E-08, 3.175E-08, 4.000E-08, 5.040E-08, 6.350E-08, 8.000E-08, 1.008E-07, 1.270E-07, 1.600E-07, 2.016E-07, 2.540E-07, 3.200E-07, 4.032E-07, 5.080E-07, 6.400E-07, 8.063E-07, 1.016E-06, 1.280E-06, 1.613E-06, 2.032E-06, 2.560E-06, 3.225E-06, 4.064E-06, 5.120E-06, 6.451E-06, 8.127E-06, 1.024E-05])

#Plot the smallest smoke partice radius
plt.figure()
plt.pcolormesh(lat, z, outnd[:, :, 1, 1].T)
plt.colorbar()
plt.xlabel('Latitude')
plt.ylabel('Altitude [km]')


figure, axes = plt.subplots(2, 2, figsize=(10, 5))
irad = 1 #smallest smoke particle radius
iday=30
im=axes[0,0].pcolormesh(lat, z, outnd[:, :, irad, iday].T)
axes[0,0].set_xlabel('Latitude')
axes[0,0].set_ylabel('Altitude [km]')
axes[0,0].set_title('Smoke Number Density [cm-3], DAY = '+str(days[iday])+ 'radius = '+str(r[irad])) 
figure.colorbar(im, ax=axes[0,0])

irad=10 
im=axes[0,1].pcolormesh(lat, z, outnd[:, :, irad, iday].T)
axes[0,1].set_xlabel('Latitude')
axes[0,1].set_ylabel('Altitude [km]')
axes[0,1].set_title('Smoke Number Density [cm-3], DAY = '+str(days[iday])+ 'radius = '+str(r[irad]))
figure.colorbar(im, ax=axes[0,1])

irad=20
im=axes[1,0].pcolormesh(lat, z, outnd[:, :, irad, iday].T)
axes[1,0].set_xlabel('Latitude')
axes[1,0].set_ylabel('Altitude [km]')
axes[1,0].set_title('Smoke Number Density [cm-3], DAY = '+str(days[iday])+ 'radius = '+str(r[irad]))
figure.colorbar(im, ax=axes[1,0])

irad=27
im=axes[1,1].pcolormesh(lat, z, outnd[:, :, irad, iday].T)
axes[1,1].set_xlabel('Latitude')
axes[1,1].set_ylabel('Altitude [km]')
axes[1,1].set_title('Smoke Number Density [cm-3], DAY = '+str(days[iday])+ 'radius = '+str(r[irad]))
figure.colorbar(im, ax=axes[1,1])
figure.tight_layout()

plt.tight_layout()

#%%


# %%
# Rayleigh scattering cross section:
lambda0=270e-9 # wavelength of laser [m]



# Ideal gas law: pV=nRT
p0=101325 # Pa
p=p0*10**(-z/16.) # pressure at altitudes [Pa]


# Define the altitudes (in km)
alts = 1e-3*np.array([-1000, 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000])
# Define the densities (in kg/m^3)
densities = np.array([1.347, 1.225, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.5900, 0.5258, 0.4671, 0.4135, 0.1948, 0.08891, 0.04008, 0.01841, 0.003996, 0.001027, 0.0003097, 0.00008283, 0.00001846])
# Interpolate the densities at the altitudes given by z
interpolated_densities = np.interp(z, alts, densities)

number_densities_air = interpolated_densities*6.022e23/28.97 # number density of air molecules [m^-3]




#ND_air = 1e6*np.array([1e15, 1e13, 1e11])
smoke=1e6*outnd # smoke number density [m^-3]
r_air=0.2e-9 # approximate radius of N2 molecule [m]
n=1
nfactor=(n**2-1)**2/((n**2+2)**2)

outdata=np.zeros([len(lat), len(z), ndays])
iday=2 #day of the year
ilat=1 #latitude
for iday in range(0, ndays):
    for ilat in range(0,len(lat)):
        for iz in range(0, len(z)): 
            #print('alt in  km =', z[iz])
            #Caluclate smoke scattering cross section
            sigma=0
            for irad, radius in enumerate(r):
                sigma_part = smoke[ilat, iz, irad, iday]*(2*np.pi**5/3)*2*(r[irad]**6)/lambda0**4
                sigma=sigma+sigma_part
            sigmaref = number_densities_air[iz]*(2*np.pi**5/3)*2*(r[irad]**6)/lambda0**4

                
            sigmarel=sigma/sigmaref
            #print('alt=',  z[iz], 'lat=', lat[ilat], 'sigmarel=', sigmarel)  
            outdata[ilat, iz, iday]=sigmarel




# %%
import matplotlib.colors as colors
norm = colors.LogNorm()  # Logarithmic colorscale
figure, axes = plt.subplots(2, 2, figsize=(10, 5))


iday=5
im=axes[0,0].pcolormesh(lat, z, (outdata[:, :, iday].T), norm=norm)
axes[0,0].set_xlabel('Latitude')
axes[0,0].set_ylabel('Altitude [km]')
axes[0,0].set_title('Sigmarel, DAY = '+str(days[iday])) 
figure.colorbar(im, ax=axes[0,0])



print('max sigmarel=', outdata.max())

# %%
