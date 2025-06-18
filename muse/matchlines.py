import numpy as np
import matplotlib.pyplot as plt

# wavelength, flux, error, continuum = np.loadtxt('data/quasar.dat', usecols=[0,1,2,3]).transpose()
wavelength, flux, error, continuum = np.loadtxt('data/spec_HE0153-4520_LA.dat', usecols=[0,1,2,3]).transpose()
normalised_flux = flux/continuum
plt.plot(wavelength, normalised_flux, 'k-')

# a,b = 1548.2, 1550.78 # CIV
# a,b = 2796.35, 2803.53 # MgII
# a,b = 3934.78, 3969.59 # CaII
# a,b = 5891.58, 5897.56 # NaI
# a, b = 1238.82, 1242.8 # NV
# a, b = 1393.76, 1402.77 # SiIV
# a, b = 1031.9261, 1037.6167 # OVI
# a, b = 1190.4, 1193.3 # SiII

# plt.plot(wavelength*(a/b), normalised_flux, 'r--')
# plt.plot(wavelength*(b/a), normalised_flux, 'b--')
plt.show()
