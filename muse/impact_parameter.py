import numpy as np
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u

def angular_separation(ra1, dec1, ra2, dec2): # RA and Dec in degrees
    x = np.cos(dec1) * np.sin(dec2) - np.sin(dec1) * np.cos(dec2) * np.cos(ra2 - ra1)
    y = np.cos(dec2) * np.sin(ra2 - ra1)
    z = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra2 - ra1)

    phi = np.arctan2(np.sqrt(x*x + y*y), z)

    return np.degrees(phi)

class Coordinate:
    def __init__(self, RA, Dec):
        self.RA = RA
        self.Dec = Dec
        self.initialise()

    def initialise(self):
        self.RA['degrees'] = (self.RA['h']*15) + (self.RA['m']/4) + (self.RA['s']/240)
        self.Dec['degrees'] = self.Dec['d'] + self.Dec['m']/60 + self.Dec['s']/3600

    def degrees(self):
        return [self.RA['degrees'], self.Dec['degrees']]
    
    def radians(self):
        return np.deg2rad([self.RA['degrees'], self.Dec['degrees']])

p1 = Coordinate(RA={'h':3, 'm':36, 's':27}, Dec={'d':-20, 'm': 19, 's': 39})
p2 = Coordinate(RA={'h':3, 'm':36, 's':25}, Dec={'d':-20, 'm': 19, 's': 31.9})

phi = angular_separation(*p1.radians(), *p2.radians()) # in degrees
print(f'Angular separation: {phi:.2e} deg = {phi*60:.2e} arcmin')

z = 0.502
d = cosmo.kpc_proper_per_arcmin(z)*phi*60*u.arcmin
print(f'Impact parameter  : {d:.1f}')