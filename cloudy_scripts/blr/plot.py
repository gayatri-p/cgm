import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import pandas as pd
from cloudy_fit_lib import *
from scipy.interpolate import RegularGridInterpolator
import pickle
import scienceplots

rootdir = '/usr3/project/gayatri.p/Documents/cloudy_scripts/blr/grids/'
plt.style.use(['science', 'notebook', 'no-latex'])
plt.rcParams.update({
    # "font.family": "serif",   # specify font family here
    # "font.serif": ["Times"],  # specify font here
    "font.size":12}) 

hm12_z_grid, hm12_wav_grid, hm12_J_nu_grid = read_uvb('', 'hm12_galaxy.ascii')
z_test = 0.226
hm12_J_nu_test = fetch_sed(z_test, hm12_z_grid, hm12_J_nu_grid)

# [X, Y] = np.meshgrid(hm12_wav_grid, hm12_z_grid)
# fig, ax = plt.subplots(1, 1, figsize=(10, 6))
# img = plt.contourf(X, Y, hm12_J_nu_grid)
# plt.xlabel(r'Wavelength [$\AA$]')
# plt.ylabel(r'Redshift [z]')
# plt.title(r'UVB spectrum based on HM12 ($J_\nu$)')
# plt.colorbar(img)
# plt.show()

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# plt.plot(hm12_wav_grid, hm12_J_nu_test)
# plt.xlabel(r'Wavelength [$\AA$]')
# plt.ylabel(r'$J_\nu$ [ergs/s/cm$^2$/Hz/sr]')
# plt.title(f'UVB spectrum based on HM12 (at $z = {z_test}$)')
# plt.show()

## plotting
logN_HI_min = 12
logN_HI_max = 12.25
logN_HI_step = 0.25

logN_HI_arr = np.array([15, 16]) #np.arange(logN_HI_min, logN_HI_max+logN_HI_step, logN_HI_step)

file_list = create_grid_file_list(logN_HI_arr)

# Specify grid points for densities and metallicities
log_hdens_min = -4.25
log_hdens_max = 1
log_hdens_step = 0.25
log_hdens_arr = np.arange(log_hdens_min, log_hdens_max+log_hdens_step, log_hdens_step)
logU_arr = np.log10(calc_U(hm12_wav_grid, hm12_J_nu_test, 10**log_hdens_arr))

log_metals_min = -0.25
log_metals_max = 1
log_metals_step = 0.25
log_metals_arr = np.arange(log_metals_min, log_metals_max+log_metals_step, log_metals_step)

# Fetch $\log N$ for each species, and iterpolate the grid points.
pkl_file = open(rootdir+'final_grid.pkl', 'rb')
species_logN_samples = pickle.load(pkl_file)
pkl_file.close()
species_logN_interp = {}

for s in list(species_logN_samples.keys()):
    species_logN_samples[s][species_logN_samples[s]==-np.inf] = -99
    species_logN_interp[s] = RegularGridInterpolator((logN_HI_arr, log_hdens_arr, log_metals_arr), 
                                                     species_logN_samples[s])

# Calculate $\log U$ from $n_\text{HI}$.
# fig, ax = plt.subplots(1, 1, figsize=(6, 4))
log_hdens_plot = np.arange(log_hdens_min, log_hdens_max+.01, .1)
logU_plot = np.log10(calc_U(hm12_wav_grid, hm12_J_nu_test, 10**log_hdens_plot))
# plt.plot(log_hdens_plot, logU_plot)
# plt.xlabel('log $n_{HI}$')
# plt.ylabel('log $U$')
# plt.show()

# Consider a particular $N_\text{HI} = 13$ and metallicity 0.
logN_HI_test = 16
log_m_test = 0
plot_points = [[logN_HI_test, np.round(log_hdens,2), log_m_test] for log_hdens in log_hdens_plot]
logN_HI_idx = np.where(logN_HI_arr==logN_HI_test)[0][0]
log_m_idx = np.where(log_metals_arr==log_m_test)[0][0]
fig, ax = plt.subplots(1,1, figsize=(10,7))
  
ax.text(x=0.02,y=0.9,s=r'$\mathrm{[X/H]}$ = ' + '{}'.format(log_m_test), fontsize=14, transform=ax.transAxes)
ax.text(x=0.02,y=0.85,s=r'$\log N$(HI) = ' + '{}'.format(logN_HI_test), fontsize=14, transform=ax.transAxes)
ax.axhline(logN_HI_test, linestyle='--', lw=1, color='k', label='HI')

# Plot column densities
species = ['CII', 'CIII', 'NII', 'NIII', "SiII", 'SiIII', 'SiIV']
for s in species:
    ax.plot(logU_plot, species_logN_interp[ion_species_dict[s]](plot_points), label=s)
    ax.plot(logU_arr, species_logN_samples[ion_species_dict[s]][logN_HI_idx,:,log_m_idx], marker='.', linestyle='', color='k')

ax.set_xlim(logU_plot[-1], logU_plot[0])
ax_copy = ax.twiny()
new_ticks = np.arange(-5,2,1, dtype='float')
ax_copy.set_xlim(ax.get_xlim())
ax_copy.set_xticks(np.log10(calc_U(hm12_wav_grid, hm12_J_nu_test, np.power(10, new_ticks))))
ax_copy.set_xticklabels(np.int_(new_ticks))

ax_label = fig.add_subplot(111, frameon=False)
ax_label.set_xticks([])
ax_label.set_yticks([])
ax_label.set_xlabel(r'$\log U$', labelpad=20)
ax_label.set_ylabel(r'$\log [N_{\mathrm{ion}} (\mathrm{cm}^{-2})]$', labelpad=35)

ax_copy_label = fig.add_subplot(111, frameon=False)
ax_copy_label.set_xticks([])
ax_copy_label.set_yticks([])
ax_copy_label.set_xlabel(r'$\log [n_\mathrm{H} (\mathrm{cm}^{-2})]$', labelpad=35)
ax_copy_label.xaxis.set_label_position('top')

ax.legend()
#plt.tight_layout()
plt.savefig(f'plots/{logN_HI_test}.png')
plt.show()
