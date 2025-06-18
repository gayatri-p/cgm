import numpy as np
from cloudy_fit_lib import *
import pickle

# change this line as required
rootdir = '/home/gayatri/Documents/cloudy_scripts/blr/grids/'

# edit all nHIs
logN_HI_min_main = 12
logN_HI_max_main = 12.25
logN_HI_step_main = 0.25
logN_HI_arr_main = np.array([15, 16]) #np.arange(logN_HI_min_main, logN_HI_max_main+logN_HI_step_main, logN_HI_step_main)
file_list_main = create_grid_file_list(logN_HI_arr_main)

# edit hden variation
log_hdens_min_main = -4.25
log_hdens_max_main = 1
log_hdens_step_main = 0.25
log_hdens_arr_main = np.arange(log_hdens_min_main, log_hdens_max_main+log_hdens_step_main, log_hdens_step_main)

# edit metallicity variation
log_metals_min_main = -0.25
log_metals_max_main = 1
log_metals_step_main = 0.25
log_metals_arr_main = np.arange(log_metals_min_main, log_metals_max_main+log_metals_step_main, log_metals_step_main)

species_logN_samples_final = {}
logT_grid_final = np.zeros((len(logN_HI_arr_main), len(log_hdens_arr_main), len(log_metals_arr_main)))

for i in range(len(logN_HI_arr_main)):
    
    logN_HI = logN_HI_arr_main[i]
    
    # The filename corresponding to the current stopping HI column density
    filename = file_list_main[i]
    
    # Get list of densities and metallicities for this stopping HI column density 
    log_hdens_grid, log_metals_grid = read_grd_file(rootdir, filename)

    # Get average (log) HI temperatures for all grid points
    log_temps_grid = read_avr_file(rootdir, filename)

    # Get column densities for all species
    species_names, log_col_dens_grid = read_col_file(rootdir, filename)
    
    # iterate over HI density
    for j in range(len(log_hdens_arr_main)):
        log_hdens = log_hdens_arr_main[j]
        
        # iterate over mettalicity
        for k in range(len(log_metals_arr_main)):
            log_metals = log_metals_arr_main[k]
            
            # Get grid index number for the current n_H and metallicity
            idx = np.intersect1d(np.where(log_hdens_grid==log_hdens)[0], np.where(log_metals_grid==log_metals)[0])[0]

            # Isolate the average temperature and column density for all species
            log_temp = log_temps_grid[idx]
            logT_grid_final[i,j,k] = log_temp
            log_col_dens = log_col_dens_grid[idx]
            
            # For each species
            for l in range(len(species_names)):
                s = species_names[l]

                if s not in species_logN_samples_final.keys():
                    species_logN_samples_final[s] = -99.*np.ones((len(logN_HI_arr_main), 
                                                            len(log_hdens_arr_main),
                                                            len(log_metals_arr_main)))
                
                # Check for converged logN(HI)
                if np.round(log_col_dens[0], 2) == logN_HI:
                    species_logN_samples_final[s][i,j,k] = log_col_dens[l]
    
    print(f'Successfully compiled for logN_HI = {logN_HI}.')

print('\nWriting final grids.')
output = open(rootdir+'final_grid.pkl', 'wb')
pickle.dump(species_logN_samples_final, output)
output.close()

np.savetxt(rootdir+'final_flat_logT.dat', logT_grid_final.flatten())

species_fracs_samples_final = {s:np.zeros((len(logN_HI_arr_main),len(log_hdens_arr_main),len(log_metals_arr_main))) for s in species_names_ions}

# H fractions
N_H_grid_final = 10**species_logN_samples_final['#column density H']+10**species_logN_samples_final['H+']
species_fracs_samples_final['#column density H'] = 10**species_logN_samples_final['#column density H']/N_H_grid_final
species_fracs_samples_final['H+'] = 10**species_logN_samples_final['H+']/N_H_grid_final

N_He_grid_final = 10**species_logN_samples_final['He']+10**species_logN_samples_final['He+']+10**species_logN_samples_final['He+2']
species_fracs_samples_final['He'] = 10**species_logN_samples_final['He']/N_He_grid_final
species_fracs_samples_final['He+'] = 10**species_logN_samples_final['He+']/N_He_grid_final
species_fracs_samples_final['He+2'] = 10**species_logN_samples_final['He+2']/N_He_grid_final

# Grid of metallicities
log_metals_grid_final = np.zeros((len(logN_HI_arr_main),len(log_hdens_arr_main),len(log_metals_arr_main)))

for i in range(len(logN_HI_arr_main)):
    for j in range(len(log_hdens_arr_main)):
        log_metals_grid_final[i,j,:] = log_metals_arr_main
        
# starting from Li
for i in range(5,len(species_names_ions)):
    s = species_names_ions[i]
    species_fracs_samples_final[s] = 10**species_logN_samples_final[s]/(N_H_grid_final*10**log_metals_grid_final)

species_logf_samples_final = {s:np.log10(species_fracs_samples_final[s]) for s in list(species_fracs_samples_final.keys())}

output = open(rootdir+'final_grid_logf.pkl', 'wb')
pickle.dump(species_logf_samples_final, output)
output.close()
print(f'Final grids generated for {log_metals_grid_final.size} configurations.')