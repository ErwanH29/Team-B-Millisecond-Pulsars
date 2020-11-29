from evol import neut_gal_evol
from plotters import plot_neut_pos, plot_mw_iso, plot_MC, plot_hist
import matplotlib.pyplot as plt
from coordinates import get_final_coord
from file_logistics import fileworker
from bias import get_xi2
import numpy as np

print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
Neutron star launcher interface. Here, neutron star sets can be generated and plotted. Also, older sets can be plotted. \n \
If you want to plot older sets instead of generating a new one, just load old simulation data using this interface. \n \
If a new set is generated, these files will be built automatically and placed in working_data. \n \
After generating, you can save the files. This moves the files into a default (or user defined) folder. \n \
Aside from plotting, one can also generate statistical info for all saved datasets. \n \
If one would like to change the path to the simulation folders, change "default_path=False" in fileworker inputs. \n \
------------------------------------------------------------\
------------------------------------------------------------')

fl = fileworker(default_path=True)
number_of_workers = 6 # enter an appropriate value for your system

build_neut_set = input('generate a new set of neutron stars (y|n):')

if build_neut_set == 'n':
    load_data = input('load copy of old simulation data (y|n):')
    
    if load_data == 'y':
        fl.move_data_to_working()
        
plot_data = input('plot data (y|n):')
calc_encounters_dist = input('calculate close encounters distance table (y|n):')

if build_neut_set == 'y':
    pos = input('old or current position (old|current):')
    save_data = input('save generated data (y|n):')   
    neut_line = neut_gal_evol(number_of_workers=number_of_workers, old_pos=pos)
      
if plot_data =='y':
    fix = input('fix axis (x|y|z):')
    neut_plot = input('plot neutron stars (y|n):')
    
    if fix == 'z':
        MW_plot = input('plot milky way isopotential lines ((x,y) only) (y|n):')
    else:
        MW_plot = 'n'
        
    MC_plot = input('plot Magellanic Cloud trajectories (y|n):')
    
    if neut_plot =='y':
        use_all = input('plot all neutron stars (y|n):')
        lines = plot_neut_pos()
        lines.get_dataframe_from_pkl(dat_str='working_data/neut_stars_positions.pkl')
        check = open('working_data/check.txt', 'r')
        lines.plot(use_all=use_all, check = check, fix=fix)
            
    if MW_plot == 'y':
        plot_mw_iso()

    if MC_plot == 'y':
        plot_MC(fix=fix)
        plt.legend(loc='lower right')
        
    if use_all=='y' :
        plt.title('Ejected MSPs from LMC')
        
    if use_all=='n': 
        plt.title('Ejected MSPs close encounters')
       
    if fix == 'z':
        plane = 'xy'
        plt.xlabel(r"$x$ coordinate (kpc)")
        plt.ylabel(r"$y$ coordinate (kpc)")
        plt.savefig("neutron_star_trajectories_xy", dpi=300)
        
    if fix == 'x':
        plane = 'yz'
        plt.xlabel(r"$y$ coordinate (kpc)")
        plt.ylabel(r"$z$ coordinate (kpc)")
        plt.savefig("neutron_star_trajectories_yz", dpi=300)
        
    if fix == 'y':
        plane = 'xz'
        plt.xlabel(r"$x$ coordinate (kpc)")
        plt.ylabel(r"$z$ coordinate (kpc)")
        plt.savefig("neutron_star_trajectories_xz", dpi=300)
    
    print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
a figure "neutron_stars_trajectories_{}.png" has been generated. \n \
------------------------------------------------------------\
------------------------------------------------------------'.format(plane))
    
    plt.show()

if calc_encounters_dist == 'y':
    
    check = open('working_data/check.txt', 'r')
    get_final_coord(check=check)
    
if build_neut_set == 'y':
    if save_data == 'y':
            fl.move_data_to_sim_results()
    
get_statistics = input('get close encounter statistics. note: only takes the simulations \
in the "sim_results" folder into account. (y|n):')

if get_statistics == 'y':
    n_sim = fl.get_enc_rate()
    
    get_hist = input('generate final position histogram for the set of simulations (y|n):')
    if get_hist == 'y':
        pos = open('statistics/position_distribution_{}sims.txt'.format(n_sim), 'r')   
        data = plot_hist(pos=pos, width=1, n_sim=n_sim)
        print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
a figure "final_pos_prob_{}sims.png" has been generated. \n \
 ------------------------------------------------------------\
------------------------------------------------------------ '.format(n_sim))
        min_value, unbiased = get_xi2(data, 100)
        print('\
 unbiased estimator(xi^2 reach minimum): {} \n\
 mean value: {} \n\
 median value: {} \n\
 xi^2: {}\n \
------------------------------------------------------------\
------------------------------------------------------------'.format(unbiased, np.mean(data), np.median(data), min_value))
        
clear_working = input('clear working_data directory. WARNING: make sure the data in \n \
"working_data" is a copy of simulated and stored data or make sure the data is saved! (y|n):')

if clear_working == 'y':
    fl.clear_working_data()
