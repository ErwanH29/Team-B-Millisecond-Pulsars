from evol import neut_gal_evol
from plotters import plot_neut_pos, plot_mw_iso, plot_MC, plot_hist
import matplotlib.pyplot as plt
from coordinates import get_final_coord

print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
Neutron star launcher interface. Here, neutron star sets can be generated and plotted. Also, older sets can be plotted. \n \
If you want to plot older sets instead of generating a new one, make sure a "neut_stars_positions.pkl" \n \
and its corresponding "check.txt" and "gal_line.pickle" are present in the current working directory. \n \
If a new set is generated, these files will be built automatically. \n \
------------------------------------------------------------\
------------------------------------------------------------')

build_neut_set = input('generate a new set of neutron stars (y|n):')
plot_data = input('plot data (y|n):')
plot_encounters_hist = input('calculate close encounters distance table (y|n):')

if build_neut_set == 'y':
    pos = input('old or current position (old|current):')
    neut_gal_evol(old_pos=pos)

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
        lines.get_dataframe_from_pkl(dat_str='neut_stars_positions.pkl')
        check = open('check.txt', 'r')
        lines.plot(use_all=use_all, check = check, fix=fix)
            
    if MW_plot == 'y':
        plot_mw_iso()

    if MC_plot == 'y':
        plot_MC(fix=fix)
        
    if use_all=='y' :
        plt.title('Ejected MSPs from LMC')
        
    if use_all=='n': 
        plt.title('Ejected MSPs close encounters')
    
    plt.legend(loc='lower right')
    
    if fix == 'z':
        plt.xlabel(r"$x$ coordinate (kpc)")
        plt.ylabel(r"$y$ coordinate (kpc)")
        plt.savefig("neutron_star_trajectories_xy", dpi=300)
        
    if fix == 'x':
        plt.xlabel(r"$y$ coordinate (kpc)")
        plt.ylabel(r"$z$ coordinate (kpc)")
        plt.savefig("neutron_star_trajectories_yz", dpi=300)
        
    if fix == 'y':
        plt.xlabel(r"$x$ coordinate (kpc)")
        plt.ylabel(r"$z$ coordinate (kpc)")
        plt.savefig("neutron_star_trajectories_xz", dpi=300)
    
    plt.show()

if plot_encounters_hist == 'y':
    
    check = open('check.txt', 'r')
    get_final_coord(check=check)
    
    print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
distances.txt has been generated, containing the galactocentric distances of the neutron star set in kpc. \n \
------------------------------------------------------------\
------------------------------------------------------------')
