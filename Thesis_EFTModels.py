import os
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import utilities.utilities as utilities
import utilities.color_utilities as col
from utilities.log_lin_scale import *

# Import CAMB:
here = './'
camb_path = os.path.realpath(os.path.join(os.getcwd(),here+'/../'))

sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower
from camb.sources import GaussianSourceWindow, SplinedSourceWindow

np.printoptions(precision=3, suppress=True)
camb.set_feedback_level(0)
high_precision = True
main_fontsize = 9

indiv_effects = False #Simulates the individual effects. Set as false unless you want full results (lengthy)

#########################################################################################
# General simulation function
#########################################################################################
def model_parameters(model_setup):
    """
    Function which defines certain constants to set up the different cosmologies
    Input:
    model_setup: The EFTCAMB parameters
    Output:
    model: The updated cosmology for which the simulation runs on"""
    model = camb.set_params(lmax=2500, As=np.exp(3.040)/1e10,
                                ns=0.9647, H0=66.88,
                                ombh2=0.02212, omch2=0.1206,
                                mnu=0.12, tau=0.0522,
                                num_massive_neutrinos=1,
                                nnu=3.046, EFTCAMB_params = model_setup
                                )
    return model

def model_simulation(model_sim, val_eff):
    """
    Bulk of the code. Simulates the luminosity distance angular power spectrum 
    for a range of parameters and models.
    
    Input: 
    model_sim: The cosmological model to be simulated.
    val_eff:   The effect (total/doppler/lensing/SW) wanting to be simulated. 
    Output:    Luminosity distance angular power spectrum for GW only
    """

    results = []
    print("...Processing " + sim_effect[val_eff] + "...")
    for setting in [model_sim]:
        setting.set_for_lmax(1000)
        setting.NonLinear = camb.model.NonLinear_none
        setting.SourceTerms.limber_phi_lmin = 20
        setting.Transfer.high_precision = True
        setting.SourceWindows = windows
        if high_precision:
            setting.SourceTerms.limber_windows = False
            setting.Accuracy.SourceLimberBoost = 1.0
            setting.Accuracy.LimberBoost = 1.0
            setting.Accuracy.TimeStepBoost = 1.1
            setting.Accuracy.IntkAccuracyBoost = 1.1
            setting.Accuracy.KmaxBoost = 0.6
        else:
            setting.SourceTerms.limber_windows = True
            setting.Accuracy.TimeStepBoost = 1
            setting.Accuracy.IntkAccuracyBoost = 0.5
            setting.Accuracy.KmaxBoost = 1

    for setting in [model_sim]:
        results_temp = []
        for ind in range(1):
            
            setting.SourceTerms.gw_velocity = sim_constraint[val_eff][0]
            setting.SourceTerms.gw_sw = sim_constraint[val_eff][1]
            setting.SourceTerms.gw_volume = sim_constraint[val_eff][2]
            setting.SourceTerms.gw_sh = sim_constraint[val_eff][3]
            setting.SourceTerms.gw_lensing = sim_constraint[val_eff][4]
            setting.SourceTerms.gw_isw = sim_constraint[val_eff][5]
            setting.SourceTerms.sn_velocity = sim_constraint[val_eff][0]
            setting.SourceTerms.sn_sw = sim_constraint[val_eff][1]
            setting.SourceTerms.sn_volume = sim_constraint[val_eff][2]
            setting.SourceTerms.sn_sh = sim_constraint[val_eff][3]
            setting.SourceTerms.sn_lensing = sim_constraint[val_eff][4]
            setting.SourceTerms.sn_isw = sim_constraint[val_eff][5]
            
            temp_results = camb.get_results(setting)
            ell, Cls = utilities.get_Cls(temp_results, raw_cl=False)
            results_temp.append(copy.deepcopy(Cls))
        results.append(copy.deepcopy(results_temp))
    results = np.array(results)

    # separate effects and build joint estimator:
    index_gw = [ ind for ind, sour in enumerate(setting.SourceWindows) if sour.source_type == 'gw' ]
    results_gw1 = results[:,:,:,index_gw,:][:,:,:,:,index_gw]
    return results_gw1, ell

def plotting_function(model_labels, model_results, fig_name, fig_fixed, fig_name_fixed):
    """
    Main function to process results and plot the graphs
    input:
    model_labels: The labels of the models simulated
    model_results: The different models you want to test holding the needed background parameters
    fig_name: Set of labels to produce cleaner graph outputs
    fig_fixed: A fixed string as caption of the plot
    returns: various plots of the models simulated.
    """
    results = []
    Results_DiffModels = [[], [], [], [], [], [], [], []]  #List which holds the result of total signal for all models
    results_list = []

    for i in range(len(model_labels)):
        print("Simulating cosmology:", model_labels[i])
        
        for j in range(4):
            results_list.append(camb.get_results(model_results[i]))
            tot_effect, ell = model_simulation(model_results[i], j)
            Results_DiffModels[i].append(tot_effect)
            results.append(model_results[i])
        
        print("--------------------------------------------")

    print("...Plotting...")

    for i in range(len(redshifts)): #Plot for each redshift the total signal depending on cosmological model
        fig = plt.gcf()

        x_size = 13
        y_size = 8.0 +0.5
        fig.set_size_inches(x_size/2.54, y_size/2.54 )
        gs = gridspec.GridSpec(1,1)
        ax1 = plt.subplot(gs[0,0])

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r'Monopole Number [$\ell$]', fontsize=main_fontsize, labelpad=4)
        ax1.set_ylabel(r'$\ell(\ell+1)C^{GW}_\ell/2\pi$', fontsize=main_fontsize)
        ax1.set_xlim([2, 1000])
        ax1.set_ylim([1e-8, 1e-3])

        ticks  = [ 2, 10, 100, 1000 ]
        ax1.set_xticks(ticks);
        ax1.set_xticklabels( ticks, horizontalalignment='center', 
                            fontsize=0.9*main_fontsize);
        ax1.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right');
        ax1.xaxis.get_majorticklabels()[0].set_horizontalalignment('left');
        ax1.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)

        for j in range(len(model_labels)):
            ax1.plot(ell, Results_DiffModels[j][0][0,0,:,i,i], lw=1., color = colors[j], 
                    linestyle = line[j], label = model_labels[j])
            ax1.text(0.04, 0.04, str(fig_fixed)+'   Redshift: '+str(redshift_labels[i]), verticalalignment='bottom', 
                    horizontalalignment='left', transform=ax1.transAxes,
                    fontsize=main_fontsize, bbox=dict(facecolor='white',
                    boxstyle='square',edgecolor='white'))
        fig.legend( labels = model_labels, fontsize = 0.9*main_fontsize, frameon   = True,
                    fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                    borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

        plt.savefig('Thesis_Images/EFT_Gen_Plots/Total_Signal/Model_'+str(fig_name_fixed)+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
        plt.close('')

    for i in redshift_vals: #Plotting the total and individual effects of each model for defined redshifts
        for j in range(len(model_labels)): #Number of models plotting
            fig = plt.gcf()

            x_size = 13
            y_size = 8.0 +0.5
            fig.set_size_inches(x_size/2.54, y_size/2.54 )
            gs = gridspec.GridSpec(1,1)
            ax1 = plt.subplot(gs[0,0])
                        
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlabel(r'Monopole Number [$\ell$]', fontsize=main_fontsize, labelpad=4)
            ax1.set_ylabel(r'$\ell(\ell+1)C^{GW}_\ell/2\pi$', fontsize=main_fontsize)
            ax1.set_xlim([2, 1000])
            ax1.set_ylim([1e-11, 1e-3])

            ticks  = [ 2, 10, 100, 1000 ]
            ax1.set_xticks(ticks);
            ax1.set_xticklabels( ticks, horizontalalignment='center', fontsize=0.9*main_fontsize);
            ax1.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right');
            ax1.xaxis.get_majorticklabels()[0].set_horizontalalignment('left');
            ax1.tick_params(axis='both', which='both', direction='in',
                                bottom=True, top=True, left=True, right=True)

            for k in range(4): #Plotting individual effects
                ax1.plot(ell, Results_DiffModels[j][k][0,0,:,i,i], lw=1., color = colors[k], linestyle = line[k], label = sim_effect[k])
                ax1.text(0.04, 0.04, "Model: "+str(model_labels[j])+"   Redshift: "+str(redshift_labels[i]), 
                        verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                        fontsize=main_fontsize, bbox=dict(facecolor='white',boxstyle='square',edgecolor='white'))
                
            fig.legend( labels = sim_effect, fontsize = 0.9*main_fontsize, frameon   = True,
                        fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                        borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

            plt.savefig('Thesis_Images/EFT_Gen_Plots/Indiv_Effects/Model_'+str(fig_name[j])+'_redshift'+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
            plt.close()

    return 

#########################################################################################
# General model parameters: Stability options, redshifts
#                           window function and probed effects
#########################################################################################

common_options = {
                 'EFT_ghost_math_stability': True,
                 'EFT_mass_math_stability' : True,
                 'EFT_ghost_stability'     : False,
                 'EFT_gradient_stability'  : False,
                 'EFT_mass_stability'      : False,
                 'EFT_additional_priors'   : False,
                 'feedback_level': 0,
                 }

redshifts = np.array([0.1, 0.2, 0.5, 1.0])
redshift_labels = [r'$z = 0.10$', r'$z = 0.20$', r'$z = 0.50$', r'$z = 1.00$']
redshift_vals = [0, -1]

redshift_width = 0.01
windows = [ GaussianSourceWindow(redshift = i, source_type='gw', sigma=redshift_width) for i in redshifts ]
windows = windows + [ GaussianSourceWindow(redshift = i, source_type='sn', sigma=redshift_width) for i in redshifts ]

sim_effect = ['Total Signal', 'Doppler Effect', 'Lensing Effect', 'Sachs-Wolfe Effect']
rs_fig_name = ['_z=0.10', '_z=0.20', '_z=0.50', '_z=1.0', '_z=2.0']
colors = ['black', 'C3', 'C0', 'C4', 'C1', 'C5']
line = ['-', '--', '--', '--', '--', '--']

sim_constraint = [[True, True, True, True, True, True], 
                  [True, False, False, False, False, False],
                  [False, False, False, False, True, False],
                  [False, True, False, False, False, True],]

GRmodel = model_parameters({'EFTFlag':0})
Results_DiffModels = [[], [], [], [], [], [], [], []] 
print(np.shape(Results_DiffModels), type(Results_DiffModels))


############################################################################################
# Omega0 linear function models
#############################################################################################

model_labels_lin = [r'$\Lambda$CDM', r'$\Omega_0 = 0.05$', r'$\Omega_0 = 0.10$', r'$\Omega_0 = 0.50$', r'$\Omega_0 = 1.00$'] 
fig_name_lin = ['LCDM', 'Omega0=0.05', 'Omega0=0.10', 'Omega0=0.50', 'Omega0=1.00']
fig_fixed_lin = '$\Omega(a) = \Omega_0 a$'
fig_name_fixed_lin = 'Omega_lin_'
Omega_param_lin = []
Omega_results_lin = [GRmodel, ]
Omega0_vals_lin = [0.05, 0.10, 0.30, 0.50, 1.00]

for i in range(len(Omega0_vals_lin)):
    print("Computing linear Omega model: ", model_labels_lin[i])
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0,
                  'PureEFTmodelOmega':2,'EFTOmega0': Omega0_vals_lin[i]
                  }
    Omega_param_lin.append(parameters)
    Omega_results_lin.append(model_parameters(parameters))

print("Computing linear models:")
plotting_function(model_labels_lin, Omega_results_lin, fig_name_lin, fig_fixed_lin, fig_name_fixed_lin)

#########################################################################################
# Model 1: Standard GR
# Model 2: Designer f(R) based on Alicia's code
#########################################################################################
model_labels_grl = [r'$\Lambda$CDM', r'$f(R)$, $B_0 = 0.01$',r'$f(R)$, $B_0 = 0.1$', r'$f(R)$, $B_0 = 1.0$', r'$f(R)$, $B_0 = 2.0$']
fig_fixed_grl = '$f(R)$ gravity'
fig_name_fixed_fR = 'fR_grav_'
fig_name_grl = ['LCDM', 'f(R)_B0.01', 'f(R)_B.1', 'f(R)_B1', 'f(R)_B2']
B0_parameters = [0.01, 0.1, 1, 2]

general_results = [GRmodel, ]

for i in range(len(B0_parameters)):
    parameters = {'EFTflag' :3,
                'DesignerEFTmodel' :1,
                'EFTwDE': 0,
                'EFTB0': B0_parameters[i],
                }
    parameters.update(common_options) 
    fRmodel = model_parameters(parameters)
    general_results.append(fRmodel)

print(general_results)
print(type(general_results))

print('Computing general models:')
plotting_function(model_labels_grl, general_results, fig_name_grl, fig_fixed_grl, fig_name_fixed_fR)

#############################################################################################
# Omega0 power law function models
#############################################################################################

model_labels_powerlaw = [r'$\gamma_1^0 = -10$', r'$\gamma_1^0 = -1.00$', r'$\gamma_1^0 = -0.10$'] 
fig_name_powerlaw = ['gamma=-10', 'gamma=-1', 'gamma=-0.1', 'gamma=-0.01']
fig_fixed_powerlaw = '$\Omega(a) = \Omega_0 a^{\gamma_1^0}$'
Omega_param_powerlaw = []
Omega_results_powerlaw = [ ]
fig_name_fixed_power = 'Power_law_'
Omega0_vals_powerlaw = [-10, -1, -0.1, -0.01]
Omega0_PL = [0.75]

for i in range(len(Omega0_vals_powerlaw)):
    print("Computing power law Omega model (Omega0 = 1.0): ",  Omega0_vals_powerlaw[i])    
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0,
                  'PureEFTmodelOmega':2,'EFTOmega0': 1,'EFTOmegaExp': Omega0_vals_powerlaw[i],
                  'PureEFTmodelGamma1': 1
                  }
    Omega_param_powerlaw.append(parameters)
    Omega_results_powerlaw.append(model_parameters(parameters))

print('Computing Omega power law models (Omega0 = 1.0):')
plotting_function(model_labels_powerlaw, Omega_results_powerlaw, fig_name_powerlaw, fig_fixed_powerlaw, fig_name_fixed_power)

#############################################################################################
# Omega0 exponential function models
#############################################################################################

model_labels_exponential = [r'$\gamma_1^0 = -0.01$', r'$\gamma_1^0 = -1.00$', r'$\gamma_1^0 = -100$'] 
fig_name_exponential = ['Omega1_exponential_neg', 'alpha=-0.001', 'alpha=-0.01', 'alpha=-0.1', 'alpha=-100']
fig_fixed_exponential = '$\Omega(a) = \Omega_0 a^{\gamma_1^0}$'
Omega_param_exponential = []
fig_name_fixed_exp = 'Exponential_'
Omega_results_exponential = []
Omega0_vals_exponential = [-0.001, -0.01, -1, -100]

for i in range(len(Omega0_vals_exponential)):
    print("Computing exponential law Omega model (Omega0 = 1.0): ", i)
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0,
                  'PureEFTmodelOmega':3,'EFTOmega0': 0.50,
                  'EFTOmegaExp': Omega0_vals_exponential[i]
                  }
    Omega_param_exponential.append(parameters)
    Omega_results_exponential.append(model_parameters(parameters))

print('Computing Omega exp. law models (Omega0 = 1.0):')
plotting_function(model_labels_exponential, Omega_results_exponential, fig_name_exponential, fig_fixed_exponential, fig_name_fixed_exp)

print('Simulation done.')
exit()