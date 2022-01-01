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

def model_params(i):
    """
    Function shich defines the background parameters.
    
    Input:
    i:      The ith iteration of the lists w0vals and wavals which define w0 and wa 
            values of CPL function and the different background parameters.
    
    Output: The defined model which is needed for model_simulation function.
    """
    value = flag[i]
    if value == 1:
        parameters = {'EFTflag':1,
                      'PureEFTmodel':1,
                      'EFTwDE':2,'EFTw0': w0vals[i],'EFTwa': wavals[i],
                      'EFTwn':0,'EFTwat':0, 'EFTw2':0, 'EFTw3':0,
                      }
        parameters.update(common_options)   

        model = camb.set_params(lmax=2500, As=np.exp(3.040)/1e10,
                                ns=0.9647, H0=66.88,
                                ombh2=0.02212, omch2=0.1206,
                                mnu=0.12, tau=0.0522,
                                num_massive_neutrinos=1,
                                nnu=3.046, EFTCAMB_params = parameters
                                )

    if value == 0:
        model = camb.set_params(lmax=2500, As=np.exp(3.040)/1e10,
                         ns=0.9647, H0=66.88,
                         ombh2=0.02212, omch2=0.1206,
                         mnu=0.12, tau=0.0522,
                         num_massive_neutrinos=1,
                         nnu=3.046, EFTCAMB_params={'EFTFlag':0}
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

#########################################################################################
# General model parameters: Stability options, redshifts
#                           window function and probed effects
#########################################################################################

common_options = {
                 'EFT_ghost_math_stability': False,
                 'EFT_mass_math_stability' : False,
                 'EFT_ghost_stability'     : False,
                 'EFT_gradient_stability'  : False,
                 'EFT_mass_stability'      : False,
                 'EFT_additional_priors'   : True,
                 'feedback_level': 0,
                 }

redshifts = np.array([0.1, 0.25, 0.50, 0.75, 1.0, 2.0, 3.0])
redshift_labels = [r'$z = 0.10$', r'$z = 0.25$', r'$z = 0.50$', r'$z = 0.75$', r'$z = 1.00$', r'$z = 2.00$', r'$z = 3.00$']
rs_fig_name = ['_z=0.10', '_z=0.25', '_z=0.50', '_z=0.75', '_z=1.00', '_z=2.00', '_z=3.00']
redshift_width = 0.01
windows = [ GaussianSourceWindow(redshift = i, source_type='gw', sigma=redshift_width) for i in redshifts ]
windows = windows + [ GaussianSourceWindow(redshift = i, source_type='sn', sigma=redshift_width) for i in redshifts ]

sim_constraint = [[True, True, True, True, True, True], 
                  [True, False, False, False, False, False],
                  [False, False, False, False, True, False],
                  [False, True, False, False, False, True],]

colors = ['black', 'C3', 'C0', 'C4', 'C1', 'C5']
line = ['-', '--', '--', '--', '--', '--']
sim_effect = ['Total Signal', 'Doppler Effect', 'Lensing Effect', 'Sachs-Wolfe Effect']

#########################################################################################
# CPL Cosmology - Parameters on w0, wa are based on arXiv: 1807.06209
# The models below look at parameters from Table 6 of the same paper for Planck+SNe+BAO
# Model 1: w0 = -1.000, wa = 0
# Model 2: w0 = -0.877, wa = 0.04
# Model 3: w0 = -1.037, wa = -0.55
#########################################################################################

model_labels = [r'$\Lambda$CDM', 'CPL Model 1', 'CPL Model 2']

flag = [0, 1, 1, 1, 1, 1]
w0vals = [-1.000, -0.877, -1.037]
wavals = [0,  0.04, -0.55]

results = []
Results_DiffModels = [[], [], [], [], [], []]  #List which holds the result of total signal for all models
results_list = []

for i in range(3):
    print("Simulating background cosmology:", model_labels[i])
    vals = model_params(i)
    results_list.append(camb.get_results(vals))
    for j in range(4):
        tot_effect, ell = model_simulation(vals, j)
        Results_DiffModels[i].append(tot_effect)
        results.append(vals)
    print("--------------------------------------------")

################################################################
# Plotting
################################################################
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

    for j in range(3):
        ax1.plot(ell, Results_DiffModels[j][0][0,0,:,i,i], lw=1., color = colors[j], 
                 linestyle = line[j], label = model_labels[j])
        ax1.text(0.04, 0.04, str(sim_effect[0]) + "   Redshift: " + str(redshift_labels[i]), verticalalignment='bottom', 
                 horizontalalignment='left', transform=ax1.transAxes,
                 fontsize=main_fontsize, bbox=dict(facecolor='white',
                 boxstyle='square',edgecolor='white'))

    fig.legend( labels = model_labels, fontsize = 0.9*main_fontsize, frameon   = True,
                fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

    plt.savefig('Thesis_Images/CPL_Plots/Total_Signal_wCDM'+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
    plt.close('')

redshift_vals = [0, -1]
for i in redshift_vals:
    for j in range(4):
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

        for k in range(3): #Plotting individual effects
            ax1.plot(ell, Results_DiffModels[k][j][0,0,:,i,i], lw=1., color = colors[k], linestyle = line[k], label = model_labels[k])
            ax1.text(0.04, 0.04, str(sim_effect[j]+'   Redshift: '+str(redshift_labels[i])), verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                    fontsize=main_fontsize, bbox=dict(facecolor='white',boxstyle='square',edgecolor='white'))
                
        fig.legend( labels = model_labels, fontsize = 0.9*main_fontsize, frameon   = True,
                    fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                    borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

        plt.savefig('Thesis_Images/CPL_Plots/Indiv_Effects_Compare_wCDM_model_redshift_'+str(rs_fig_name[i])+'_effect_'+str(sim_effect[j])+'.pdf', bbox_inches='tight')
        plt.close()

for i in redshift_vals: #Plotting the total and individual effects of each model for defined redshifts
    for j in range(3): #Number of models plotting
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

        plt.savefig('Thesis_Images/CPL_Plots/Indiv_Effects_wCDM_'+str(model_labels[j])+'_redshift'+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
        plt.close()

print('Simulation done.')
exit()