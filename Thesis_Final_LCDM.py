import copy
import os
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

import utilities.color_utilities as col
import utilities.utilities as utilities
from utilities.log_lin_scale import *


def model_parameters(param):
    """
    Function which defines the particular cosmology.
    input:
    param:  Parameters which define the cosmology.
    output: The updated cosmology for which the simulation runs on.
    """
    
    model = camb.set_params(lmax=2500, As=np.exp(3.040)/1e10,
                                ns=0.9647, H0=66.88,
                                ombh2=0.02212, omch2=0.1206,
                                mnu=0, tau=0.0522,
                                num_massive_neutrinos=0,
                                nnu=3.046, **param
                                )

    return model

def list_former(variables, figure_label, file_label, math_name, name):
    """
    Function which initialises the models and provides labels.
    input:
    variables:    An array of different parameters which define the model.
    figure_label: An empty array to store strings uses in figure output.
    file_label:   An empty array to store strings used in plot legends.
    math_name:    Fixed string of the general gravitational theory in LaTex form.
    name:         Fixed string of the general gravitational theory.
    output:
    figure_label: Array of strings with the various labels defining the model used for figures.
    file_label:   Array of strings with the various labels defining the model used for output files.
    """

    for i in variables:
        name1 = str(math_name)+str.format('{:.3f}',(i))
        figure_label.append(str(name1))
        name2 = str(name)+str.format('{:.3f}',(i))
        file_label.append(str(name2))
    
    return figure_label, file_label
    

def file_write(model_initial, parameters, file_label, model_fixed, background, it):
    """
    A function which stores information of model parameters in a text file.
    input:
    model_initial: An array holding information on the model parameters.
    parameters:    Array holding information on defined parameters of the model.
    file_label:    A string to provide cleaner output.
    model_fixed:   The gravitational group being looked at.
    background:    Simulated background cosmology (LCDM, wCDM or CPL).
    it:            Fixed number to prevent data from overwriting.
    output:        A text file with stored information on parameters looked at
    """
    
    num_models = len(model_initial)
    textfile = open("Thesis_Data/Final_Data/Parameters_"+str(background)+"_"+str(model_fixed)+"_"+str(it)+".txt", "w")

    for i in range(num_models):
        if i == 0:
            textfile.write("-------------------------------------------- \n")
            textfile.write(str(file_label[i]) + "\n")
            textfile.write("-------------------------------------------- \n")
            textfile.write(str(model_initial[i]) + "\n")
            textfile.write(str((parameters))+"\n\n\n")
        else:
            textfile.write("-------------------------------------------- \n")
            textfile.write(str(file_label[i]) + "\n")
            textfile.write("-------------------------------------------- \n")
            textfile.write(str((parameters[i-1]))+"\n\n\n")
    textfile.close()

    return

def model_simulation(model_initial, effect_ind):
    """
    Function which processes results from a given cosmology
    input:
    model_initial: Array consisting of the initial model parameters.
    effect_ind:    Index which defines the effect wanting to simulate.
    output:
    results_gw:    Array holding results of the GW luminosity distance.
    results_delta: Array holding results of dark energy clustering.
    ell:           Maximum value of ell, used for x-axis in figures.
    """

    results = []
    print("...Processing " + sim_effect[effect_ind] + "...", flush = True)
    for setting in [model_initial]:
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

    for setting in [model_initial]:
        results_temp = []
        for ind in range(1):
            
            setting.SourceTerms.gw_velocity = sim_constraint[effect_ind][0]
            setting.SourceTerms.gw_sw = sim_constraint[effect_ind][1]
            setting.SourceTerms.gw_volume = sim_constraint[effect_ind][2]
            setting.SourceTerms.gw_sh = sim_constraint[effect_ind][3]
            setting.SourceTerms.gw_lensing = sim_constraint[effect_ind][4]
            setting.SourceTerms.gw_isw = sim_constraint[effect_ind][5]
            setting.SourceTerms.sn_velocity = sim_constraint[effect_ind][0]
            setting.SourceTerms.sn_sw = sim_constraint[effect_ind][1]
            setting.SourceTerms.sn_volume = sim_constraint[effect_ind][2]
            setting.SourceTerms.sn_sh = sim_constraint[effect_ind][3]
            setting.SourceTerms.sn_lensing = sim_constraint[effect_ind][4]
            setting.SourceTerms.sn_isw = sim_constraint[effect_ind][5]
            
            temp_results = camb.get_results(setting)
            ell, Cls = utilities.get_Cls(temp_results, raw_cl=False)
            results_temp.append(copy.deepcopy(Cls))
        results.append(copy.deepcopy(results_temp))
    results = np.array(results)

    index_sn = [ ind for ind, sour in enumerate(setting.SourceWindows) if sour.source_type == 'sn' ]
    index_gw = [ ind for ind, sour in enumerate(setting.SourceWindows) if sour.source_type == 'gw' ]

    results_sn = results[:,:,:,index_sn,:][:,:,:,:,index_sn]
    results_gw = results[:,:,:,index_gw,:][:,:,:,:,index_gw]
    results_sngw = results[:,:,:,index_sn,:][:,:,:,:,index_gw]
    results_delta = results_sn +results_gw - 2.*results_sngw

    return results_gw, results_delta, ell

def plotting_function(results_gw, results_DE, figure_label, figure_fixed, fixed, iteration, bckg, ell, text_x1, text_y1, text_x2, text_y2, text_x3, text_y3):
    """
    Function which processes all the plots
    
    input:
    results_gw:    Power spectrum data for GW luminosity distance.
    results_DE:    Power spectrum data for dark energy clustering.
    figure_label:  Array of strings holding labels for final plots in LaTex form.
    figure_fixed:  String denoting the general model of gravitational theory plotted in LaTex form.
    fixed:         String denoting the general model of gravitational theory plotted.
    file_label:    String denoting individual models for cleaner plot outputs.
    iteration:     The iteration number ranging from 1-3 (simulations of theories are grouped as 5 with 15 different parameters)
    ylim:          Value to scale the y-axis for GW-SN plots.
    bckg:          Background cosmology for plot saving.
    ell:           The x-axis for the plots.
    text_xB:       x position to place the fixed text along the Bth independent plot.
    text_yB:       y position to place the fixed text along the Bth independent plot.
    
    output:        Saved plots, all stored in /Thesis_Images/
    """
    
    print("...Plotting...")

    for i in range(len(redshifts)): #Plot total signal (GW) of all models for fixed redshifts
        fig = plt.gcf()

        x_size = 30
        y_size = 8.0 +0.5
        fig.set_size_inches( x_size/2.7, y_size/2.7 )
        gs = gridspec.GridSpec(1,2)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])

        ax1.set_ylim([1e-11, 1e-3])
        ax2.set_ylim([1e-18, 1e-3])
        ax1.set_ylabel(r'$\ell(\ell+1)C^{GW}_\ell/2\pi$', fontsize=main_fontsize)

        for _ax in [ax1, ax2]:        
            _ax.set_xscale('log')
            _ax.set_yscale('log')
            _ax.set_xlabel(r'Monopole Number [$\ell$]', fontsize=main_fontsize, labelpad=4)
            _ax.set_xlim([2, 1000])

            ticks  = [ 2, 10, 100, 1000 ]
            _ax.set_xticks(ticks);
            _ax.set_xticklabels( ticks, horizontalalignment='center', fontsize=0.9*main_fontsize);
            _ax.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right');
            _ax.xaxis.get_majorticklabels()[0].set_horizontalalignment('left');
            _ax.tick_params(axis='both', which='both', direction='in',
                                    bottom=True, top=True, left=True, right=True)

        for j in range(len(figure_label)):
            ax1.plot(ell, results_gw[j][0][0,0,:,i,i], lw=1., color = colors[j], linestyle = line[j], label = figure_label[j])
            ax2.plot(ell, results_gw[j][0][0,0,:,i,i], lw=1., color = colors[j], linestyle = line[j])
            ax2.plot(ell, results_gw[j][1][0,0,:,i,i], lw=1., color = colors[j], linestyle = ":")                       

        fig.legend( labels = figure_label, fontsize = 0.75*main_fontsize, frameon   = True,
                    fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',  bbox_to_anchor=(0.5, 1.02),
                    borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

        ax1.text(text_x1, text_y1, "Total Signal:   "+str(figure_fixed)+'\nRedshift:         '+str(redshift_labels[i]), 
                 verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                 fontsize=main_fontsize, bbox=dict(facecolor='white', boxstyle='square',edgecolor='white')) 

        ax2.text(text_x2, text_y2, 'Dark Energy Clustering', verticalalignment='bottom',  horizontalalignment='left', 
                 transform=ax2.transAxes, fontsize=main_fontsize, bbox=dict(facecolor='white',
                 boxstyle='square',edgecolor='white')) 
                 
        plt.savefig('Thesis_Images/Final_Plots/Total_Signal/Total_Signal_'+str(fixed)+"_Plot_GW_"+str(bckg)+'.'+str(iteration)+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
        plt.close('')

    for i in range(len(redshifts)): #Plot total signal (GW-SN) of all models for fixed redshifts
        fig = plt.gcf()

        x_size = 7.5
        y_size = 8.0 +0.5
        fig.set_size_inches( x_size/1.4, y_size/2.7 )
        gs = gridspec.GridSpec(1,1)
        ax1 = plt.subplot(gs[0,0])
        ax1.set_ylim([1e-18, 1e-4])
        ax1.set_ylabel(r'$\ell(\ell+1)C^{\Delta\varphi}_\ell/2\pi$', fontsize=main_fontsize)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r'Monopole Number [$\ell$]', fontsize=main_fontsize, labelpad=4)
        ax1.set_xlim([2, 1000])

        ticks  = [ 2, 10, 100, 1000 ]
        ax1.set_xticks(ticks);
        ax1.set_xticklabels( ticks, horizontalalignment='center', fontsize=0.9*main_fontsize);
        ax1.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right');
        ax1.xaxis.get_majorticklabels()[0].set_horizontalalignment('left');
        ax1.tick_params(axis='both', which='both', direction='in',
                                    bottom=True, top=True, left=True, right=True)

        for j in range(len(figure_label)):
            index_maximum = np.where(np.max(results_DE[j][0][0,0,:,i,i]) == results_DE[j][0][0,0,:,i,i])
            ax1.plot(ell, results_DE[j][0][0,0,:,i,i], lw=1., color = colors[j], linestyle = '-', label = figure_label[j])
            ax1.scatter(ell[index_maximum], results_DE[j][0][0,0,index_maximum,i,i], color = colors[j], s = 1.5)

        fig.legend( labels = figure_label, fontsize = 0.75*main_fontsize, frameon   = True,
                    fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',  bbox_to_anchor=(0.5, 1.02),
                    borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)
        
        for j in range(len(figure_label)):
            ax1.plot(ell, results_DE[j][1][0,0,:,i,i], lw=1., color = colors[j], linestyle = ":")                       

        ax1.text(text_x3, text_y3, "Total Signal:   "+str(figure_fixed)+'\n      Redshift:   '+str(redshift_labels[i]), 
                 verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                 fontsize=main_fontsize, bbox=dict(facecolor='white', boxstyle='square',edgecolor='white')) 
                 
        plt.savefig('Thesis_Images/Final_Plots/Total_Signal/Total_Signal_'+str(fixed)+"_Plot_GW-SN_"+str(bckg)+'.'+str(iteration)+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
        plt.close('')

    return 

# Import CAMB:
here = './'
camb_path = os.path.realpath(os.path.join(os.getcwd(),here+'/../'))

sys.path.insert(0,camb_path)
import camb
from camb import initialpower, model, eftcamb
from camb.sources import GaussianSourceWindow, SplinedSourceWindow

np.printoptions(precision=3, suppress=True)
camb.set_feedback_level(0)
high_precision = True
main_fontsize = 9

#########################################################################################
# Global parameters
#########################################################################################

stability_flag = {
                  'EFT_ghost_math_stability'   : True,
                  'EFT_mass_math_stability'    : False,
                  'EFT_ghost_stability'        : False,
                  'EFT_gradient_stability'     : False,
                  'EFT_mass_stability'         : False,
                  'EFT_additional_priors'      : False,
                  }

sim_constraint = [[True, True, True, True, True, True], 
                  [False,False,False,False,False,False,True],
                  [True, False, False, False, False, False],
                  [False, False, False, False, True, False],
                  [False, True, False, False, False, True],]

sim_effect = ['Total Signal', 'Dark Energy Clustering', 'Doppler Effect', 'Lensing Effect', 'Sachs-Wolfe Effect']

redshifts = np.array([0.05, 0.1, 0.5, 1.0])
redshift_vals = [0, -1]
redshift_labels = [r'$z = 0.05$', r'$z = 0.10$', r'$z = 0.50$', r'$z = 1.00$']
rs_fig_name = ['_z=0.05', '_z=0.10', '_z=0.50', '_z=1.0']

colors = ['black', 'C3', 'C0', 'C4', 'C1', 'C2']
line = ['-', '--', '--', '--', '--', '--']
line2 = ['-', ':', ':', ':', ':', ':']

redshift_width = 0.01
windows = [ GaussianSourceWindow(redshift = i, source_type='gw', sigma=redshift_width) for i in redshifts ]
windows = windows + [ GaussianSourceWindow(redshift = i, source_type='sn', sigma=redshift_width) for i in redshifts ]

Standard_GR = model_parameters({'EFTflag': 0})
LCDM_bckg = 'LCDM'
iteration = 0

#########################################################################################
####### Horava Gravity ##################################################################
#########################################################################################

Horava_test_fixed = 'Horava_test_lambda_grav'
Horava_test_figure_fixed = 'Horava_test Gravity'
Horava_test_figure_label = ['GR', ]
Horava_test_file_label = ['General Relativity', ]

Horava_test_math_name = r'$\lambda = $'
Horava_test_name = 'Horava_test_lambda_'

Horava_test_models = [Standard_GR, ]
Horava_test_lambda = [1e-3, 1e-2, 1e-1]
Horava_test_results_gw = np.empty((len(Horava_test_lambda)+1, 0)).tolist()
Horava_test_results_DE = np.empty((len(Horava_test_lambda)+1, 0)).tolist()
Horava_test_parameters_string = [ ]
Horava_test_figure_label, Horava_test_file_label = list_former(Horava_test_lambda, Horava_test_figure_label, Horava_test_file_label, Horava_test_math_name, Horava_test_name)
    
for i in range(len(Horava_test_lambda)):
    parameters = {'EFTflag':4,
                    'EFTflag':4, 'dark_energy_model' : 'EFTCAMB',
                    'FullMappingEFTmodel':2,
                    }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Horava_test_models.append(model_temp)
    Horava_test_parameters_string.append(parameters)
    
print("Iteration: ", iteration)
file_write(Horava_test_models, Horava_test_parameters_string, Horava_test_file_label, Horava_test_fixed, LCDM_bckg, iteration)

for i in range(len(Horava_test_models)):
    print("Simulating cosmology: ", Horava_test_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Horava_test_models[i], j)
        Horava_test_results_gw[i].append(final_res_gw)
        Horava_test_results_DE[i].append(final_res_DE)

plotting_function(Horava_test_results_gw, Horava_test_results_DE, Horava_test_figure_label, Horava_test_figure_fixed, Horava_test_fixed, 
                  iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)
#

Horava_fixed = 'Horava_lambda_grav'
Horava_figure_fixed = 'Horava Gravity'
Horava_figure_label = ['GR', ]
Horava_file_label = ['General Relativity', ]

Horava_math_name = r'$\lambda = $'
Horava_name = 'Horava_lambda_'

Horava_models = [Standard_GR, ]
Horava_lambda = [1.00001, 1.0000001, 1.000000001]
Horava_results_gw = np.empty((len(Horava_lambda)+1, 0)).tolist()
Horava_results_DE = np.empty((len(Horava_lambda)+1, 0)).tolist()
Horava_parameters_string = [ ]
Horava_figure_label, Horava_file_label = list_former(Horava_lambda, Horava_figure_label, Horava_file_label, Horava_math_name, Horava_name)
    
for i in range(len(Horava_lambda)):
    parameters = {'EFTflag':4,
                    'FullMappingEFTmodel':1,
                    'Horava_eta'   :    10**-4.5, 'dark_energy_model' : 'EFTCAMB',
                    'Horava_lambda':    Horava_lambda[i],
                    'Horava_xi'    :   0.99,
                    }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Horava_models.append(model_temp)
    Horava_parameters_string.append(parameters)
    
print("Iteration: ", iteration)
file_write(Horava_models, Horava_parameters_string, Horava_file_label, Horava_fixed, LCDM_bckg, iteration)

for i in range(len(Horava_models)):
    print("Simulating cosmology: ", Horava_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Horava_models[i], j)
        Horava_results_gw[i].append(final_res_gw)
        Horava_results_DE[i].append(final_res_DE)

plotting_function(Horava_results_gw, Horava_results_DE, Horava_figure_label, Horava_figure_fixed, Horava_fixed, 
                  iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

# Eta
Horava_fixed = 'Horava_eta_grav'
Horava_figure_fixed = 'Horava Gravity'
Horava_figure_label = ['GR', ]
Horava_file_label = ['General Relativity', ]

Horava_math_name = r'$\eta = $'
Horava_name = 'Horava_eta_'

Horava_models = [Standard_GR, ]
Horava_eta = [1e-11]
Horava_results_gw = np.empty((len(Horava_eta)+1, 0)).tolist()
Horava_results_DE = np.empty((len(Horava_eta)+1, 0)).tolist()
Horava_parameters_string = [ ]
Horava_figure_label, Horava_file_label = list_former(Horava_eta, Horava_figure_label, Horava_file_label, Horava_math_name, Horava_name)
    
for i in range(len(Horava_eta)):
    parameters = {'EFTflag':4,
                    'FullMappingEFTmodel':1, 'dark_energy_model' : 'EFTCAMB',
                    'Horava_eta'   :    Horava_eta[i],
                    'Horava_lambda':    1.00001,
                    'Horava_xi'    :    0.99,
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Horava_models.append(model_temp)
    Horava_parameters_string.append(parameters)
    
print("Iteration: ", iteration)
file_write(Horava_models, Horava_parameters_string, Horava_file_label, Horava_fixed, LCDM_bckg, iteration)

for i in range(len(Horava_models)):
    print("Simulating cosmology: ", Horava_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Horava_models[i], j)
        Horava_results_gw[i].append(final_res_gw)
        Horava_results_DE[i].append(final_res_DE)

plotting_function(Horava_results_gw, Horava_results_DE, Horava_figure_label, Horava_figure_fixed, Horava_fixed, 
                      iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

#########################################################################################
####### Linear Omega ####################################################################
#########################################################################################

Omega_Lin_fixed = 'Omega_Lin'
Omega_Lin_figure_fixed = '$\Omega(a) = \Omega_0 a$'
Omega_Lin_figure_label = ['GR', ]
Omega_Lin_file_label = ['General Relativity', ]

Omega_Lin_math_name = '$\Omega_0 = $'
Omega_Lin_name = 'Omega0_Lin_'

Omega_Lin_models = [Standard_GR, ]
Omega_Lin_parameters = [0.005, 0.025, 0.05, 0.075]
Omega_Lin_results_gw = np.empty((len(Omega_Lin_parameters)+1, 0)).tolist()
Omega_Lin_results_DE = np.empty((len(Omega_Lin_parameters)+1, 0)).tolist()
Omega_Lin_parameters_string = [ ]
Omega_Lin_figure_label, Omega_Lin_file_label = list_former(Omega_Lin_parameters, Omega_Lin_figure_label, 
                                                               Omega_Lin_file_label, Omega_Lin_math_name, Omega_Lin_name)

for i in range(len(Omega_Lin_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': Omega_Lin_parameters[i],
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_Lin_models.append(model_temp)
    Omega_Lin_parameters_string.append(parameters)

file_write(Omega_Lin_models, Omega_Lin_parameters_string, Omega_Lin_file_label, Omega_Lin_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_Lin_models)):
    print("Simulating cosmology: ", Omega_Lin_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_Lin_models[i], j)
        Omega_Lin_results_gw[i].append(final_res_gw)
        Omega_Lin_results_DE[i].append(final_res_DE)

plotting_function(Omega_Lin_results_gw, Omega_Lin_results_DE, Omega_Lin_figure_label, Omega_Lin_figure_fixed, 
                  Omega_Lin_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)
                  
#########################################################################################
####### Linear Omega (Negative) #########################################################
#########################################################################################

Omega_Neg_Lin_fixed = 'Omega_Neg_Lin'
Omega_Neg_Lin_figure_fixed = '$\Omega(a) = \Omega_0 a$'
Omega_Neg_Lin_figure_label = ['GR', ]
Omega_Neg_Lin_file_label = ['General Relativity', ]

Omega_Neg_Lin_math_name = '$\Omega_0 = $'
Omega_Neg_Lin_name = 'Omega0_Neg_Lin_'

Omega_Neg_Lin_models = [Standard_GR, ]
Omega_Neg_Lin_parameters = [-0.012, -0.006, 0.006, 0.012]
Omega_Neg_Lin_results_gw = np.empty((len(Omega_Neg_Lin_parameters)+1, 0)).tolist()
Omega_Neg_Lin_results_DE = np.empty((len(Omega_Neg_Lin_parameters)+1, 0)).tolist()
Omega_Neg_Lin_parameters_string = [ ]
Omega_Neg_Lin_figure_label, Omega_Neg_Lin_file_label = list_former(Omega_Neg_Lin_parameters, Omega_Neg_Lin_figure_label, 
                                                                       Omega_Neg_Lin_file_label, Omega_Neg_Lin_math_name, Omega_Neg_Lin_name)

for i in range(len(Omega_Neg_Lin_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': Omega_Neg_Lin_parameters[i],
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_Neg_Lin_models.append(model_temp)
    Omega_Neg_Lin_parameters_string.append(parameters)

file_write(Omega_Neg_Lin_models, Omega_Neg_Lin_parameters_string, Omega_Neg_Lin_file_label, Omega_Neg_Lin_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_Neg_Lin_models)):
    print("Simulating cosmology: ", Omega_Neg_Lin_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_Neg_Lin_models[i], j)
        Omega_Neg_Lin_results_gw[i].append(final_res_gw)
        Omega_Neg_Lin_results_DE[i].append(final_res_DE)

plotting_function(Omega_Neg_Lin_results_gw, Omega_Neg_Lin_results_DE, Omega_Neg_Lin_figure_label, Omega_Neg_Lin_figure_fixed, 
                  Omega_Neg_Lin_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)

#########################################################################################
# k-Mouflage (Fig. 6 arXiv:1809.09958)
#########################################################################################

##### Gamma_A
kMOuflage_max_gammaA_fixed = 'kMOuflage_max_gammaA1'
kMOuflage_max_gammaA_figure_fixed = 'kMOuflage_max_gammaA'
kMOuflage_max_gammaA_figure_label = ['GR', ]
kMOuflage_max_gammaA_file_label = ['General Relativity', ]

kMOuflage_max_gammaA_math_name = r'$\gamma_A = $'
kMOuflage_max_gammaA_name = 'kMOuflage_max_gammaA_'

kMOuflage_max_gammaA_models = [Standard_GR, ]
kMOuflage_max_gammaA_parameters = [0.1, 1, 10]
kMOuflage_max_gammaA_results_gw = np.empty((len(kMOuflage_max_gammaA_parameters)+1, 0)).tolist()
kMOuflage_max_gammaA_results_DE = np.empty((len(kMOuflage_max_gammaA_parameters)+1, 0)).tolist()
kMOuflage_max_gammaA_parameters_string = [ ]
kMOuflage_max_gammaA_figure_label, kMOuflage_max_gammaA_file_label = list_former(kMOuflage_max_gammaA_parameters, kMOuflage_max_gammaA_figure_label, 
                                                                       kMOuflage_max_gammaA_file_label, kMOuflage_max_gammaA_math_name, kMOuflage_max_gammaA_name)

for i in range(len(kMOuflage_max_gammaA_parameters)):
    parameters = {'EFTflag': 4, 'FullMappingEFTmodel': 3, 'dark_energy_model' : 'EFTCAMB',
                  'Kmimic': False, 'alphaU': 0.6, 'gammaU': 10.0, 
                  'm': 1.5, 'eps2_0': -0.037, 'gammaA': kMOuflage_max_gammaA_parameters[i],
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    kMOuflage_max_gammaA_models.append(model_temp)
    kMOuflage_max_gammaA_parameters_string.append(parameters)

file_write(kMOuflage_max_gammaA_models, kMOuflage_max_gammaA_parameters_string, kMOuflage_max_gammaA_file_label, kMOuflage_max_gammaA_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(kMOuflage_max_gammaA_models)):
    print("Simulating cosmology: ", kMOuflage_max_gammaA_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(kMOuflage_max_gammaA_models[i], j)
        kMOuflage_max_gammaA_results_gw[i].append(final_res_gw)
        kMOuflage_max_gammaA_results_DE[i].append(final_res_DE)

plotting_function(kMOuflage_max_gammaA_results_gw, kMOuflage_max_gammaA_results_DE, kMOuflage_max_gammaA_figure_label, kMOuflage_max_gammaA_figure_fixed, kMOuflage_max_gammaA_fixed, 
                  iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

##### Epsilon_2,0
kMOuflage_max_eps_fixed = 'kMOuflage_max_eps1'
kMOuflage_max_eps_figure_fixed = 'kMOuflage_max_eps'
kMOuflage_max_eps_figure_label = ['GR', ]
kMOuflage_max_eps_file_label = ['General Relativity', ]

kMOuflage_max_eps_math_name = '$\epsilon_{2,0} = $'
kMOuflage_max_eps_name = 'kMOuflage_max_eps_'

kMOuflage_max_eps_models = [Standard_GR, ]
kMOuflage_max_eps_parameters = [-0.015, -0.025, -0.035]
kMOuflage_max_eps_results_gw = np.empty((len(kMOuflage_max_eps_parameters)+1, 0)).tolist()
kMOuflage_max_eps_results_DE = np.empty((len(kMOuflage_max_eps_parameters)+1, 0)).tolist()
kMOuflage_max_eps_parameters_string = [ ]
kMOuflage_max_eps_figure_label, kMOuflage_max_eps_file_label = list_former(kMOuflage_max_eps_parameters, kMOuflage_max_eps_figure_label, 
                                                                       kMOuflage_max_eps_file_label, kMOuflage_max_eps_math_name, kMOuflage_max_eps_name)

for i in range(len(kMOuflage_max_eps_parameters)):
    parameters = {'EFTflag': 4, 'FullMappingEFTmodel': 3, 'dark_energy_model' : 'EFTCAMB',
                  'Kmimic': False, 'alphaU': 1, 'gammaU': 10.0, 
                  'm': 1.5, 'eps2_0': kMOuflage_max_eps_parameters[i], 'gammaA': 0.1,
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    kMOuflage_max_eps_models.append(model_temp)
    kMOuflage_max_eps_parameters_string.append(parameters)

file_write(kMOuflage_max_eps_models, kMOuflage_max_eps_parameters_string, kMOuflage_max_eps_file_label, kMOuflage_max_eps_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(kMOuflage_max_eps_models)):
    print("Simulating cosmology: ", kMOuflage_max_eps_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(kMOuflage_max_eps_models[i], j)
        kMOuflage_max_eps_results_gw[i].append(final_res_gw)
        kMOuflage_max_eps_results_DE[i].append(final_res_DE)

plotting_function(kMOuflage_max_eps_results_gw, kMOuflage_max_eps_results_DE, kMOuflage_max_eps_figure_label, kMOuflage_max_eps_figure_fixed, kMOuflage_max_eps_fixed, 
                  iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)


#########################################################################################
# k-Mouflage (Table 2 arXiv:1809.09958)
#########################################################################################
##### Epsilon_2,0
kMouflage_eps_fixed = 'kMouflage_eps1'
kMouflage_eps_figure_fixed = 'kMouflage_eps'
kMouflage_eps_figure_label = ['GR', ]
kMouflage_eps_file_label = ['General Relativity', ]

kMouflage_eps_math_name = '$\epsilon_{2,0} = $'
kMouflage_eps_name = 'kMouflage_eps_'

kMouflage_eps_models = [Standard_GR, ]
kMouflage_eps_parameters = [-1, -0.1, -0.01]
kMouflage_eps_results_gw = np.empty((len(kMouflage_eps_parameters)+1, 0)).tolist()
kMouflage_eps_results_DE = np.empty((len(kMouflage_eps_parameters)+1, 0)).tolist()
kMouflage_eps_parameters_string = [ ]
kMouflage_eps_figure_label, kMouflage_eps_file_label = list_former(kMouflage_eps_parameters, kMouflage_eps_figure_label, 
                                                                       kMouflage_eps_file_label, kMouflage_eps_math_name, kMouflage_eps_name)

for i in range(len(kMouflage_eps_parameters)):
    parameters = {'EFTflag': 4, 'FullMappingEFTmodel': 3, 'dark_energy_model' : 'EFTCAMB',
                  'Kmimic': False, 'alphaU': 0.1, 'gammaU': 2.0, 
                  'm': 5.0, 'eps2_0': kMouflage_eps_parameters[i], 'gammaA': 2,
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    kMouflage_eps_models.append(model_temp)
    kMouflage_eps_parameters_string.append(parameters)

file_write(kMouflage_eps_models, kMouflage_eps_parameters_string, kMouflage_eps_file_label, kMouflage_eps_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(kMouflage_eps_models)):
    print("Simulating cosmology: ", kMouflage_eps_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(kMouflage_eps_models[i], j)
        kMouflage_eps_results_gw[i].append(final_res_gw)
        kMouflage_eps_results_DE[i].append(final_res_DE)

plotting_function(kMouflage_eps_results_gw, kMouflage_eps_results_DE, kMouflage_eps_figure_label, kMouflage_eps_figure_fixed, kMouflage_eps_fixed, 
                  iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

##### Gamma_A
kMouflage_gammaA_fixed = 'kMouflage_gammaA1'
kMouflage_gammaA_figure_fixed = 'kMouflage_gammaA'
kMouflage_gammaA_figure_label = ['GR', ]
kMouflage_gammaA_file_label = ['General Relativity', ]

kMouflage_gammaA_math_name = r'$\gamma_A = $'
kMouflage_gammaA_name = 'kMouflage_gammaA_'

kMouflage_gammaA_models = [Standard_GR, ]
kMouflage_gammaA_parameters = [0.2, 2, 20]
kMouflage_gammaA_results_gw = np.empty((len(kMouflage_gammaA_parameters)+1, 0)).tolist()
kMouflage_gammaA_results_DE = np.empty((len(kMouflage_gammaA_parameters)+1, 0)).tolist()
kMouflage_gammaA_parameters_string = [ ]
kMouflage_gammaA_figure_label, kMouflage_gammaA_file_label = list_former(kMouflage_gammaA_parameters, kMouflage_gammaA_figure_label, 
                                                                       kMouflage_gammaA_file_label, kMouflage_gammaA_math_name, kMouflage_gammaA_name)

for i in range(len(kMouflage_gammaA_parameters)):
    parameters = {'EFTflag': 4, 'FullMappingEFTmodel': 3, 
                  'Kmimic': False, 'alphaU': 0.1, 'gammaU': 2.0, 
                  'm': 5.0, 'eps2_0': -0.1, 'gammaA': kMouflage_gammaA_parameters[i],
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    kMouflage_gammaA_models.append(model_temp)
    kMouflage_gammaA_parameters_string.append(parameters)

file_write(kMouflage_gammaA_models, kMouflage_gammaA_parameters_string, kMouflage_gammaA_file_label, kMouflage_gammaA_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(kMouflage_gammaA_models)):
    print("Simulating cosmology: ", kMouflage_gammaA_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(kMouflage_gammaA_models[i], j)
        kMouflage_gammaA_results_gw[i].append(final_res_gw)
        kMouflage_gammaA_results_DE[i].append(final_res_DE)

plotting_function(kMouflage_gammaA_results_gw, kMouflage_gammaA_results_DE, kMouflage_gammaA_figure_label, kMouflage_gammaA_figure_fixed, kMouflage_gammaA_fixed, 
                  iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

#########################################################################################
####### Hu-Sawicki: Linear - PL - Exp.###################################################
#########################################################################################

# Varying fR0 with Omega = 0.05, n = 4
m0 = (2.18e-11)**2                 # Planck Mass squared
m2 = 66.88**2*0.3166               # Equation 3 (arXiv:1601.07536) and Table 2 (arXiv:1807.06209)
c1c2 = 6*(0.6834/0.3166)           # Table 2 (arXiv:1807.06209)
k2rho = 1.228e30*m2                # Equation 5 (arXiv:0705.1158)
R = k2rho+2*c1c2*m2                # Equation 25 (arXiv:0705.1158)

HuSawLin_fR0_fixed = 'HuSawLin_fR0_grav'
HuSawLin_fR0_figure_fixed = 'Hu-Sawicki Gravity'
HuSawLin_fR0_figure_label = ['GR', ]
HuSawLin_fR0_file_label = ['General Relativity', ]

HuSawLin_fR0_math_name = r'$f_{R0} = $'
HuSawLin_fR0_name = 'HuSawLin_fR0_'

HuSawLin_fR0_models = [Standard_GR, ]
HuSawLin_fR0_parameters = [1e-5, 1e-4, 1e-3, 1e-2] #for Omega = 
HuSawLin_fR0_results_gw = np.empty((len(HuSawLin_fR0_parameters)+1, 0)).tolist()
HuSawLin_fR0_results_DE = np.empty((len(HuSawLin_fR0_parameters)+1, 0)).tolist()
HuSawLin_fR0_parameters_string = [ ]
HuSawLin_fR0_figure_label, HuSawLin_fR0_file_label = list_former(HuSawLin_fR0_parameters, HuSawLin_fR0_figure_label, HuSawLin_fR0_file_label, HuSawLin_fR0_math_name, HuSawLin_fR0_name)
    
for i in range(len(HuSawLin_fR0_parameters)):
    parameters = {'EFTflag':2, 'AltParEFTmodel': 2,
                  'OLLambdamodel': 2, 'OLOmegamodel': 3,
                  'OLOmega0': 0.05, 'dark_energy_model' : 'EFTCAMB',
                  'OLLambda0': m0**2/2*(-m2*c1c2 + HuSawLin_fR0_parameters[i]*m2*(m2/R)**4 - R*0.05)
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    HuSawLin_fR0_models.append(model_temp)
    HuSawLin_fR0_parameters_string.append(parameters)
    
print("Iteration: ", iteration)
file_write(HuSawLin_fR0_models, HuSawLin_fR0_parameters_string, HuSawLin_fR0_file_label, HuSawLin_fR0_fixed, LCDM_bckg, iteration)

for i in range(len(HuSawLin_fR0_models)):
    print("Simulating cosmology: ", HuSawLin_fR0_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(HuSawLin_fR0_models[i], j)
        HuSawLin_fR0_results_gw[i].append(final_res_gw)
        HuSawLin_fR0_results_DE[i].append(final_res_DE)

plotting_function(HuSawLin_fR0_results_gw, HuSawLin_fR0_results_DE, HuSawLin_fR0_figure_label, HuSawLin_fR0_figure_fixed, HuSawLin_fR0_fixed, 
                      iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

# Varying n with fR0 = 0.01, Omega = 0.05
HuSawLin_n_fixed = 'HuSawLin_n_grav'
HuSawLin_n_figure_fixed = 'Hu-Sawicki Gravity'
HuSawLin_n_figure_label = ['GR', ]
HuSawLin_n_file_label = ['General Relativity', ]

HuSawLin_n_math_name = '$\Omega = $'
HuSawLin_n_name = 'HuSawLin_n_'

HuSawLin_n_models = [Standard_GR, ]
HuSawLin_n_parameters = [1, 2, 3, 4, 5] #for Omega = 
HuSawLin_n_results_gw = np.empty((len(HuSawLin_n_parameters)+1, 0)).tolist()
HuSawLin_n_results_DE = np.empty((len(HuSawLin_n_parameters)+1, 0)).tolist()
HuSawLin_n_parameters_string = [ ]
HuSawLin_n_figure_label, HuSawLin_n_file_label = list_former(HuSawLin_n_parameters, HuSawLin_n_figure_label, HuSawLin_n_file_label, HuSawLin_n_math_name, HuSawLin_n_name)
    
for i in range(len(HuSawLin_n_parameters)):
    parameters = {'EFTflag':2, 'AltParEFTmodel': 2,
                  'OLLambdamodel': 2, 'OLOmegamodel': 3,
                  'OLOmega0': 0.05, 'dark_energy_model' : 'EFTCAMB',
                  'OLLambda0': m0**2/2*(-m2*c1c2 + 0.01*m2*(m2/R)**(HuSawLin_n_parameters[i])-R*0.05)
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    HuSawLin_n_models.append(model_temp)
    HuSawLin_n_parameters_string.append(parameters)
    
print("Iteration: ", iteration)
file_write(HuSawLin_n_models, HuSawLin_n_parameters_string, HuSawLin_n_file_label, HuSawLin_n_fixed, LCDM_bckg, iteration)

for i in range(len(HuSawLin_n_models)):
    print("Simulating cosmology: ", HuSawLin_n_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(HuSawLin_n_models[i], j)
        HuSawLin_n_results_gw[i].append(final_res_gw)
        HuSawLin_n_results_DE[i].append(final_res_DE)

plotting_function(HuSawLin_n_results_gw, HuSawLin_n_results_DE, HuSawLin_n_figure_label, HuSawLin_n_figure_fixed, HuSawLin_n_fixed, 
                      iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)


# Varying Omega with fR0 = 0.01, n = 4
fR = -m2*c1c2 + 0.01*m2*(m2/R)**4  # Equation 4 (arXiv:1601.07536) 

HuSawLin_Omega_fixed = 'HuSawLin_Omega_grav'
HuSawLin_Omega_figure_fixed = 'Hu-Sawicki Gravity'
HuSawLin_Omega_figure_label = ['GR', ]
HuSawLin_Omega_file_label = ['General Relativity', ]

HuSawLin_Omega_math_name = '$\Omega = $'
HuSawLin_Omega_name = 'HuSawLin_Omega_'

HuSawLin_Omega_models = [Standard_GR, ]
HuSawLin_Omega_parameters = [0.005, 0.025, 0.05, 0.075] #for Omega = 
HuSawLin_Omega_results_gw = np.empty((len(HuSawLin_Omega_parameters)+1, 0)).tolist()
HuSawLin_Omega_results_DE = np.empty((len(HuSawLin_Omega_parameters)+1, 0)).tolist()
HuSawLin_Omega_parameters_string = [ ]
HuSawLin_Omega_figure_label, HuSawLin_Omega_file_label = list_former(HuSawLin_Omega_parameters, HuSawLin_Omega_figure_label, HuSawLin_Omega_file_label, HuSawLin_Omega_math_name, HuSawLin_Omega_name)
    
for i in range(len(HuSawLin_Omega_parameters)):
    parameters = {'EFTflag':2, 'dark_energy_model' : 'EFTCAMB',
                  'AltParEFTmodel': 2, 'OLLambdamodel': 2, 'OLOmegamodel': 3,
                  'OLOmega0': HuSawLin_Omega_parameters[i],
                  'OLLambda0': m0**2/2*(fR-R*HuSawLin_Omega_parameters[i])
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    HuSawLin_Omega_models.append(model_temp)
    HuSawLin_Omega_parameters_string.append(parameters)
for i in range(len(HuSawLin_Omega_parameters)):
    print(m0**2/2*(fR-R*HuSawLin_Omega_parameters[i]))

print("Iteration: ", iteration)
file_write(HuSawLin_Omega_models, HuSawLin_Omega_parameters_string, HuSawLin_Omega_file_label, HuSawLin_Omega_fixed, LCDM_bckg, iteration)

for i in range(len(HuSawLin_Omega_models)):
    print("Simulating cosmology: ", HuSawLin_Omega_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(HuSawLin_Omega_models[i], j)
        HuSawLin_Omega_results_gw[i].append(final_res_gw)
        HuSawLin_Omega_results_DE[i].append(final_res_DE)

plotting_function(HuSawLin_Omega_results_gw, HuSawLin_Omega_results_DE, HuSawLin_Omega_figure_label, HuSawLin_Omega_figure_fixed, HuSawLin_Omega_fixed, 
                      iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

#########################################################################################
# PL Omega + EXP (Slide 22 - https://indico.cern.ch/event/438475/contributions/1090964/attachments/1150896/1652058/DE_Tue_Raveri__0108.pdf)
#########################################################################################

Omega_PL_fixed = 'Omega_PL1'
Omega_PL_figure_fixed = '$\Omega(a) = \Omega_0 a^s$'
Omega_PL_figure_label = ['GR', ]
Omega_PL_file_label = ['General Relativity', ]

Omega_PL_math_name = '$s = $'
Omega_PL_name = 'Omega0_PL1_'

Omega_PL_models = [Standard_GR, ]
Omega_PL_parameters = [0, 0.7, 1.4]
Omega_PL_results_gw = np.empty((len(Omega_PL_parameters)+1, 0)).tolist()
Omega_PL_results_DE = np.empty((len(Omega_PL_parameters)+1, 0)).tolist()
Omega_PL_parameters_string = [ ]
Omega_PL_figure_label, Omega_PL_file_label = list_former(Omega_PL_parameters, Omega_PL_figure_label, 
                                                           Omega_PL_file_label, Omega_PL_math_name, Omega_PL_name)

for i in range(len(Omega_PL_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':3,'EFTOmega0': 0.05,
                  'EFTOmegaExp': Omega_PL_parameters[i]
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_PL_models.append(model_temp)
    Omega_PL_parameters_string.append(parameters)

file_write(Omega_PL_models, Omega_PL_parameters_string, Omega_PL_file_label, Omega_PL_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_PL_models)):
    print("Simulating cosmology: ", Omega_PL_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_PL_models[i], j)
        Omega_PL_results_gw[i].append(final_res_gw)
        Omega_PL_results_DE[i].append(final_res_DE)

plotting_function(Omega_PL_results_gw, Omega_PL_results_DE, Omega_PL_figure_label, Omega_PL_figure_fixed, 
                  Omega_PL_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)


Omega_Exp_fixed = 'Omega_Exp1'
Omega_Exp_figure_fixed = '$\Omega(a) = \exp{\Omega_0 a^s}-1$'
Omega_Exp_figure_label = ['GR', ]
Omega_Exp_file_label = ['General Relativity', ]

Omega_Exp_math_name = '$s = $'
Omega_Exp_name = 'Omega0_Exp1_'

Omega_Exp_models = [Standard_GR, ]
Omega_Exp_parameters = [0, 0.7, 1.4]
Omega_Exp_results_gw = np.empty((len(Omega_Exp_parameters)+1, 0)).tolist()
Omega_Exp_results_DE = np.empty((len(Omega_Exp_parameters)+1, 0)).tolist()
Omega_Exp_parameters_string = [ ]
Omega_Exp_figure_label, Omega_Exp_file_label = list_former(Omega_Exp_parameters, Omega_Exp_figure_label, 
                                                           Omega_Exp_file_label, Omega_Exp_math_name, Omega_Exp_name)

for i in range(len(Omega_Exp_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':4,'EFTOmega0': 0.05,
                  'EFTOmegaExp': Omega_Exp_parameters[i]
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_Exp_models.append(model_temp)
    Omega_Exp_parameters_string.append(parameters)

file_write(Omega_Exp_models, Omega_Exp_parameters_string, Omega_Exp_file_label, Omega_Exp_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_Exp_models)):
    print("Simulating cosmology: ", Omega_Exp_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_Exp_models[i], j)
        Omega_Exp_results_gw[i].append(final_res_gw)
        Omega_Exp_results_DE[i].append(final_res_DE)

plotting_function(Omega_Exp_results_gw, Omega_Exp_results_DE, Omega_Exp_figure_label, Omega_Exp_figure_fixed, 
                  Omega_Exp_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)



#########################################################################################
####### Acoustic Dark Energy ############################################################
#########################################################################################

AcousticDE_fixed = 'AcousticDE_cs2_grav'
AcousticDE_figure_fixed = 'AcousticDE Gravity'
AcousticDE_figure_label = ['GR', ]
AcousticDE_file_label = ['General Relativity', ]

AcousticDE_math_name = r'$c_s^2 = $'
AcousticDE_name = 'AcousticDE_cs2_'

AcousticDE_models = [Standard_GR, ]
AcousticDE_cs2 = [0.001, 0.01, 0.1]
AcousticDE_results_gw = np.empty((len(AcousticDE_cs2)+1, 0)).tolist()
AcousticDE_results_DE = np.empty((len(AcousticDE_cs2)+1, 0)).tolist()
AcousticDE_parameters_string = [ ]
AcousticDE_figure_label, AcousticDE_file_label = list_former(AcousticDE_cs2, AcousticDE_figure_label, AcousticDE_file_label, 
                                                             AcousticDE_math_name, AcousticDE_name)
    
for i in range(len(AcousticDE_cs2)):
    parameters = {'EFTflag': 4, 
                  'FullMappingEFTmodel': 2, 'dark_energy_model' : 'EFTCAMB',
                  'cs2': AcousticDE_cs2[i], 'Log_ac': -0.3, 
                  'f_ac': 0.2, 'p': 0.2, 'wf': 0.2
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    AcousticDE_models.append(model_temp)
    AcousticDE_parameters_string.append(parameters)

file_write(AcousticDE_models, AcousticDE_parameters_string, AcousticDE_file_label, AcousticDE_fixed, LCDM_bckg, iteration)

for i in range(len(AcousticDE_models)):
    print("Simulating cosmology: ", AcousticDE_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(AcousticDE_models[i], j)
        AcousticDE_results_gw[i].append(final_res_gw)
        AcousticDE_results_DE[i].append(final_res_DE)

plotting_function(AcousticDE_results_gw, AcousticDE_results_DE, AcousticDE_figure_label, AcousticDE_figure_fixed, 
                  AcousticDE_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

AcousticDE_fixed = 'AcousticDE_log_grav'
AcousticDE_figure_fixed = 'AcousticDE Gravity'
AcousticDE_figure_label = ['GR', ]
AcousticDE_file_label = ['General Relativity', ]

AcousticDE_math_name = r'$\log{(a_c)} = $'
AcousticDE_name = 'AcousticDE_log_'

AcousticDE_models = [Standard_GR, ]
AcousticDE_log= [-0.1, -0.3, -0.5]
AcousticDE_results_gw = np.empty((len(AcousticDE_log)+1, 0)).tolist()
AcousticDE_results_DE = np.empty((len(AcousticDE_log)+1, 0)).tolist()
AcousticDE_parameters_string = [ ]
AcousticDE_figure_label, AcousticDE_file_label = list_former(AcousticDE_log, AcousticDE_figure_label, AcousticDE_file_label, 
                                                             AcousticDE_math_name, AcousticDE_name)
    
for i in range(len(AcousticDE_log)):
    parameters = {'EFTflag': 4, 
                  'FullMappingEFTmodel': 2, 'dark_energy_model' : 'EFTCAMB',
                  'cs2': 0.01, 'Log_ac': AcousticDE_log[i], 
                  'f_ac': 0.2, 'p': 0.2, 'wf': 0.2
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    AcousticDE_models.append(model_temp)
    AcousticDE_parameters_string.append(parameters)

file_write(AcousticDE_models, AcousticDE_parameters_string, AcousticDE_file_label, AcousticDE_fixed, LCDM_bckg, iteration)

for i in range(len(AcousticDE_models)):
    print("Simulating cosmology: ", AcousticDE_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(AcousticDE_models[i], j)
        AcousticDE_results_gw[i].append(final_res_gw)
        AcousticDE_results_DE[i].append(final_res_DE)

plotting_function(AcousticDE_results_gw, AcousticDE_results_DE, AcousticDE_figure_label, AcousticDE_figure_fixed, 
                  AcousticDE_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

AcousticDE_fixed = 'AcousticDE_fac_grav'
AcousticDE_figure_fixed = 'AcousticDE Gravity'
AcousticDE_figure_label = ['GR', ]
AcousticDE_file_label = ['General Relativity', ]

AcousticDE_math_name = r'$f_{ac} = $'
AcousticDE_name = 'AcousticDE_fac_'

AcousticDE_models = [Standard_GR, ]
AcousticDE_ac= [0.02, 0.2, 0.5]
AcousticDE_results_gw = np.empty((len(AcousticDE_log)+1, 0)).tolist()
AcousticDE_results_DE = np.empty((len(AcousticDE_log)+1, 0)).tolist()
AcousticDE_parameters_string = [ ]
AcousticDE_figure_label, AcousticDE_file_label = list_former(AcousticDE_log, AcousticDE_figure_label, AcousticDE_file_label, 
                                                             AcousticDE_math_name, AcousticDE_name)
    
for i in range(len(AcousticDE_ac)):
    parameters = {'EFTflag': 4, 
                  'FullMappingEFTmodel': 2, 'dark_energy_model' : 'EFTCAMB',
                  'cs2': 0.01, 'Log_ac': -0.1, 
                  'f_ac': AcousticDE_ac[i], 'p': 0.2, 'wf': 0.2
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    AcousticDE_models.append(model_temp)
    AcousticDE_parameters_string.append(parameters)

file_write(AcousticDE_models, AcousticDE_parameters_string, AcousticDE_file_label, AcousticDE_fixed, LCDM_bckg, iteration)

for i in range(len(AcousticDE_models)):
    print("Simulating cosmology: ", AcousticDE_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(AcousticDE_models[i], j)
        AcousticDE_results_gw[i].append(final_res_gw)
        AcousticDE_results_DE[i].append(final_res_DE)

plotting_function(AcousticDE_results_gw, AcousticDE_results_DE, AcousticDE_figure_label, AcousticDE_figure_fixed, 
                  AcousticDE_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

#########################################################################################
####### f(R) Gravity ####################################################################
####### Five parameters are chosen per list to not clutter plots. #######################
#########################################################################################

fR_fixed = 'fR_grav'
fR_figure_fixed = '$f(R)$ Gravity'
fR_figure_label = ['GR', ]
fR_file_label = ['General Relativity', ]

fR_math_name = '$B_0 = $'
fR_name = 'f(R)_B_'

fR_models = [Standard_GR, ]
B0_parameters = np.linspace(1, 4, 4)
fR_results_gw = np.empty((len(B0_parameters)+1, 0)).tolist()
fR_results_DE = np.empty((len(B0_parameters)+1, 0)).tolist()
fR_parameters_string = [ ]
fR_figure_label, fR_file_label = list_former(B0_parameters, fR_figure_label, fR_file_label, fR_math_name, fR_name)
    
for i in range(len(B0_parameters)):
    parameters = {'EFTflag' :3,
                  'DesignerEFTmodel' :1,
                  'EFTwDE': 0, 'dark_energy_model' : 'EFTCAMB',
                  'EFTB0': B0_parameters[i],
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    fR_models.append(model_temp)
    fR_parameters_string.append(parameters)
    
print("Iteration: ", iteration)
file_write(fR_models, fR_parameters_string, fR_file_label, fR_fixed, LCDM_bckg, iteration)

for i in range(len(fR_models)):
    print("Simulating cosmology: ", fR_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(fR_models[i], j)
        fR_results_gw[i].append(final_res_gw)
        fR_results_DE[i].append(final_res_DE)

plotting_function(fR_results_gw, fR_results_DE, fR_figure_label, fR_figure_fixed, fR_fixed, 
                      iteration, LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)

#########################################################################################
####### Linear Omega (Negative) #########################################################
#########################################################################################

Omega_Neg_Lin_fixed = 'Omega_Neg_Lin'
Omega_Neg_Lin_figure_fixed = '$\Omega(a) = \Omega_0 a$'
Omega_Neg_Lin_figure_label = ['GR', ]
Omega_Neg_Lin_file_label = ['General Relativity', ]

Omega_Neg_Lin_math_name = '$\Omega_0 = $'
Omega_Neg_Lin_name = 'Omega0_Neg_Lin_'

Omega_Neg_Lin_models = [Standard_GR, ]
Omega_Neg_Lin_parameters = [-0.012, -0.006, 0.006, 0.012]
Omega_Neg_Lin_results_gw = np.empty((len(Omega_Neg_Lin_parameters)+1, 0)).tolist()
Omega_Neg_Lin_results_DE = np.empty((len(Omega_Neg_Lin_parameters)+1, 0)).tolist()
Omega_Neg_Lin_parameters_string = [ ]
Omega_Neg_Lin_figure_label, Omega_Neg_Lin_file_label = list_former(Omega_Neg_Lin_parameters, Omega_Neg_Lin_figure_label, 
                                                                       Omega_Neg_Lin_file_label, Omega_Neg_Lin_math_name, Omega_Neg_Lin_name)

for i in range(len(Omega_Neg_Lin_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': Omega_Neg_Lin_parameters[i],
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_Neg_Lin_models.append(model_temp)
    Omega_Neg_Lin_parameters_string.append(parameters)

file_write(Omega_Neg_Lin_models, Omega_Neg_Lin_parameters_string, Omega_Neg_Lin_file_label, Omega_Neg_Lin_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_Neg_Lin_models)):
    print("Simulating cosmology: ", Omega_Neg_Lin_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_Neg_Lin_models[i], j)
        Omega_Neg_Lin_results_gw[i].append(final_res_gw)
        Omega_Neg_Lin_results_DE[i].append(final_res_DE)

plotting_function(Omega_Neg_Lin_results_gw, Omega_Neg_Lin_results_DE, Omega_Neg_Lin_figure_label, Omega_Neg_Lin_figure_fixed, 
                  Omega_Neg_Lin_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)

#########################################################################################
# Linear Omega + Constant Gamma1 ########################################################
#########################################################################################

LinO_ConstGamma1_fixed = 'LinO_ConstGamma1'
LinO_ConstGamma1_figure_fixed = '$\Omega(a) = \Omega_0 a, \gamma_1 (a) = \gamma_1^0$\n                      $\Omega_0 = 0.05$'
LinO_ConstGamma1_figure_label = ['GR', ]
LinO_ConstGamma1_file_label = ['General Relativity', ]

LinO_ConstGamma1_math_name = '$\gamma_1^0 = $'
LinO_ConstGamma1_name = 'Lin0_LinGamma1_'

LinO_ConstGamma1_models = [Standard_GR, ]
LinO_ConstGamma1_parameters = [1e-3, 1e-2, 5e-2]
LinO_ConstGamma1_results_gw = np.empty((len(LinO_ConstGamma1_parameters)+1, 0)).tolist()
LinO_ConstGamma1_results_DE = np.empty((len(LinO_ConstGamma1_parameters)+1, 0)).tolist()
LinO_ConstGamma1_parameters_string = [ ]
LinO_ConstGamma1_figure_label, LinO_ConstGamma1_file_label = list_former(LinO_ConstGamma1_parameters, LinO_ConstGamma1_figure_label, 
                                                                         LinO_ConstGamma1_file_label, LinO_ConstGamma1_math_name, 
                                                                         LinO_ConstGamma1_name
                                                                         )

for i in range(len(LinO_ConstGamma1_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': 0.05,
                  'PureEFTmodelGamma1':1, 'EFTGamma10': LinO_ConstGamma1_parameters[i],
                  }
    parameters.update(stability_flag)
    model_temp = model_parameters(parameters)
    LinO_ConstGamma1_models.append(model_temp)
    LinO_ConstGamma1_parameters_string.append(parameters)

file_write(LinO_ConstGamma1_models, LinO_ConstGamma1_parameters_string, LinO_ConstGamma1_file_label, 
           LinO_ConstGamma1_fixed, LCDM_bckg, iteration)

for i in range(len(LinO_ConstGamma1_models)):
    print("Simulating cosmology: ", LinO_ConstGamma1_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(LinO_ConstGamma1_models[i], j)
        LinO_ConstGamma1_results_gw[i].append(final_res_gw)
        LinO_ConstGamma1_results_DE[i].append(final_res_DE)

plotting_function(LinO_ConstGamma1_results_gw, LinO_ConstGamma1_results_DE, LinO_ConstGamma1_figure_label, 
                  LinO_ConstGamma1_figure_fixed, LinO_ConstGamma1_fixed, iteration, 
                  LCDM_bckg, ell, 0.04, 0.04, 0.04, 0.04, 0.55, 0.75
                  )

#########################################################################################
# Linear Omega + Linear Gamma1 ##########################################################
#########################################################################################

LinO_LinGamma1_fixed = 'LinO_LinGamma1'
LinO_LinGamma1_figure_fixed = '$\Omega(a) = \Omega_0 a, \gamma_1 (a) = \gamma_1^0 a$\n                     $\Omega_0 = 0.05$'
LinO_LinGamma1_figure_label = ['GR', ]
LinO_LinGamma1_file_label = ['General Relativity', ]

LinO_LinGamma1_math_name = '$\gamma_1^0 = $'
LinO_LinGamma1_name = 'Lin0_LinGamma1_'

LinO_LinGamma1_models = [Standard_GR, ]
LinO_LinGamma1_parameters = [1e-3, 1e-2, 5e-2]
LinO_LinGamma1_results_gw = np.empty((len(LinO_LinGamma1_parameters)+1, 0)).tolist()
LinO_LinGamma1_results_DE = np.empty((len(LinO_LinGamma1_parameters)+1, 0)).tolist()
LinO_LinGamma1_parameters_string = [ ]
LinO_LinGamma1_figure_label, LinO_LinGamma1_file_label = list_former(LinO_LinGamma1_parameters, LinO_LinGamma1_figure_label, 
                                                                     LinO_LinGamma1_file_label, LinO_LinGamma1_math_name, 
                                                                     LinO_LinGamma1_name
                                                                     )

for i in range(len(LinO_LinGamma1_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': 0.05,
                  'PureEFTmodelGamma1':2, 'EFTGamma10': LinO_LinGamma1_parameters[i],
                  }
    parameters.update(stability_flag)
    model_temp = model_parameters(parameters)
    LinO_LinGamma1_models.append(model_temp)
    LinO_LinGamma1_parameters_string.append(parameters)

file_write(LinO_LinGamma1_models, LinO_LinGamma1_parameters_string, LinO_LinGamma1_file_label, 
           LinO_LinGamma1_fixed, LCDM_bckg, iteration
           )

for i in range(len(LinO_LinGamma1_models)):
    print("Simulating cosmology: ", LinO_LinGamma1_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(LinO_LinGamma1_models[i], j)
        LinO_LinGamma1_results_gw[i].append(final_res_gw)
        LinO_LinGamma1_results_DE[i].append(final_res_DE)

plotting_function(LinO_LinGamma1_results_gw, LinO_LinGamma1_results_DE, LinO_LinGamma1_figure_label, 
                  LinO_LinGamma1_figure_fixed, LinO_LinGamma1_fixed, iteration, LCDM_bckg, ell, 
                  0.04, 0.04, 0.04, 0.04, 0.55, 0.75)


#########################################################################################
# Linear Omega + Constant Gamma2 ########################################################
#########################################################################################

LinO_ConstGamma2_fixed = 'LinO_ConstGamma2'
LinO_ConstGamma2_figure_fixed = '$\Omega(a) = \Omega_0 a, \gamma_2 (a) = \gamma_2^0$'
LinO_ConstGamma2_figure_label = ['GR', ]
LinO_ConstGamma2_file_label = ['General Relativity', ]

LinO_ConstGamma2_math_name = '$\gamma_2^0 = $'
LinO_ConstGamma2_name = 'Lin0_LinGamma2_'

LinO_ConstGamma2_models = [Standard_GR, ]
LinO_ConstGamma2_parameters = [1e-3, 1e-2, 5e-2]
LinO_ConstGamma2_results_gw = np.empty((len(LinO_ConstGamma2_parameters)+1, 0)).tolist()
LinO_ConstGamma2_results_DE = np.empty((len(LinO_ConstGamma2_parameters)+1, 0)).tolist()
LinO_ConstGamma2_parameters_string = [ ]
LinO_ConstGamma2_figure_label, LinO_ConstGamma2_file_label = list_former(LinO_ConstGamma2_parameters, LinO_ConstGamma2_figure_label, 
                                                                         LinO_ConstGamma2_file_label, LinO_ConstGamma2_math_name, 
                                                                         LinO_ConstGamma2_name
                                                                         )

for i in range(len(LinO_ConstGamma2_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': 0.05,
                  'PureEFTmodelGamma2':1, 'EFTGamma10': LinO_ConstGamma2_parameters[i],
                  }
    parameters.update(stability_flag)
    model_temp = model_parameters(parameters)
    LinO_ConstGamma2_models.append(model_temp)
    LinO_ConstGamma2_parameters_string.append(parameters)

file_write(LinO_ConstGamma2_models, LinO_ConstGamma2_parameters_string, LinO_ConstGamma2_file_label, 
           LinO_ConstGamma2_fixed, LCDM_bckg, iteration
           )

for i in range(len(LinO_ConstGamma2_models)):
    print("Simulating cosmology: ", LinO_ConstGamma2_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(LinO_ConstGamma2_models[i], j)
        LinO_ConstGamma2_results_gw[i].append(final_res_gw)
        LinO_ConstGamma2_results_DE[i].append(final_res_DE)

plotting_function(LinO_ConstGamma2_results_gw, LinO_ConstGamma2_results_DE, LinO_ConstGamma2_figure_label, 
                  LinO_ConstGamma2_figure_fixed, LinO_ConstGamma2_fixed, iteration, 
                  LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)

#########################################################################################
# Linear Omega + Linear Gamma2 ##########################################################
#########################################################################################

LinO_LinGamma2_fixed = 'LinO_LinGamma2'
LinO_LinGamma2_figure_fixed = '$\Omega(a) = \Omega_0 a, \gamma_2 (a) = \gamma_2^0 a$'
LinO_LinGamma2_figure_label = ['GR', ]
LinO_LinGamma2_file_label = ['General Relativity', ]

LinO_LinGamma2_math_name = '$\gamma_1^0 = $'
LinO_LinGamma2_name = 'Lin0_LinGamma2_'

LinO_LinGamma2_models = [Standard_GR, ]
LinO_LinGamma2_parameters = [1e-3, 1e-2, 5e-2]
LinO_LinGamma2_results_gw = np.empty((len(LinO_LinGamma2_parameters)+1, 0)).tolist()
LinO_LinGamma2_results_DE = np.empty((len(LinO_LinGamma2_parameters)+1, 0)).tolist()
LinO_LinGamma2_parameters_string = [ ]
LinO_LinGamma2_figure_label, LinO_LinGamma2_file_label = list_former(LinO_LinGamma2_parameters, LinO_LinGamma2_figure_label, 
                                                                     LinO_LinGamma2_file_label, LinO_LinGamma2_math_name, 
                                                                     LinO_LinGamma2_name
                                                                     )

for i in range(len(LinO_LinGamma2_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':2,'EFTOmega0': 0.05,
                  'PureEFTmodelGamma2':2, 'EFTGamma10': LinO_LinGamma2_parameters[i],
                  }
    parameters.update(stability_flag)
    model_temp = model_parameters(parameters)
    LinO_LinGamma2_models.append(model_temp)
    LinO_LinGamma2_parameters_string.append(parameters)

file_write(LinO_LinGamma2_models, LinO_LinGamma2_parameters_string, LinO_LinGamma2_file_label, 
           LinO_LinGamma2_fixed, LCDM_bckg, iteration
           )

for i in range(len(LinO_LinGamma2_models)):
    print("Simulating cosmology: ", LinO_LinGamma2_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(LinO_LinGamma2_models[i], j)
        LinO_LinGamma2_results_gw[i].append(final_res_gw)
        LinO_LinGamma2_results_DE[i].append(final_res_DE)

plotting_function(LinO_LinGamma2_results_gw, LinO_LinGamma2_results_DE, LinO_LinGamma2_figure_label, 
                  LinO_LinGamma2_figure_fixed, LinO_LinGamma2_fixed, iteration, LCDM_bckg, 
                  ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)


#########################################################################################
# PL Omega
#########################################################################################

Omega_PL_fixed = 'Omega_PL'
Omega_PL_figure_fixed = '$\Omega(a) = \Omega_0 a^s$'
Omega_PL_figure_label = ['GR', ]
Omega_PL_file_label = ['General Relativity', ]

Omega_PL_math_name = '$s = $'
Omega_PL_name = 'Omega0_PL_'

Omega_PL_models = [Standard_GR, ]
Omega_PL_parameters = [0.1, 1, 10]
Omega_PL_results_gw = np.empty((len(Omega_PL_parameters)+1, 0)).tolist()
Omega_PL_results_DE = np.empty((len(Omega_PL_parameters)+1, 0)).tolist()
Omega_PL_parameters_string = [ ]
Omega_PL_figure_label, Omega_PL_file_label = list_former(Omega_PL_parameters, Omega_PL_figure_label, 
                                                           Omega_PL_file_label, Omega_PL_math_name, Omega_PL_name)

for i in range(len(Omega_PL_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':3,'EFTOmega0': 0.05,
                  'EFTOmegaExp': Omega_PL_parameters[i]
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_PL_models.append(model_temp)
    Omega_PL_parameters_string.append(parameters)

file_write(Omega_PL_models, Omega_PL_parameters_string, Omega_PL_file_label, Omega_PL_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_PL_models)):
    print("Simulating cosmology: ", Omega_PL_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_PL_models[i], j)
        Omega_PL_results_gw[i].append(final_res_gw)
        Omega_PL_results_DE[i].append(final_res_DE)

plotting_function(Omega_PL_results_gw, Omega_PL_results_DE, Omega_PL_figure_label, Omega_PL_figure_fixed, 
                  Omega_PL_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)

#########################################################################################
# Exp. Omega
#########################################################################################
Omega_Exp_fixed = 'Omega_Exp'
Omega_Exp_figure_fixed = '$\Omega(a) = \exp{\Omega_0 a^s}-1$'
Omega_Exp_figure_label = ['GR', ]
Omega_Exp_file_label = ['General Relativity', ]

Omega_Exp_math_name = '$s = $'
Omega_Exp_name = 'Omega0_Exp_'

Omega_Exp_models = [Standard_GR, ]
Omega_Exp_parameters = [0.1, 1.0, 10]
Omega_Exp_results_gw = np.empty((len(Omega_Exp_parameters)+1, 0)).tolist()
Omega_Exp_results_DE = np.empty((len(Omega_Exp_parameters)+1, 0)).tolist()
Omega_Exp_parameters_string = [ ]
Omega_Exp_figure_label, Omega_Exp_file_label = list_former(Omega_Exp_parameters, Omega_Exp_figure_label, 
                                                           Omega_Exp_file_label, Omega_Exp_math_name, Omega_Exp_name)

for i in range(len(Omega_Exp_parameters)):
    parameters = {'EFTflag':1,
                  'PureEFTmodel':1, 'EFTwDE':0, 'dark_energy_model' : 'EFTCAMB',
                  'PureEFTmodelOmega':4,'EFTOmega0': 0.05,
                  'EFTOmegaExp': Omega_Exp_parameters[i]
                  }
    parameters.update(stability_flag) 
    model_temp = model_parameters(parameters)
    Omega_Exp_models.append(model_temp)
    Omega_Exp_parameters_string.append(parameters)

file_write(Omega_Exp_models, Omega_Exp_parameters_string, Omega_Exp_file_label, Omega_Exp_fixed, LCDM_bckg, iteration)
print("Iteration: ", iteration)

for i in range(len(Omega_Exp_models)):
    print("Simulating cosmology: ", Omega_Exp_file_label[i])
    for j in range(2):
        final_res_gw, final_res_DE, ell = model_simulation(Omega_Exp_models[i], j)
        Omega_Exp_results_gw[i].append(final_res_gw)
        Omega_Exp_results_DE[i].append(final_res_DE)

plotting_function(Omega_Exp_results_gw, Omega_Exp_results_DE, Omega_Exp_figure_label, Omega_Exp_figure_fixed, 
                  Omega_Exp_fixed, iteration, LCDM_bckg, ell, 0.04, 0.04, 0.55, 0.04, 0.55, 0.75)


