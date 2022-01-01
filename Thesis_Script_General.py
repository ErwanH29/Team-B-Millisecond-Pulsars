import sys, os, copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import utilities.utilities as utilities
import utilities.color_utilities as col
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
                                mnu=0.12, tau=0.0522,
                                num_massive_neutrinos=1,
                                nnu=3.046, EFTCAMB_params = param
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
    textfile = open("Thesis_Data/Parameters_"+str(background)+"_"+str(model_fixed)+"_"+str(it)+".txt", "w")

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

def plotting_function(results_gw, results_DE, figure_label, figure_fixed, fixed, file_label, iteration, ylim, bckg, ell):
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

        ax1.set_ylim([1e-11, 1e-4])
        ax2.set_ylim([ylim, 1e-4])
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
                    fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                    borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

        ax1.text(0.04, 0.04, "Total Signal:   "+str(figure_fixed)+'\nRedshift:         '+str(redshift_labels[i]), 
                 verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                 fontsize=main_fontsize, bbox=dict(facecolor='white', boxstyle='square',edgecolor='white')) 

        ax2.text(0.04, 0.04, 'Dark Energy Clustering', verticalalignment='bottom',  horizontalalignment='left', 
                 transform=ax2.transAxes, fontsize=main_fontsize, bbox=dict(facecolor='white',
                 boxstyle='square',edgecolor='white')) 
                 
        plt.savefig('Thesis_Images/'+str(bckg)+'_Plots/Total_Signal/Total_Signal_'+str(fixed)+"_Plot_GW_"+str(iteration)+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
        plt.close('')

    for i in range(len(redshifts)): #Plot total signal (GW-SN) of all models for fixed redshifts
        fig = plt.gcf()

        x_size = 30
        y_size = 8.0 +0.5
        fig.set_size_inches( x_size/2.7, y_size/2.7 )
        gs = gridspec.GridSpec(1,2)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax1.set_ylim([ylim, 1e-4])
        ax2.set_ylim([ylim, 1e-4])
        ax1.set_ylabel(r'$\ell(\ell+1)C^{\Delta\varphi}_\ell/2\pi$', fontsize=main_fontsize)

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
            index_maximum = np.where(np.max(results_DE[j][0][0,0,:,i,i]) == results_DE[j][0][0,0,:,i,i])
            ax1.plot(ell, results_DE[j][0][0,0,:,i,i], lw=1., color = colors[j], linestyle = '-', label = figure_label[j])
            ax1.scatter(ell[index_maximum], results_DE[j][0][0,0,index_maximum,i,i], color = colors[j], s = 1.5)
            ax2.plot(ell, results_DE[j][0][0,0,:,i,i], lw=1., color = colors[j], linestyle = '-')
            ax2.plot(ell, results_DE[j][1][0,0,:,i,i], lw=1., color = colors[j], linestyle = ":")                       

        fig.legend( labels = figure_label, fontsize = 0.75*main_fontsize, frameon   = True,
                    fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                    borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

        ax1.text(0.04, 0.04, "Total Signal:   "+str(figure_fixed)+'\n      Redshift:   '+str(redshift_labels[i]), 
                 verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                 fontsize=main_fontsize, bbox=dict(facecolor='white', boxstyle='square',edgecolor='white')) 

        ax2.text(0.04, 0.04, 'Dark Energy Clustering', verticalalignment='bottom',  horizontalalignment='left', 
                 transform=ax2.transAxes, fontsize=main_fontsize, bbox=dict(facecolor='white',
                 boxstyle='square',edgecolor='white')) 
                 
        plt.savefig('Thesis_Images/'+str(bckg)+'_Plots/Total_Signal/Total_Signal_'+str(fixed)+"_Plot_GW-SN_"+str(iteration)+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
        plt.close('')

    for i in redshift_vals: #Plotting the total and individual effects of each model for defined redshifts
        for j in range(len(figure_label)): #Number of models plotting
            fig = plt.gcf()

            x_size = 30
            y_size = 8.0 +0.5
            fig.set_size_inches( x_size/2.7, y_size/2.7 )
            gs = gridspec.GridSpec(1,2)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0,1])
            ax1.set_ylabel(r'$\ell(\ell+1)C^{GW}_\ell/2\pi$', fontsize=main_fontsize)
            ax2.set_ylabel(r'$\ell(\ell+1)C^{\Delta\varphi}_\ell/2\pi$', fontsize=main_fontsize)
            ax1.set_ylim([1e-12, 1e-4])
            ax2.set_ylim([ylim, 1e-4])

            for _ax in [ax1, ax2]:        
                _ax.set_xscale('log')
                _ax.set_yscale('log')
                _ax.set_xlim([2, 1000])
                _ax.set_xlabel(r'Monopole Number [$\ell$]', fontsize=main_fontsize, labelpad=4)

                ticks  = [ 2, 10, 100, 1000 ]
                _ax.set_xticks(ticks);
                _ax.set_xticklabels( ticks, horizontalalignment='center', fontsize=0.9*main_fontsize);
                _ax.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right');
                _ax.xaxis.get_majorticklabels()[0].set_horizontalalignment('left');
                _ax.tick_params(axis='both', which='both', direction='in',
                                    bottom=True, top=True, left=True, right=True)

            for k in range(5): #Plotting individual effects
                ax1.plot(ell, results_gw[j][k][0,0,:,i,i], lw=1., color = colors[k], linestyle = line[k], label = sim_effect[k])

                ax1.text(0.04, 0.04, "Model: "+str(figure_label[j])+"   Redshift: "+str(redshift_labels[i]), 
                         verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes,
                         fontsize=main_fontsize, bbox=dict(facecolor='white',boxstyle='square',edgecolor='white'))

                ax2.plot(ell, results_DE[j][k][0,0,:,i,i], lw=1., color = colors[k], linestyle = line2[k])
                
            fig.legend( labels = sim_effect, fontsize = 0.75*main_fontsize, frameon   = True,
                        fancybox = False, edgecolor = 'k', ncol = 6, loc = 'upper center',
                        borderaxespad = 0.7, columnspacing = 1.0, handlelength = 1.5)

            plt.savefig('Thesis_Images/'+str(bckg)+'_Plots/Indiv_Effects/Indiv_Effects_Model_'+str(file_label[j])+str(iteration)+'_redshift'+str(rs_fig_name[i])+'.pdf', bbox_inches='tight')
            plt.close()

    return 

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

#########################################################################################
# Global parameters
#########################################################################################

stability_flag = {
                  'EFT_ghost_math_stability'   : False,
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

redshifts = np.array([0.1, 0.2, 0.5, 1.0])
redshift_vals = [0, -1]
redshift_labels = [r'$z = 0.10$', r'$z = 0.20$', r'$z = 0.50$', r'$z = 1.00$']
rs_fig_name = ['_z=0.10', '_z=0.20', '_z=0.50', '_z=1.0']

colors = ['black', 'C3', 'C0', 'C4', 'C1', 'C2']
line = ['-', '--', '--', '--', '--', '--']
line2 = ['-', ':', ':', ':', ':', ':']

redshift_width = 0.01
windows = [ GaussianSourceWindow(redshift = i, source_type='gw', sigma=redshift_width) for i in redshifts ]
windows = windows + [ GaussianSourceWindow(redshift = i, source_type='sn', sigma=redshift_width) for i in redshifts ]

Standard_GR = model_parameters({'EFTflag': 0})
LCDM_bckg = 'LCDM'
wCDM_bckg = 'wCDM'
CPL_bckg = 'CPL_EFT'