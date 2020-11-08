import numpy as np
from statistics import mean, stdev
import os
import ast
import json

class fileworker(object):
    def __init__(self, default_path=True):
        self.sim_files_path = 'sim_results/'
        self.dir_sim = os.listdir(self.sim_files_path)
        self.n_sim =len(self.dir_sim)
        
        self.specific_name = 'neut_star_{}' 
        self.sim = self.sim_files_path + self.specific_name
        
        if default_path == False:
            print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
please specify your path from "gravsolve" to your simulation folders. \n \
the path should have the following form: "path/to/your/simulation_folders_{}", with the "{}" in the string. \n \
it is important that the simulation pickles and text files are in folders indexed as such: \n \
simulation 1: "simulation_folder_1", simulation 2: "simulation_folder_2", etc. \n \
------------------------------------------------------------\
------------------------------------------------------------')
            self.sim = input('give path:')
        
        

    def move_data_to_working(self):
        if self.n_sim == 0:
            print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
no simulation files present. Please generate a new set to create files. \n \
------------------------------------------------------------\
------------------------------------------------------------')
            
        else:
            i = input('you have {} different simulations to choose from. \n \
please enter the number of the simulation folder you would like to work with:'.format(self.n_sim))
            sim_select = self.sim.format(i)
            os.system('cp -a {}/. working_data'.format(sim_select))
        print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
loaded data from: "{}" \n \
------------------------------------------------------------\
------------------------------------------------------------'.format(sim_select))
                          
    def move_data_to_sim_results(self):
        self.n_sim += 1
        i = self.n_sim
        sim_select = self.sim.format(i)
        os.system('mkdir {}'.format(sim_select))
        os.system('cp -a working_data/. {}/'.format(sim_select))
        print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
data has been saved to: "{}" \n \
------------------------------------------------------------\
------------------------------------------------------------'.format(sim_select))

    def clear_working_data(self):
            os.system('rm working_data/*')
            print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
"working_data" directory has been cleared \n \
------------------------------------------------------------\
------------------------------------------------------------')
            
    def get_enc_rate(self):
        tot_enc_r = []
        tot_total = []
        pos_dist_total = []
        if self.n_sim == 0:
            print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
no simulation files present. Please generate a new set to create files. \n \
------------------------------------------------------------\
------------------------------------------------------------')
        else:
            for i in range(self.n_sim):
                sim_select = self.sim.format(i+1)
                check = open('{}/check.txt'.format(sim_select))
                check = check.read()
                check = check.strip()
                check = ast.literal_eval(check)
                
                total = open('{}/total.txt'.format(sim_select))
                total = total.read()
                total = total.strip()
                total = ast.literal_eval(total)
                enc_r = len(check) / total
                tot_total.append(total)
                tot_enc_r.append(enc_r)
                
                pos_dist = open('{}/distances.txt'.format(sim_select))
                pos_dist = pos_dist.read()
                pos_dist = pos_dist.strip()
                pos_dist = ast.literal_eval(pos_dist)
                pos_dist_total.extend(pos_dist)
                    
            mean_r = mean(tot_enc_r)
            stdev_r = stdev(tot_enc_r)
            pos_dist_total.sort()
            
            mean_total = mean(tot_total)
            stdev_total = stdev(tot_total)
            s = open('statistics/encounter_rate_{}sims.txt'.format(self.n_sim), 'w')
            string = 'mean close encounter rate = {} and standard deviation = {} for {} simulations. \
Mean amount of neutron stars formed = {} and standard deviation = {}'.format(mean_r, stdev_r, self.n_sim, mean_total, stdev_total)
            json.dump(string, s)
            s.close()
            p = open('statistics/position_distribution_{}sims.txt'.format(self.n_sim), 'w')
            json.dump(pos_dist_total, p)
            p.close()
            print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
close encounter rate statistics have been saved in: statistics/encounter_rate_{}sims.txt \n \
position distribution statistics have been saved in: statistics/position_distribution_{}.txt \n \
------------------------------------------------------------\
------------------------------------------------------------'.format(self.n_sim, self.n_sim))
        return self.n_sim