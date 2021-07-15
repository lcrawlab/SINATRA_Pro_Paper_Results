#!/bin/python3

import numpy as np
import MDAnalysis as mda
from sklearn.linear_model import SGDClassifier
import os
from matplotlib import pyplot as plt

directory_original = '../WT_R164S_65_213/'
n_frame = 100

directories_perturb = [
                       'simulation_perturb_omega_loop_vector_0.5_0.5_0.5',
                       'simulation_perturb_omega_loop_vector_1.0_1.0_1.0',
                       'simulation_perturb_omega_loop_vector_2.0_2.0_2.0',
                       'simulation_perturb_omega_loop_sphere_0.5',
                       'simulation_perturb_omega_loop_sphere_1.0',
                       'simulation_perturb_omega_loop_sphere_2.0',
                      ]

filenames = [
             'perturb_0.5_vector',
             'perturb_1.0_vector',
             'perturb_2.0_vector',
             'perturb_0.5_sphere',
             'perturb_1.0_sphere',
             'perturb_2.0_sphere',
            ]

for offset in range(0,50,10):
    for directory_perturb, filename in zip(directories_perturb,filenames):
        print(filename, offset)
        positions = []
        y = []
        for i_frame in range(n_frame):
            pdbfile = directory_perturb + '/pdb/WT_offset_%d/'%offset + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            protein = u.select_atoms('protein')
            protein_pos = protein.positions
            positions.append(protein_pos.flatten())
            y.append(0)

        for i_frame in range(n_frame):
            pdbfile = directory_perturb + '/pdb/offset_%d_perturbed/'%(offset+50) + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            protein = u.select_atoms('protein')
            protein_pos = protein.positions
            positions.append(protein_pos.flatten())
            y.append(1)

        positions = np.array(positions) 
        en = SGDClassifier(loss="log", penalty="elasticnet")
        en.fit(positions,y)
        np.savetxt('elasticnetlog_coef_%s_offset_%d_%d.txt'%(filename,offset,offset+50),en.coef_[0])

