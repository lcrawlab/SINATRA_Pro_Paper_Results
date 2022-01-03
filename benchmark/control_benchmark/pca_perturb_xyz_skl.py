#!/bin/python3

import MDAnalysis as mda
import numpy as np
import os
from sklearn.decomposition import PCA

directory_original = '../WT_R164S_65_213/'
n_frame = 100
"""
directories_perturb = ['../simulation_perturb_omega_loop_new_0.5_0.5_0.0/',
                       '../simulation_perturb_omega_loop_new_0.5_0.5_0.5/',
                       '../simulation_perturb_omega_loop_new_0.5_0.1_0.5/',
                       '../simulation_perturb_omega_loop_vector_0.5_0.5_0.5_0.5/',
                       '../simulation_perturb_omega_loop_vector_0.5_0.1_0.1_0.1/']

filenames = ['perturb_0.5_0.5_0.0',
             'perturb_0.5_0.5_0.5',
             'perturb_0.5_0.1_0.5',
             'perturb_0.5_0.5_vector',
             'perturb_0.5_0.1_vector']


labels = ['r += 0.5 $\mathrm{\AA}$',
          'r += 0.5 $\pm$ 0.5 $\mathrm{\AA}$',
          'r += 0.1 $\pm$ 0.5 $\mathrm{\AA}$',
          'xyz += [0.5,0.5,0.5] $\mathrm{\AA}$',
          'xyz += [0.1,0.1,0.1] $\mathrm{\AA}$']
"""

directory = "../../python_script/"
directories_perturb = ['simulation_perturb_omega_loop_vector_0.0_0.5_0.5_0.5_offset',
                       'simulation_perturb_omega_loop_vector_0.0_1.0_1.0_1.0_offset',
                       'simulation_perturb_omega_loop_vector_0.0_2.0_2.0_2.0_offset',
                       'simulation_perturb_omega_loop_sphere_0.5_0.0_offset',
                       'simulation_perturb_omega_loop_sphere_1.0_0.0_offset',
                       'simulation_perturb_omega_loop_sphere_2.0_0.0_offset']


filenames = [
             'perturb_0.0_0.5_vector',
             'perturb_0.0_1.0_vector',
             'perturb_0.0_2.0_vector',
             'perturb_0.5_0.0_sphere',
             'perturb_1.0_0.0_sphere',
             'perturb_2.0_0.0_sphere',
            ]

n_components = 100

for offset in range(0,50,10):

    positions = []
    for i_frame in range(n_frame):
        pdbfile = directory + directories_perturb[0] + '/pdb/WT_offset_%d/'%offset + 'WT_frame%d.pdb'%i_frame
        u = mda.Universe(pdbfile)
        #CA = u.select_atoms('name CA')
        CA = u.select_atoms('protein')
        CA_pos = CA.positions
        CA_pos -= CA.center_of_mass()
        positions.append(CA_pos.flatten())
    positions = np.array(positions)
    pca = PCA(n_components=n_components)
    pca.fit(positions)
    w = pca.singular_values_
    v = pca.components_
    ev = pca.explained_variance_
    np.savetxt('pca_aa_original_%d_w_offset_%d.txt'%(n_components,offset),w)
    np.savetxt('pca_aa_original_%d_v_offset_%d.txt'%(n_components,offset),v)
    np.savetxt('pca_aa_original_%d_ev_offset_%d.txt'%(n_components,offset),ev)

    for directory_perturb, filename in zip(directories_perturb,filenames):

        positions = []
        for i_frame in range(n_frame):
            pdbfile = directory + directory_perturb + '/pdb/offset_%d_perturbed/'%(offset+50) + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            #CA = u.select_atoms('name CA')
            CA = u.select_atoms('protein')
            CA_pos = CA.positions
            CA_pos -= CA.center_of_mass()
            positions.append(CA_pos.flatten())
        positions = np.array(positions)
        pca = PCA(n_components=n_components)
        pca.fit(positions)
        w = pca.singular_values_
        v = pca.components_
        ev = pca.explained_variance_
        np.savetxt('pca_aa_%s_%d_w_offset_%d.txt'%(filename,n_components,offset+50),w)
        np.savetxt('pca_aa_%s_%d_v_offset_%d.txt'%(filename,n_components,offset+50),v)
        np.savetxt('pca_aa_%s_%d_ev_offset_%d.txt'%(filename,n_components,offset+50),ev)
