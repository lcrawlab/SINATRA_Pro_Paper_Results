#!/bin/python3

import MDAnalysis as mda
import numpy as np
import os
from sklearn.decomposition import PCA

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
     
    positions = []
    for i_frame in range(n_frame):
        pdbfile = directories_perturb[0] + '/pdb/WT_offset_%d/'%offset + 'WT_frame%d.pdb'%i_frame
        u = mda.Universe(pdbfile)
        CA = u.select_atoms('protein')
        CA_pos = CA.positions
        CA_pos -= CA.center_of_mass()
        positions.append(CA_pos.flatten())
    positions = np.array(positions)
    pca = PCA(n_components=2)
    pca.fit(positions)
    w = pca.singular_values_
    v = pca.components_
    np.savetxt('pca_aa_original_w_offset_%d.txt'%offset,w)
    np.savetxt('pca_aa_original_v_offset_%d.txt'%offset,v)
     
    for directory_perturb, filename in zip(directories_perturb,filenames):
       
        positions = []
        for i_frame in range(n_frame):
            pdbfile = directory_perturb + '/pdb/offset_%d_perturbed/'%(offset+50) + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            CA = u.select_atoms('protein')
            CA_pos = CA.positions    
            CA_pos -= CA.center_of_mass()
            positions.append(CA_pos.flatten())
        positions = np.array(positions)
        pca = PCA(n_components=10)
        pca.fit(positions)
        w = pca.singular_values_
        v = pca.components_
        np.savetxt('pca_aa_%s_w_offset_%d.txt'%(filename,offset+50),w)
        np.savetxt('pca_aa_%s_v_offset_%d.txt'%(filename,offset+50),v)


