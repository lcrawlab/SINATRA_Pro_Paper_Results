#!/bin/python3

import MDAnalysis as mda
import numpy as np
import os
from sklearn.decomposition import PCA

directory_original = '../WT_R164S_65_213/'
n_frame = 100

directory = "../../sinatra_pro/"

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

for offset in range(0,50,10):

    positions_A = []
    for i_frame in range(n_frame):
        pdbfile = directory + '/' + directories_perturb[0] + '/pdb/WT_offset_%d/'%offset + 'WT_frame%d.pdb'%i_frame
        u = mda.Universe(pdbfile)
        #CA = u.select_atoms('name CA')
        CA = u.select_atoms('protein')
        CA_pos = CA.positions
        CA_pos -= CA.center_of_mass()
        positions_A.append(CA_pos.flatten())
    positions_A = np.array(positions_A)
    positions_A -= np.mean(positions_A, axis=0)
    positions_A /= np.std(positions_A, axis=0)

    for directory_perturb, filename in zip(directories_perturb,filenames):

        positions_B = []
        for i_frame in range(n_frame):
            pdbfile = directory + '/' + directory_perturb + '/pdb/offset_%d_perturbed/'%(offset+50) + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            #CA = u.select_atoms('name CA')
            CA = u.select_atoms('protein')
            CA_pos = CA.positions
            CA_pos -= CA.center_of_mass()
            positions_B.append(CA_pos.flatten())
        positions_B = np.array(positions_B)
        positions_B -= np.mean(positions_B, axis=0)
        positions_B /= np.std(positions_B, axis=0)

        cov = np.einsum('ki,kj->ij',positions_A,positions_B) / positions_A.shape[0]
                
        pca = PCA(n_components=100)
        pca.fit(cov)
        w = pca.singular_values_
        v = pca.components_
        ev = pca.explained_variance_
        np.savetxt('pca_aa_%s_%d_w_offset_%d.txt'%(filename,n_components,offset),w)
        np.savetxt('pca_aa_%s_%d_v_offset_%d.txt'%(filename,n_components,offset),v)
        np.savetxt('pca_aa_%s_%d_ev_offset_%d.txt'%(filename,n_components,offset),ev)

