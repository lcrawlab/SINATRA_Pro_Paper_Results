#!/bin/python3

import numpy as np, os, sys
import MDAnalysis as mda
from sklearn.neighbors import NearestNeighbors
from scipy import stats

method = "DECT"
protA = "WT"
protB = "R164S"
c = 20
d = 8
theta = 0.80
l = 120

offsets = range(0,100,10)

datasets = [
            "WT_R164S_whole",
            "WT_R164S_65_230",
            "WT_R164S_65_213",
           ]

selections = [
              'resid 163:178',
              'resid 210:230',
             ]

n_offset = len(offsets)
for sm_radius in [2.0,4.0,6.0]:
    for selection in selections:
        for dataset in datasets: 
            u = mda.Universe('../%s/pdb/%s_offset_0/%s_frame0.pdb'%(dataset,protA,protA))
            protein = u.select_atoms('protein')
            roi = u.select_atoms('protein and (%s)'%selection)
            if len(roi) == 0:
                continue
            n_atom = len(protein)
            volume = len(roi)
            nbrs = NearestNeighbors(n_neighbors=volume, algorithm='ball_tree').fit(protein.positions)
            distances, indices = nbrs.kneighbors(protein.positions)

            rate_vert = np.loadtxt('../%s_%.1f_new/vert_prob_%s_%s_%s_%d_%d_%.2f_%d.txt'%(dataset,sm_radius,method,protA,protB,c,d,theta,l))
            roi_rate = np.sum(rate_vert[roi.ix])
            nn_rates = np.zeros(n_atom,dtype=float)
            n_nreg = 0
            for i in range(n_atom):
                if any(i in roi.ix for i in indices[i]):
                    nn_rates[i] = 0
                else:
                    nn_rates[i] = np.sum(rate_vert[indices[i]])
                    n_nreg += 1
            p_value = (np.sum(nn_rates > roi_rate)+1)/(n_nreg+1)
            bayes_factor = -1.0/(np.exp(1)*p_value*np.log(p_value))
            print(dataset,sm_radius,selection,p_value,bayes_factor)

