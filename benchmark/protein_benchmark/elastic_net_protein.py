#!/bin/python3

import numpy as np
import MDAnalysis as mda
from sklearn.linear_model import SGDClassifier
import os
from matplotlib import pyplot as plt

n_frame = 100

protA = 'WT'
protB = 'R164S'
directories = [
               'WT_R164S_whole',
               'WT_R164S_65_230',
               'WT_R164S_65_213',
              ]

#protA = 'WT'
#protB = 'I50V'
#directories = [
#                 'hiv_A',
#                 'hiv_B',
#              ]

#protA = 'WT'
#protB = 'M290A'
#directories = [
#                 'abl1_WT_M290A_227_512',
#                 'abl1_WT_M290A_242_502',
#                 'abl1_WT_M290A_242_315',
#              ]

#protA = 'gtp'
#protB = 'gdp'
#directories = [
#                 '1ttt_gtp_gdp_whole',
#                 '1ttt_gtp_gdp_208_308',
#                 '1ttt_gtp_gdp_311_405',
#              ]

#protA = 'ibb'
#protB = 'apo'
#directories = ['imp_ibb_apo']

#selection = 'protein and not type H'
selection = 'protein and not (resid 164 and not backbone) and not type H'

filename = directory #+ '_noh'

for offset in range(0,80,8):
    print(filename, offset)
    positions = []
    y = []
    for i_frame in range(n_frame):
        pdbfile = '../' + directory + '/pdb/%s_offset_%d/'%(protA,offset) + '%s_frame%d.pdb'%(protA,i_frame)
        u = mda.Universe(pdbfile)
        protein = u.select_atoms(selection)
        protein_pos = protein.positions
        positions.append(protein_pos.flatten())
        y.append(0)
    for i_frame in range(n_frame):
        pdbfile = '../' + directory + '/pdb/%s_offset_%d/'%(protB,offset) + '%s_frame%d.pdb'%(protB,i_frame)
        u = mda.Universe(pdbfile)
        protein = u.select_atoms(selection)
        protein_pos = protein.positions
        positions.append(protein_pos.flatten())
        y.append(1)
    positions = np.array(positions) 
    en = SGDClassifier(loss="log", penalty="elasticnet")
    en.fit(positions,y)
    np.savetxt('elasticnetlog_coef_%s_offset_%d.txt'%(filename,offset),en.coef_[0])

