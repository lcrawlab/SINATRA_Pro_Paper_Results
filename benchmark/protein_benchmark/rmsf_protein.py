#!/bin/python3

import MDAnalysis as mda
import numpy as np
import os, sys

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

n_frame = 1000

for directory in directories:
     
    positions = []
    CA_res = []
    for i_frame in range(n_frame):
        pdbfile = directory + '/pdb/' + protA + '/' + '%s_frame%d.pdb'%(protA,i_frame)
        u = mda.Universe(pdbfile)
        CA = u.select_atoms('name CA')
        CA_pos = CA.positions
        if len(CA_res) == 0:
            CA_res = CA.resids
        else:
            if np.any(np.not_equal(CA_res,CA.resids)):
                print("residue not matching!")
                exit()
        positions.append(CA_pos)
        sys.stdout.write("%s frame %d...\r"%(protA,i_frame))
        sys.stdout.flush()
    positions = np.array(positions)
    avg_positions = np.average(positions,axis=0)
    deviations = positions - avg_positions
    sf = np.sum(np.square(deviations),axis=2)
    rmsfA = np.sqrt(np.average(sf,axis=0))

    positions = []
    for i_frame in range(n_frame):
        pdbfile = directory + '/pdb/' + protB + '/' + '%s_frame%d.pdb'%(protB,i_frame)
        u = mda.Universe(pdbfile)
        CA = u.select_atoms('name CA')
        CA_pos = CA.positions
        positions.append(CA_pos)
        if len(CA_res) == 0:
            CA_res = CA.resids
        else:
            if np.any(np.not_equal(CA_res,CA.resids)):
                print("residue not matching!")
                exit()    
        sys.stdout.write("%s frame %d...\r"%(protA,i_frame))
        sys.stdout.flush()

    positions = np.array(positions)
    avg_positions = np.average(positions,axis=0)
    deviations = positions - avg_positions
    sf = np.sum(np.square(deviations),axis=2)
    rmsfB = np.sqrt(np.average(sf,axis=0))
    
    drmsf = rmsfB - rmsfA
    
    np.savetxt('rmsf_%s.txt'%directory,np.c_[CA_res,rmsfA,rmsfB,drmsf])
    
