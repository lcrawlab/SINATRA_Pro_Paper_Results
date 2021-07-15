#!/bin/python3

import os
from traj_reader import *
from control_simulation import *

import multiprocessing

n_core = multiprocessing.cpu_count()
print("Detected %d cores"%n_core)

prot = "WT"
protA = "original"
protB = "perturbed"

## Simplicies construction parameters
selection = 'resid 65:213'
sm_radius = 1.0

## perturb parameter
region_selection = 'resid 163:178'
dr_vector = [0.5,0.5,0.5] # displacement vector

## Input files
struct_file = 'data/WT/md_0_1.gro'
traj_file = 'data/WT/md_0_1_noPBC.xtc'

directory = "simulation_perturb_omega_loop_static_%.1f_%.1f_%.1f"%(dr_vector[0],dr_vector[1],dr_vector[2]) # output folder

if not os.path.exists(directory):
    os.mkdir(directory)
if not os.path.exists(directory + '/pdb'):
    os.mkdir(directory + '/pdb')
if not os.path.exists(directory + '/msh'):
    os.mkdir(directory + '/msh')

for offset in range(0,50,10):
    ## Original set of frames
    convert_traj_pdb_aligned(protA = prot, protB = prot, struct_file_A = struct_file, traj_file_A = traj_file, struct_file_B = struct_file, traj_file_B = traj_file, align_frame = offset, n_sample = n_sample, offset = offset, selection = selection, directory = directory) 
    ## Pick different set of frames from trajectory
    convert_traj_pdb_aligned(protA = prot, protB = prot, struct_file_A = struct_file, traj_file_A = traj_file, struct_file_B = struct_file, traj_file_B = traj_file, align_frame = offset, n_sample = n_sample, offset = offset+50, selection = selection, directory = directory)   
    ## Perturb pdb positions
    perturb_protein_region_static(selection=region_selection, dr_vector = dr_vector, prot = prot, n_sample = n_sample, directory_original = "%s/pdb/WT_offset_%d"%(directory,offset+50), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+50), directory_new = directory, parallel = True, n_core = n_core) 
    ## Convert PDB to mesh
    convert_pdb_mesh("original", "perturbed", n_sample = n_sample, sm_radius=sm_radius, directory_pdb_A = "%s/pdb/WT_offset_%d"%(directory,offset), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+50), directory_mesh = "%s/msh_offset_%d/"%(directory,offset), parallel = True, n_core = n_core)


