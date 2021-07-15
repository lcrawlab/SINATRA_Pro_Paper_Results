#!/bin/python3

import os
from traj_reader import *
from control_simulation import *
import multiprocessing
from mesh import *

n_core = multiprocessing.cpu_count()
print("Detected %d cores"%n_core)

n_core = 4

n_sample = 100
prot = "WT"
protA = "original"
protB = "perturbed"

## Simplicies construction parameters
selection = 'resid 65:213'
sm_radius = 1.0

## perturb parameter
region_selection = 'resid 163:178'
r = 0.5

## Input files
struct_file = 'data/WT/md_0_1.gro'
traj_file = 'data/WT/md_0_1_noPBC.xtc'

directory = "simulation_perturb_omega_loop_sphere_%.1f"%r # output folder

if not os.path.exists(directory):
    os.mkdir(directory)
if not os.path.exists(directory + '/pdb'):
    os.mkdir(directory + '/pdb')

for offset in range(0,50,10):
    ## Original set of frames
    convert_traj_pdb_aligned(protA = prot, protB = prot, struct_file_A = struct_file, traj_file_A = traj_file, struct_file_B = struct_file, traj_file_B = traj_file, align_frame = offset, n_sample = n_sample, offset = offset, selection = selection, directory = directory, verbose = True)
    ## Pick different set of frames from trajectory
    convert_traj_pdb_aligned(protA = prot, protB = prot, struct_file_A = struct_file, traj_file_A = traj_file, struct_file_B = struct_file, traj_file_B = traj_file, align_frame = offset, n_sample = n_sample, offset = offset+50, selection = selection, directory = directory, verbose = True)   
    ## Perturb pdb positions
    perturb_protein_region_sphere(selection=region_selection, r = r, prot = prot, n_sample = n_sample, directory_original = "%s/pdb/WT_offset_%d"%(directory,offset+50), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+50), directory_new = directory, parallel = True, n_core = n_core) 
    ## Convert PDB to mesh
    convert_pdb_mesh("original", "perturbed", n_sample = n_sample, sm_radius = sm_radius, directory_pdb_A = "%s/pdb/WT_offset_%d"%(directory,offset), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+50), directory_mesh = "%s/msh_offset_%d/"%(directory,offset), parallel = True, n_core = n_core)

