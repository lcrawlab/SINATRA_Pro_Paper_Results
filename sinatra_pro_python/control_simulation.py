#!/bin/python3

import os, sys
import numpy as np
import MDAnalysis as mda
from traj_reader import *
import multiprocessing
from joblib import Parallel, delayed

def calc_radius_write_pdb_perturb_static(directory_original, prot, selection, directory_pdb_B, i_sample, dr_vector = [0.5,0.5,0.5]):
    original_pdb_file = '%s/%s_frame%d.pdb'%(directory_original,prot,i_sample)
    u = mda.Universe(original_pdb_file).select_atoms('protein')
    original_pos = u.atoms.positions
    r1 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    u_new = mda.Universe(original_pdb_file)
    region = u_new.select_atoms(selection).ix
    for i in range(len(u.atoms)):
        if u_new.atoms[i].ix in region:
            u_new.atoms[i].position += dr_vector
    u_new.atoms.write('%s/%s_frame%d.pdb'%(directory_pdb_B,prot,i_sample))
    r2 = np.amax(np.linalg.norm(u_new.atoms.positions,axis=1))
    return [r1,r2]

def perturb_protein_region_static(selection='resid 163:178', dr_vector = [0.5,0.5,0.5], prot = "WT", n_sample = 101, directory_original = 'WT_R164S/pdb/WT/', directory_new = "simulation", directory_pdb_B = None, parallel = False, n_core = -1):   
    if parallel:
        if n_core == -1:
            n_core = multiprocessing.cpu_count()   
    if not os.path.exists(directory_new):
        os.mkdir(directory_new)    
    if not os.path.exists('%s/pdb'%directory_new):
        os.mkdir('%s/pdb'%directory_new)
    if directory_pdb_B == None:
        directory_pdb_B = '%s/pdb/perturbed/'%directory_new
    if not os.path.exists(directory_pdb_B):
        os.mkdir(directory_pdb_B)
    if parallel:
        processed_list = Parallel(n_jobs=n_core)(delayed(calc_radius_write_pdb_perturb_static)(directory_original, prot, selection, directory_pdb_B, i_sample, dr_vector) for i_sample in range(n_sample))
        r = np.array(processed_list).flatten()
    else:
        r = np.array([])
        for i_sample in range(n_sample):
            output = calc_radius_write_pdb_perturb_static(directory_original, prot, selection, directory_pdb_B, i_sample, dr_vector)
            r = np.append(r,output)
    rmax = np.amax(r) 
    return rmax

def calc_radius_write_pdb_perturb_sphere(directory_original, prot, selection, directory_pdb_B, i_sample, r = 1.0):
    original_pdb_file = '%s/%s_frame%d.pdb'%(directory_original,prot,i_sample)
    u = mda.Universe(original_pdb_file).select_atoms('protein')
    original_pos = u.atoms.positions
    r1 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    u_new = mda.Universe(original_pdb_file)
    region = u_new.select_atoms(selection).ix 
    x = np.random.normal(loc=0.0,scale=1.0,size=len(region))
    y = np.random.normal(loc=0.0,scale=1.0,size=len(region))
    z = np.random.normal(loc=0.0,scale=1.0,size=len(region))
    v = np.vstack((x,y,z))
    v = v/np.linalg.norm(v,axis=0)*r
    v = v.T
    k = 0
    for i in range(len(u.atoms)):
        if u_new.atoms[i].ix in region:
            u_new.atoms[i].position += v[k]
            k += 1
    u_new.atoms.write('%s/%s_frame%d.pdb'%(directory_pdb_B,prot,i_sample))
    r2 = np.amax(np.linalg.norm(u_new.atoms.positions,axis=1))
    return [r1,r2]

def perturb_protein_region_sphere(selection='resid 163:178', r = 1.0, prot = "WT", n_sample = 101, directory_original = 'WT_R164S/pdb/WT/', directory_new = "simulation", directory_pdb_B = None, parallel = False, n_core = -1):    
    if parallel:
        if n_core == -1:
            n_core = multiprocessing.cpu_count()
    if not os.path.exists(directory_new):
        os.mkdir(directory_new)
    if not os.path.exists('%s/pdb'%directory_new):
        os.mkdir('%s/pdb'%directory_new)
    if directory_pdb_B == None:
        directory_pdb_B = '%s/pdb/perturbed/'%directory_new
    if not os.path.exists(directory_pdb_B):
        os.mkdir(directory_pdb_B)
    if parallel:
        processed_list = Parallel(n_jobs=n_core)(delayed(calc_radius_write_pdb_perturb_sphere)(directory_original, prot, selection, directory_pdb_B, i_sample, r) for i_sample in range(n_sample))
        r = np.array(processed_list).flatten()
    else:
        r = np.array([])
        for i_sample in range(n_sample):
            output = calc_radius_write_pdb_perturb_new(directory_original, prot, selection, directory_pdb_B, i_sample, r)
            r = np.append(r,output)
    rmax = np.amax(r)
    return rmax

