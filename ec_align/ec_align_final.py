#!/bin/python3

import os, sys

import mesh
import euler
import directions
import align

from copy import deepcopy

import numpy as np
from scipy.spatial.transform import Rotation as R
import MDAnalysis as mda
from scipy import signal
from sklearn.linear_model import RANSACRegressor
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.stats import wasserstein_distance

import multiprocessing
from joblib import Parallel, delayed

from matplotlib import pyplot as plt

n_core = 48  # number of CPU cores used in parallelized calculations

"""
sm_radius = 1.0  # radius cutoff in simplicial complex construction

### Fine alignment parameters ###
n_iter_fine = 1  # number of iterations for fine alignment
nd = 100  # number of directions to calculate DECT
w = 0.5  # Bandwidth of Gaussian kernel applied on DEC for smoothing
k_nf = 200  # number of filtraion per Angstrom in Radial DECT
ransac_n = 100  # number of draws in each iteration of RANSAC
ransac_tol = 0.9  # tolerance in dot product between aligned directions in RANSAC
ransac_n_iter = 1000  # number of iteraction in RANSAC

"""

offset = int(sys.argv[1])
sm_radius = 2.0  # radius cutoff in simplicial complex construction

### Fine alignment parameters ###
n_iter_fine = 3  # number of iterations for fine alignment
nd = 500  # number of directions to calculate DECT
w = 1.0  # Bandwidth of Gaussian kernel applied on DEC for smoothing
k_nf = 200  # number of filtraion per Angstrom in Radial DECT
ransac_n = 100  # number of draws in each iteration of RANSAC
ransac_tol = 0.9  # tolerance in dot product between aligned directions in RANSAC
ransac_n_iter = 20000  # number of iteraction in RANSAC
sigma = int(np.ceil(w*k_nf))


def convolve_single(y, kernel, kernelsum):
    y_new = signal.convolve(y, kernel, mode='same') / kernelsum
    return y_new

def convolve_parallel(ec, kernel, kernelsum, n_core):
    processed_list = Parallel(n_jobs=n_core)(delayed(convolve_single)(
        ec[idir], kernel, kernelsum) for idir in range(nd))
    return np.array(processed_list)

def ECalign(protA, protB, meshA, meshB, directions_set, n_iter_fine = 3,
    nd = 500, w = 0.5, k_nf = 200,
    ransac_n = 20, ransac_tol = 0.9, ransac_n_iter = 20000):

    """
    Alignment using directional DECT
    """

    sigma = int(np.ceil(w*k_nf))
    rot_mat = np.identity(3, dtype=float)
    
    radiusA = meshA.calc_radius()
    radiusB = meshB.calc_radius()
    if radiusA > radiusB:
        rmax = radiusA
    else:
        rmax = radiusB
    
    nf = int(np.ceil(rmax*k_nf))
    sigma = int(np.ceil(w*k_nf))
    kernel = signal.windows.gaussian(nf, sigma, sym=True)
    kernelsum = np.sum(kernel)
    
    ecA_k = np.zeros((nd, nf), dtype=float)
    ecB_k = np.zeros((nd, nf), dtype=float)

    for i in range(n_iter_fine):


        rot_mat = np.identity(3, dtype=float)
        if i == 0:
            r, ecA = euler.compute_ec_curve_parallel(
                meshA, directions_set, ball_radius=rmax, n_filtration=nf,
                ec_type="DECT", include_faces=False, n_core=n_core)
            ecA_k = convolve_parallel(ecA, kernel, kernelsum, n_core)
        r, ecB = euler.compute_ec_curve_parallel(
            meshB, directions_set, ball_radius=rmax, n_filtration=nf,
            ec_type="DECT", include_faces=False, n_core=n_core)
        ecB_k = convolve_parallel(ecB, kernel, kernelsum, n_core)

        dev = cdist(ecA_k, ecB_k, metric='sqeuclidean')
        row, col = linear_sum_assignment(dev)

        a = directions_set[row]
        b = directions_set[col]
        #match_dev = dev[row, col]

        rot, rmsd = align.ransac_directions(a, b, n=ransac_n, tol=ransac_tol, n_iter=ransac_n_iter, n_core=n_core)
        print(rot.as_euler('zxy', degrees=True))
        meshB.vertices = rot.apply(meshB.vertices)
        protB.atoms.rotate(rot.as_matrix())

    return


def ECalign_to_mesh(protA, protB, struct_file_A, traj_file_A, struct_file_B, traj_file_B, align_frame = 0,
    n_sample = 100, selection = None, directory = None, offset = 0, single = False, verbose = False,
    n_iter_fine = 3,
    nd = 500, w = 0.5, k_nf = 200,
    ransac_n = 20, ransac_tol = 0.9, ransac_n_iter = 10000):

    if not os.path.exists(directory):
        os.mkdir(directory)
    if not os.path.exists(directory+"/pdb"):
        os.mkdir(directory+"/pdb")

    directions_set = directions.generate_equidistributed_cones(
        n_cone=nd, n_direction_per_cone=1, hemisphere=False)

    radii = []
    for struct_file, traj_file in zip([struct_file_A,struct_file_B],[traj_file_A,traj_file_B]):
        u = mda.Universe(struct_file_A,traj_file_A)
        n_frame = len(u.trajectory)
        nskip =  int(n_frame/n_sample)
        frame = 0
        for ts in u.trajectory:
            if (frame-offset) % nskip == 0:
                prot = u.select_atoms(selection)
                prot.translate(prot.center_of_mass())
                radius = np.amax(np.linalg.norm(prot.atoms.positions,axis=1))
                radii.append(radius)
            frame += 1
    rmax = np.amax(radii)            

    ref = mda.Universe(struct_file_A,traj_file_A)
    ref.trajectory[align_frame]
    ref_atoms = ref.select_atoms(selection)
    ref_atoms.translate(-ref_atoms.select_atoms("backbone").center_of_mass())
    
    meshref = mesh.mesh()
    meshref.vertices = deepcopy(ref_atoms.select_atoms("backbone").positions)
    meshref.convert_vertices_to_mesh(sm_radius=sm_radius,msh_file='%s/ref_offset_%d.msh'%(directory,offset), rmax=rmax, make_faces=True)

    for prot, struct_file, traj_file in zip([protA,protB],[struct_file_A,struct_file_B],[traj_file_A,traj_file_B]):

        """
        File parsing and simplicial complex construction
        """
        if not single:
            directory_pdb = "%s/pdb/%s_offset_%d"%(directory,prot,offset)
            directory_mesh = "%s/msh_%s_%.1f_offset_%d/"%(directory,prot,sm_radius,offset)
        else:
            directory_pdb = "%s/pdb/%s"%(directory,prot)
            directory_mesh = "%s/msh_%s_%.1f/"%(directory,offset,prot,sm_radius)

        if not os.path.exists(directory_pdb):
            os.mkdir(directory_pdb)
        if not os.path.exists(directory_pdb+'/unaligned'):
            os.mkdir(directory_pdb+'/unaligned')
        if not os.path.exists(directory_pdb+'/aligned'):
            os.mkdir(directory_pdb+'/aligned')
        if not os.path.exists(directory_mesh):
            os.mkdir(directory_mesh)
        if not os.path.exists(directory_mesh+'/unaligned'):
            os.mkdir(directory_mesh+'/unaligned')
        if not os.path.exists(directory_mesh+'/aligned'):
            os.mkdir(directory_mesh+'/aligned')

        u =  mda.Universe(struct_file, traj_file)
        n_frame = len(u.trajectory)
        nskip =  int(n_frame/n_sample)
    
        frame = 0
        i_sample = 0
        for ts in u.trajectory:
            if (frame-offset) % nskip == 0:
                if verbose:
                    print("%s %d"%(prot, i_sample))
                if os.path.exists('%s/aligned/%s_frame%d.msh'%(directory_mesh,prot,i_sample)):
                     frame += 1
                     i_sample += 1
                     continue
                u_atoms = u.select_atoms(selection)
                u_atoms.translate(-u_atoms.select_atoms("backbone").center_of_mass())
                u_atoms.atoms.write('%s/unaligned/%s_frame%d.pdb'%(directory_pdb,prot,i_sample))
                meshu = mesh.mesh()
                meshu.vertices = deepcopy(u_atoms.select_atoms("backbone").positions)
                meshu.convert_vertices_to_mesh(sm_radius=sm_radius,msh_file='%s/unaligned/%s_frame%d.msh'%(directory_mesh,prot,i_sample), rmax=rmax, make_faces=True)
                ECalign(ref, u_atoms, meshref, meshu, directions_set,
                    n_iter_fine, nd, w, k_nf, ransac_n, ransac_tol, ransac_n_iter)
                meshu.write_mesh_file(filename='%s/aligned/%s_frame%d.msh'%(directory_mesh,prot,i_sample))
                u_atoms.atoms.write('%s/aligned/%s_frame%d.pdb'%(directory_pdb,prot,i_sample))
                i_sample += 1
                if i_sample == n_sample:
                    break
            frame += 1

ECalign_to_mesh("WT", "R164S",
    "../../../sinatra_pro/data/WT/md_0_1.gro",
    "../../../sinatra_pro/data/WT/md_0_1_noPBC.xtc",
    "../../../sinatra_pro/data/R164S/md_0_1.gro",
    "../../../sinatra_pro/data/R164S/md_0_1_noPBC.xtc",
    align_frame = offset,
    n_sample = 100,
    selection = "protein and resid 65:230",
    directory = "../WT_R164S_65_230_bbalign_new/",
    offset = offset,
    single = False,
    verbose = True,
    n_iter_fine = n_iter_fine,
    nd = nd, w = w, k_nf = k_nf,
    ransac_n = ransac_n, ransac_tol = ransac_tol, ransac_n_iter = ransac_n_iter
    )
