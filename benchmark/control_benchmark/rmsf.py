#!/bin/python3

import MDAnalysis as mda
import numpy as np
import os

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

    for filename, directory in zip(filenames,directories_perturb):
        """ 
        positions = []
        for i_frame in range(n_frame):
            pdbfile = directory + '/pdb/WT_offset_%d/'%offset + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            CA_pos = u.select_atoms('name CA').positions
            positions.append(CA_pos)

        positions = np.array(positions)
        avg_positions = np.average(positions,axis=0)
        deviations = positions - avg_positions
        sf = np.sum(np.square(deviations),axis=2)
        rmsf = np.sqrt(np.average(sf,axis=0))
        np.savetxt('rmsf_original_offset_%d.txt'%offset,rmsf)
        """
        positions = []
        for i_frame in range(n_frame):
            pdbfile = directory + '/pdb/offset_%d_perturbed/'%(offset+50) + 'WT_frame%d.pdb'%i_frame
            u = mda.Universe(pdbfile)
            CA_pos = u.select_atoms('name CA').positions
            positions.append(CA_pos)
        positions = np.array(positions)
        avg_positions = np.average(positions,axis=0)
        deviations = positions - avg_positions
        sf = np.sum(np.square(deviations),axis=2)
        rmsf = np.sqrt(np.average(sf,axis=0))
        np.savetxt('rmsf_%s_offset_%d.txt'%(filename,offset+50),rmsf)

"""
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

rmsf_o = np.loadtxt('rmsf_original.txt')

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 2

colors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']

fig, ax = plt.subplots(figsize=(8,6))

x = np.arange(65,65+rmsf_o.size,1)

for i, filename in enumerate(filenames):
    rmsf = np.loadtxt('rmsf_%s.txt'%filename)
    drmsf = rmsf - rmsf_o
    ax.plot(x,drmsf,color=colors[i],label=labels[i])

#ax.set_xlim()
#ax.set_ylim(0,4e-5)

ax.set_xlabel('Residue #',fontsize=30)
ax.set_ylabel('$\Delta$RMSF ($\mathrm{\AA}$)',fontsize=30)
ax.yaxis.get_offset_text().set_fontsize(24)
#yml = MultipleLocator(0.00001)
#ax.yaxis.set_minor_locator(yml)
#ax.ticklabel_format(axis='y',style='sci',scilimits=(0.1,1.0))

ax.tick_params(axis='y',labelsize=24)
ax.tick_params(axis='x',labelsize=24)
ax.tick_params(axis='both',which='major',width=2,length=10)
ax.tick_params(axis='both',which='minor',width=2,length=5)

#x = np.linspace(0,1,10)
#ax.plot(x,x,color='k')
plt.legend(fontsize=18,loc='upper left',
#        bbox_to_anchor=(0.15,1),
        frameon=False)
plt.tight_layout()
plt.savefig('rmsf_perturb.jpg',dpi=300)
plt.show()
"""
