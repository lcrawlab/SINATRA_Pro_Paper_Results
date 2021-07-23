#!/bin/python3

for offset in 0 10 20 30 40 50 60 70 80 90
do
python3 main_protein.py --protA WT --protB R164S \
                --directory "WT_R164S_65_213_6.0" \
                --n_sample 100 \
                --offset ${offset} \
                --struct_file_A "data/WT.gro" \
                --traj_file_A "data/WT.xtc" \
                --struct_file_B "data/R164S.gro" \
                --traj_file_B "data/R164S.xtc" \
                --selection "protein and resid 65:213" \
                --radius 6.0 \
                --n_cone 20 \
                --n_direction_per_cone 4 \
                --cap_radius 0.80 \
                --ec_type "DECT" \
                --n_filtration 120 \
                --n_mcmc 100000 \
                --parallel --n_core 16 --verbose
done
