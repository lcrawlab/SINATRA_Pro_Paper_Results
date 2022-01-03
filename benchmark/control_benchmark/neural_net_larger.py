#!/bin/python3

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import backend as k
import MDAnalysis as mda
import numpy as np

#tf.compat.v1.disable_eager_execution()

directory_original = '../../sinatra_pro/WT_R164S_65_213/'

directory = '../../sinatra_pro/'

directories_perturb = ['simulation_perturb_omega_loop_sphere_0.5_0.0_offset/',
                       'simulation_perturb_omega_loop_sphere_1.0_0.0_offset/',
                       'simulation_perturb_omega_loop_sphere_2.0_0.0_offset/',
                       'simulation_perturb_omega_loop_vector_0.0_0.5_0.5_0.5_offset/',
                       'simulation_perturb_omega_loop_vector_0.0_1.0_1.0_1.0_offset/',
                       'simulation_perturb_omega_loop_vector_0.0_2.0_2.0_2.0_offset/']
                        

filenames = ['perturb_0.5_0.0_sphere_offset',
             'perturb_1.0_0.0_sphere_offset',
             'perturb_2.0_0.0_sphere_offset',
             'perturb_0.0_0.5_vector_offset',
             'perturb_0.0_1.0_vector_offset',
             'perturb_0.0_2.0_vector_offset']

for directory_perturb, filename in zip(directories_perturb,filenames):

    for offset in range(0,50,10):
 
        directory_A = directory + directory_perturb + 'pdb/WT_offset_%d/'%offset
        directory_B = directory + directory_perturb + 'pdb/offset_%d_perturbed/'%(offset+50)

        x = []
        y = []
        for frame in range(100):
            pdbfile = directory_A + 'WT_frame%d.pdb'%frame
            u = mda.Universe(pdbfile)
            positions = u.select_atoms('protein').positions
            x.append(positions.flatten())
            y.append(0)
        for frame in range(100):
            pdbfile = directory_B + 'WT_frame%d.pdb'%frame
            u = mda.Universe(pdbfile)
            positions = u.select_atoms('protein').positions
            x.append(positions.flatten())
            y.append(1)
        x = np.array(x)
        y = np.array(y)

        n_sample = x.shape[0]
        dim = x.shape[1]
        mu = 0.01

        model = keras.Sequential(
            [
             keras.Input(shape=dim),
             layers.BatchNormalization(),
             layers.Dense(2048, activation="relu", kernel_regularizer=keras.regularizers.l2()),
             layers.BatchNormalization(),
             layers.Dense(2048, activation="relu", kernel_regularizer=keras.regularizers.l2()),
             layers.BatchNormalization(),
             layers.Dense(512, activation="relu", kernel_regularizer=keras.regularizers.l2()),
             layers.BatchNormalization(),
             layers.Dense(128, activation="relu", kernel_regularizer=keras.regularizers.l2()),
             layers.BatchNormalization(),
             layers.Dense(1, name = "last_block"),
             layers.Dense(1, activation="sigmoid", use_bias=False)
            ]
        )
        optimizer = keras.optimizers.Adam(learning_rate=0.001)
        model.compile(loss = "binary_crossentropy",
                      optimizer = optimizer,
                      metrics = ["accuracy"])

        model.fit(x,y,
                  batch_size = 50,
                  epochs = 400)

        with tf.GradientTape() as tape:
            x = tf.convert_to_tensor(x)
            tape.watch(x)
            predictions = model(x)
            loss = predictions
            gradient = tape.gradient(loss,x) #preds,model.trainable_weights)
            gradient = tf.reduce_max(gradient, axis=0)
            gradient = gradient.numpy()
            min_val, max_val = np.amin(gradient), np.amax(gradient)
            saliency_map = (gradient - min_val) / (max_val - min_val + k.epsilon())
            np.savetxt('larger_reg_L2/saliency_map_%s_offset_%d.txt'%(filename,offset),saliency_map)
