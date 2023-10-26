"""
Skript for Damask3 Postprocessing, including macro - micro strain and defgrad comparison for RVE with phase field Damage
For 1 - 3 Phases

Niklas Fehlemann, IMS
"""

import damask
import os
import glob
import sys
import numpy as np
import pandas as pd
# import seaborn as sns
# import pyvista as pv
import yaml
import matplotlib.pyplot as plt

if __name__ == "__main__":

    result_file = '/home/p0021070/damask/test/grid_load.hdf5'
    images = '/home/p0021070/damask/test/evaluation_images'

    result = damask.Result(result_file)
    # result.add_strain(F='F')
    # result.add_stress_Cauchy(F='F')
    # result.add_equivalent_Mises('sigma')
    # result.add_equivalent_Mises('epsilon_V^0.0(F)')


    # Averaged Quantities
    strain_eq = result.place('epsilon_V^0.0(F)')
    stress_eq = result.place('sigma')

    strain_vector = list()
    stress_vector = list()

    for _, strain in strain_eq.items():
        strain_mean = np.abs(strain[:, 2, 2].mean())
        strain_vector.append(strain_mean)

    for _, stress in stress_eq.items():
        stress_mean = np.abs(stress[:, 2, 2].mean()/1e6)
        stress_vector.append(stress_mean)


    plt.plot(strain_vector, stress_vector, label='True Stress')
    plt.legend()
    plt.xlabel("True Strain")
    plt.ylabel("True Stress (MPa)")
    plt.savefig(f"{images}/test_result.png")

    # # CRSS Distribution
    # crss_start= np.mean(result.view('increment', 2).get('xi_sl'))
    # crss_96 = np.mean(result.view('increment', 100).get('xi_sl'))
    # crss_146 = np.mean(result.view('increment', 154).get('xi_sl'))
    # crss_216 = np.mean(result.view('increment', -1).get('xi_sl'))
    # print(crss_start)

    # # Partitioned Flowcurves
    # ferrite_stress = list()
    # martensite_stress = list()
    # for key, val in result.get('sigma').items():
    #     temp = list()
    #     for subkey, subvalue in val.items():
    #         if 'Ferrite' in subkey:
    #             temp.append(np.mean(subvalue[:,0,0]))
    #         else:
    #             martensite_stress.append(np.mean(subvalue[:,0,0]))
    #     ferrite_stress.append(np.mean(temp))

    # ferrite_strain = list()
    # martensite_strain = list()
    # for key, val in result.get('epsilon_V^0.0(F)').items():
    #     temp = list()
    #     for subkey, subvalue in val.items():
    #         if 'Ferrite' in subkey:
    #             temp.append(np.mean(subvalue[:,0,0]))
    #         else:
    #             martensite_strain.append(np.mean(subvalue[:,0,0]))
    #     ferrite_strain.append(np.mean(temp))

    # fig, ax = plt.subplots(1, 1, figsize=(12,8))

    # ax.plot(ferrite_strain, ferrite_stress, label='Avg. Stress Ferrite', linewidth=4)
    # ax.plot(martensite_strain, martensite_stress, label='Avg. Stress Martensite', linewidth=4)
    # ax.legend(fontsize=24)
    # ax.set_xlabel('Strain', fontsize=22)
    # ax.set_ylabel('Stress', fontsize=22)
    # ax.tick_params(which='both', size=10, labelsize=20)

    # fig.tight_layout()
    # fig.savefig(path + r'/FlowcurvePartitioning.png')
    # plt.close()


