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
import seaborn as sns
import pyvista as pv
import yaml
import matplotlib.pyplot as plt

path = os.getcwd()

result = damask.Result(path + r'/grid_load.hdf5')
result.add_strain(F='F')
result.add_stress_Cauchy(F='F')
result.add_equivalent_Mises('sigma')
result.add_equivalent_Mises('epsilon_V^0.0(F)')


# Averaged Quantities
strain_eq = result.place('epsilon_V^0.0(F)')
stress_eq = result.place('sigma')
crss_eq = result.get('xi_sl')

strain_vector = list()
stress_vector = list()
crss_vector = list()

strain_vector = list()
stress_vector = list()
crss_vector = list()
for _, strain in strain_eq.items():  
    strain_mean = np.abs(strain[:, 0, 0].mean())
    strain_vector.append(strain_mean)  

for _, stress in stress_eq.items():
    stress_mean = np.abs(stress[:, 0, 0].mean())
    stress_vector.append(stress_mean)  
    
for key, value in crss_eq.items():
    crss_mean = np.mean(value)
    crss_vector.append(crss_mean)

flowcurve = pd.read_csv(path + r'/TrueStressStrainDP800.csv', header=None)
flowcurve = flowcurve[flowcurve[0] <= 0.05]
real_stress = flowcurve[1] * 1000000
real_strain = flowcurve[0]

plt.style.use('ggplot')
fig, ax = plt.subplots(1, 2, figsize=(20,8))

ax[0].set_title('True Stress and CRSS', fontsize=26)
ax[0].plot(strain_vector, stress_vector, linewidth=3, label='True Stress')
ax[0].scatter(strain_vector, stress_vector)
ax[0].plot(strain_vector, crss_vector, label='CRSS')
ax[0].scatter(strain_vector, crss_vector)
ax[0].set_xlabel('True Strain (-)', fontsize=22)
ax[0].set_ylabel('Stress (MPa)', fontsize=22)
ax[0].set_ylim([0, 10e8])
ax[0].set_xlim([-0.001, 0.16])
ax[0].tick_params(which='both', size=10, labelsize=20)
ax[0].legend(fontsize=26)

ax[1].set_title('Comp with DP800 Flowcurve', fontsize=26)
ax[1].plot(strain_vector, stress_vector, linewidth=3, label='True Stress RVE')
ax[1].scatter(strain_vector, stress_vector)
ax[1].plot(real_strain.to_numpy(), real_stress.to_numpy(), label='DP800 Flowcurve', linewidth=3)
ax[1].set_xlabel('True Strain (-)', fontsize=22)
ax[1].set_ylabel('Stress (MPa)', fontsize=22)
ax[1].set_ylim([0, 10e8])
ax[1].set_xlim([-0.001, 0.16])
ax[1].tick_params(which='both', size=10, labelsize=20)
ax[1].legend(fontsize=26)

fig.tight_layout()
fig.savefig(path + r'/AveragedCurves.png')


# CRSS Distribution
crss_start= np.mean(result.view('increment', 2).get('xi_sl'))
crss_96 = np.mean(result.view('increment', 100).get('xi_sl'))
crss_146 = np.mean(result.view('increment', 154).get('xi_sl'))
crss_216 = np.mean(result.view('increment', -1).get('xi_sl'))
print(crss_start)

"""crss_start = list()
for key, val in crss_dist_start.items():
    crss_start.append(np.mean(val))

crss_96 = list()
for key, val in crss_dist_96.items():
    crss_96.append(np.mean(val))

crss_146 = list()
for key, val in crss_dist_146.items():
    crss_146.append(np.mean(val))

crss_216 = list()
for key, val in crss_dist_216.items():
    crss_216.append(np.mean(val))"""

"""fig, ax = plt.subplots(1,2, figsize=(28,8))

#pd.DataFrame(crss_96).to_csv(path + '/crss_96.csv', index=False)
#pd.DataFrame(crss_146).to_csv(path + '/crss_146.csv', index=False)
#pd.DataFrame(crss_216).to_csv(path + '/crss_216.csv', index=False)

strain_0 = strain_vector[0] * 2.75
strain_96 = strain_vector[50] * 2.75
strain_146 = strain_vector[77] * 2.75
strain_216 = strain_vector[-1] * 2.75
strain_points = [strain_0, strain_96, strain_146, strain_216]
crss_points = [float(np.mean(crss_start)), float(np.mean(crss_96)), float(np.mean(crss_146)), float(np.mean(crss_216))]

mean_hardening = (crss_216 - crss_start) / (strain_216 - strain_0) * 1e-6

ax[0].scatter(strain_points, crss_points)
ax[0].plot(strain_points, crss_points)
ax[0].legend(fontsize=24)
ax[0].set_ylim([80e6, 450e6])
ax[0].set_ylabel('crss / xi_sl', fontsize=22)
ax[0].set_xlabel('Averaged Shear', fontsize=22)
ax[0].set_title(f'Mean hardening is: {mean_hardening:.2f}')
ax[0].tick_params(which='both', size=10, labelsize=20)

strains = [i*2.75 for i in strain_vector]
ax[1].plot(np.asarray(strains), np.asarray(crss_vector), linewidth=4, linestyle='-', label='Mean CRSS')
ax[1].scatter(x=np.asarray(strains), y=np.asarray(crss_vector), s=100, marker='x')
ax[1].legend(fontsize=24)
ax[1].set_ylim([80e6, 450e6])
ax[1].set_ylabel('crss / xi_sl', fontsize=22)
ax[1].set_xlabel('Averaged Shear', fontsize=22)
ax[1].tick_params(which='both', size=10, labelsize=20)

fig.tight_layout()
fig.savefig(path + r'/CRSS_Hardening.png')
"""

# Partitioned Flowcurves
ferrite_stress = list()
martensite_stress = list()
for key, val in result.get('sigma').items():
    temp = list()
    for subkey, subvalue in val.items():
        if 'Ferrite' in subkey:
            temp.append(np.mean(subvalue[:,0,0]))
        else:
            martensite_stress.append(np.mean(subvalue[:,0,0]))
    ferrite_stress.append(np.mean(temp))
    
ferrite_strain = list()
martensite_strain = list()
for key, val in result.get('epsilon_V^0.0(F)').items():
    temp = list()
    for subkey, subvalue in val.items():
        if 'Ferrite' in subkey:
            temp.append(np.mean(subvalue[:,0,0]))
        else:
            martensite_strain.append(np.mean(subvalue[:,0,0]))
    ferrite_strain.append(np.mean(temp))

fig, ax = plt.subplots(1, 1, figsize=(12,8))
    
ax.plot(ferrite_strain, ferrite_stress, label='Avg. Stress Ferrite', linewidth=4)
ax.plot(martensite_strain, martensite_stress, label='Avg. Stress Martensite', linewidth=4)
ax.legend(fontsize=24)
ax.set_xlabel('Strain', fontsize=22)
ax.set_ylabel('Stress', fontsize=22)
ax.tick_params(which='both', size=10, labelsize=20)

fig.tight_layout()
fig.savefig(path + r'/FlowcurvePartitioning.png')
plt.close()


"""# Store Values
final_strain = strain_vector[-1]
with open(path + r'/report.txt', 'w') as report:
    report.writelines('Final CRSS: \n\n')

    report.writelines(f'CRSS-mean at 0% Strain: {float(np.mean(crss_start)) / 1e6:.4} MPa \n')
    report.writelines(f'CRSS-mean at {strain_96:.2%} Strain: {float(np.mean(crss_96)) / 1e6:.4} MPa \n')
    report.writelines(f'CRSS-mean at {strain_146:.2%} Strain: {float(np.mean(crss_146)) / 1e6:.4} MPa \n')
    report.writelines(f'CRSS-mean at {strain_216:.2%} Strain: {float(np.mean(crss_216)) / 1e6:.4} MPa \n')


# Plot Comparison
real_values = {0 : 226, 9.6*2.75/100 : 292, 14.6*2.75/100 : 323, 21.6*2.75/100 : 386}
min = real_values[0]
for key, value in real_values.items():
    real_values[key] = (real_values[key] - min) #/ (max - min)

sim_values = {0 : float(np.mean(crss_start)),
             strain_96 : float(np.mean(crss_96)),
             strain_146 : float(np.mean(crss_146)),
             strain_216 : float(np.mean(crss_216))}


min = sim_values[0]
for key, value in sim_values.items():
    sim_values[key] = (sim_values[key] - min) / 1e6

print(sim_values)

plt.figure(figsize=(8,6))
plt.plot(real_values.keys(), real_values.values(), label='Experimental')
plt.scatter(real_values.keys(), real_values.values(), marker='o', s=180)

plt.plot(sim_values.keys(), sim_values.values(), label='Simulative', linestyle='--')
plt.scatter(sim_values.keys(), sim_values.values(), marker='x', s=180)

plt.xlabel('Shear Strain', fontsize=16)
plt.ylabel('Scaled CRSS', fontsize=16)
plt.ylim([0, 300])
plt.legend(fontsize=18)

plt.savefig(path + r'/StrainHardening.png')"""

