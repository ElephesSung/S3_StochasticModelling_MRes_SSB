#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 10:03:34 2023

@author: SUNG, Elephes
Research Postgraduate, DoLS, FoNS, Imperial College
"""

from numpy.random import rand 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson


k0 = 0.2
k1 = 0.01
k2 = 5
k3 = 1

stoichiometry = [1, -1, 1, -1]

def propensities(mRNA, P):
    return [k0, k1*mRNA, k2*mRNA, k3*P]

def reaction_times(mRNA, P):
    a = propensities(mRNA, P)
    aInv = [1/s if s> 0 else np.inf for s in a]
    return -np.log(rand(4))*aInv

def ssa_step(mRNA, P, tIn, tOut):
    t=tIn
    while t<tOut:
        rt = reaction_times(mRNA, P)
        min_idx = np.argmin(rt)
        tau = rt[min_idx]
        if min_idx == 0 or min_idx == 1:
            mRNA += stoichiometry[min_idx]
        elif min_idx == 2 or min_idx == 3:
            P += stoichiometry[min_idx]
        t += tau 
    return mRNA, P

dt = 0.05
mRNA = 0
P = 0
t = 0.0
end_time = 5000
num_process = int(end_time/dt)

mRNA_level = [mRNA]
Protein_level = [P]
time = [t]

for i in range(num_process):
    mRNA, P = ssa_step(mRNA, P, t, t+dt)
    t += dt
    mRNA_level.append(mRNA)
    Protein_level.append(P)
    time.append(t)

mean_mRNA_expression = np.mean(mRNA_level)
mean_Protein_expression = np.mean(Protein_level)

V_mRNA_expression = np.var(mRNA_level)
V_Protein_expression = np.var(Protein_level)


fig, ax1 = plt.subplots(figsize=(17, 8), dpi=300)

ax1.set_xlabel('Time (t)', fontsize=30)
ax1.set_ylabel('mRNA Expression Level', color='b', fontsize=27)
ax1.set_xlim(0, end_time)
ax1.set_ylim(0, 75)
ax1.plot(time, mRNA_level, color='b', label='mRNA Level')
ax1.tick_params(axis='x', labelsize=20) 
ax1.tick_params(axis='y', labelcolor='b', labelsize=20)

ax2 = ax1.twinx()
ax2.set_ylabel('Protein Expression Level', color='C1', fontsize=27)
ax2.plot(time, Protein_level, color='C1', label='Protein Level')
ax2.tick_params(axis='y', labelcolor='C1', labelsize=20)
ax2.set_ylim(0, 200)

ax1.text(0.05, 0.95, f'Mean mRNA: {mean_mRNA_expression:.2f}\nVariance mRNA: {V_mRNA_expression:.2f}', 
         transform=ax1.transAxes, fontsize=15, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

ax2.text(0.05, 0.85, f'Mean Protein: {mean_Protein_expression:.2f}\nVariance Protein: {V_Protein_expression:.2f}', 
         transform=ax2.transAxes, fontsize=15, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))


ax1.grid(axis='x', linestyle='-', alpha=0.7)
plt.grid(True)
plt.xticks(np.arange(0, 5500, 500))
plt.show()



indices_to_delete = np.where(np.array(time) <= 50)
time_SS = np.delete(time, indices_to_delete)
mRNA_level_SS = np.delete(mRNA_level, indices_to_delete)

fano = np.var(mRNA_level_SS)/np.mean(mRNA_level_SS)

plt.figure(figsize=(10, 8))
plt.hist(mRNA_level_SS, bins=range(0, max(mRNA_level_SS)), density=True, alpha=0.7, color='C0', edgecolor='black')
plt.xlim(min(mRNA_level_SS), max(mRNA_level_SS))
plt.xlabel("Number of mRNA", fontsize = 20)
plt.ylabel("Probability Density (%)", fontsize = 20)
plt.tick_params(axis='both', labelsize=16)
#plt.axvline(np.mean(mRNA_level), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(fano_list):.4f}')
#mean_label = f'Mean: {np.mean(fano_list):.2f}'
#variance_label = f'Variance: {np.var(fano_list):.2f}'
#mean_label = f'Mean: {np.mean(fano_list):.2f}'
#variance_label = f'Variance: {np.var(fano_list):.2f}'
#middle_top_x = (max(fano_list) + min(fano_list)) / 2
#middle_top_y = 300

#plt.annotate(mean_label, xytext=None, xy=(np.mean(fano_list), 30000), fontsize=12, color='red', ha='center')
#plt.annotate(variance_label, xytext=None, xy=(np.mean(fano_list),3000), fontsize=12, color='red', ha='center')


mu = k0/k1
ar = np.arange(0, max(mRNA_level_SS))
plt.plot(ar + 0.5, poisson.pmf(ar, mu), color='C1', linewidth=3, label='Poisson Distribution')

fano_label = f'Fano Factor: {fano:.2f}'
plt.annotate(fano_label, xy=(0.85, 0.6), xycoords='axes fraction', fontsize=18, color='red', ha='center')

plt.legend(fontsize = 18)
plt.grid(True)
dpi_value = 3000
plt.show()


