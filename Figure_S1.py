#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 09:22:04 2023

@author: Elephes
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


k0 = 0.2
k1 = 0.01
k2 = 5
k3 = 1

init_mRNA = 0
init_Protein = 0

initial_conditions = [init_mRNA, init_Protein]


def model(population, t):
    mRNA, Protein = population
    dmRNA_dt = k0 - k1 * mRNA
    dProtein_dt = k2 * mRNA - k3 * Protein
    return [dmRNA_dt, dProtein_dt]

t_end = int(700)


t = np.linspace(0, t_end, 10000)  


population_over_time = odeint(model, initial_conditions, t)


mRNA_over_time, Protein_over_time = population_over_time.T



plt.figure(figsize=(10, 6), dpi=3000)
plt.plot(t, mRNA_over_time, label='mRNA', linewidth=3)
plt.plot(t, Protein_over_time, label='Protein', linewidth=3)

plt.xlabel('Time (t)', fontsize=21)
plt.xlim(0, t_end)
plt.ylabel('Number of Molecules', fontsize=21)
plt.ylim(0, )
plt.legend()
plt.legend(fontsize=18)
plt.tick_params(axis='both', labelsize=16)


t_label = 600
plt.axvline(x=t_label, color='gray', linestyle='--')

# Add labels at t = 700 on the y-axis
mRNA_at_t_label = np.interp(t_label, t, mRNA_over_time)
Protein_at_t_label = np.interp(t_label, t, Protein_over_time)
plt.text(t_label, mRNA_at_t_label, f'mRNA at t={t_label}: {mRNA_at_t_label:.2f}', fontsize=15, ha='right')
plt.text(t_label, Protein_at_t_label, f'Protein at t={t_label}: {Protein_at_t_label:.2f}', fontsize=15, ha='right')

# Add horizontal dashed line at mRNA and Protein values
plt.axhline(y=mRNA_at_t_label, color='gray', linestyle='--')
plt.axhline(y=Protein_at_t_label, color='gray', linestyle='--')

plt.grid(True)


plt.show()



