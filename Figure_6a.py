#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 15:52:08 2023

@author: Elephes
"""
from numpy.random import rand 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon
from scipy.stats import poisson
from scipy.stats import norm
from scipy.optimize import curve_fit

k0 = 0.2
k1 = 0.01

stoichiometry = [1, -1]


def propensities(x):
    return [k0, k1*x]

def reaction_times(x):
    a = propensities(x)
    aInv = [1/s if s> 0 else np.inf for s in a]
    return -np.log(rand(2))*aInv

def ssa_step(x, tIn, tOut):
    t=tIn
    while t<tOut:
        rt = reaction_times(x)
        minidx = np.argmin(rt)
        tau = rt[minidx]
        x += stoichiometry[minidx]
        t += tau 
    return x 



fano_list = []

for b in range(5000):
    dt = 0.05
    x = 0
    t = 0.0

    mrna = [x]
    time = [t]
    end_time = int(10000)

    num_process = int(end_time/dt)
    
    for i in range(num_process):
        x = ssa_step(x, t, t+dt)
        t+=dt
        mrna.append(x)
        time.append(t)
    
    indices_to_delete = np.where(np.array(time) <= 50)
    time = np.delete(time, indices_to_delete)
    mrna = np.delete(mrna, indices_to_delete)

    m = np.mean(mrna)
    v = np.var(mrna)
    fano = v/m
    fano_round = round(fano, 5)
    fano_list.append(fano_round)

print(fano_list)


plt.figure(figsize=(8, 6))
plt.hist(fano_list, bins=20, density=True, alpha=0.7, color='C0', edgecolor='black')
plt.xlim(min(fano_list), max(fano_list))
plt.xlabel("Fano Factor", fontsize = 20)
plt.ylabel("Probability Density (%)", fontsize = 20)
plt.tick_params(axis='both', labelsize=15)
plt.axvline(np.mean(fano_list), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(fano_list):.4f}')
mean_label = f'Mean: {np.mean(fano_list):.2f}'
variance_label = f'Variance: {np.var(fano_list):.2f}'
mean_label = f'Mean: {np.mean(fano_list):.2f}'
variance_label = f'Variance: {np.var(fano_list):.2f}'
middle_top_x = (max(fano_list) + min(fano_list)) / 2
middle_top_y = 300

plt.annotate(mean_label, xytext=None, xy=(np.mean(fano_list), 30000), fontsize=12, color='red', ha='center')
plt.annotate(variance_label, xytext=None, xy=(np.mean(fano_list),3000), fontsize=12, color='red', ha='center')

plt.grid(True)
dpi_value = 3000
plt.show()


plt.annotate(variance_label, xytext=None, xy=(1.1, 10), fontsize=18, color='red', ha='center')
plt.grid(True)
dpi_value = 3000
plt.show()


