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

dt = 0.05
x = 0
t = 0.0

mrna = [x]
time = [t]
end_time = int(5000)

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


plt.figure(figsize=(8, 6))
plt.hist(mrna, bins=range(0, max(mrna)), density=True, alpha=0.7, color='C0', edgecolor='black')
plt.xlim(min(mrna), max(mrna))
plt.xlabel("Number of mRNA", fontsize = 20)
plt.ylabel("Probability Density (%)", fontsize = 20)
plt.tick_params(axis='both', labelsize=16)
plt.axvline(np.mean(mrna), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(mrna):.4f}')

mu = k0/k1
ar = np.arange(min(mrna), max(mrna))
plt.plot(ar+0.5, poisson.pmf(ar, mu), color='C1', linewidth=3, label = 'Poisson Distribution')
plt.legend(fontsize=15)


plt.grid(True)
dpi_value = 4000
plt.show()



