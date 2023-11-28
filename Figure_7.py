#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 09:01:18 2023

@author: Elephes
"""

import numpy as np
from numpy.random import rand 
import matplotlib.pyplot as plt
from scipy.stats import poisson

k0 = 0.2
k1 = 0.01

end_time = int(1000)
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


dt = 0.05
x = 0
t = 0.0

mrna = [x]
time = [t]
mRNA=[]

num_process = int(end_time/dt)

for bibi in range(19000):
    for i in range(num_process):
        x = ssa_step(x, t, t+dt)
        t+=dt
        mrna.append(x)
    mRNA.append(mrna[-1])

m = np.mean(mRNA)
v = np.var(mRNA)

plt.figure(figsize=(8, 6))
plt.hist(mRNA, bins=range(0, max(mRNA)), density=True, alpha=0.7, color='C0', edgecolor='black')
plt.xlim(min(mRNA), max(mRNA))
plt.xlabel("Number of mRNA at 1000s", fontsize = 20)
plt.ylabel("Probability Density (%)", fontsize = 20)
plt.tick_params(axis='both', labelsize=16)
plt.axvline(np.mean(mRNA), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(mRNA):.4f}')


mu = k0/k1
ar = np.arange(min(mRNA), max(mRNA))
plt.plot(ar+0.5, poisson.pmf(ar, mu), color='C1', linewidth=3, label = 'Poisson Distribution')
plt.legend(fontsize=15)


plt.grid(True)
dpi_value = 4000
plt.show()



