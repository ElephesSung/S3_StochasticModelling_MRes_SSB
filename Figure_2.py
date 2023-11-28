#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 08:50:48 2023

@author: SUNG, Elephes
Research Postgraduate, DoLS, FoNS, Imperial College
"""
from numpy.random import rand 
#Create an array of the given shape and populate it with random samples from a uniform distribution over [0, 1).
import numpy as np
import matplotlib.pyplot as plt


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

num_process = int(end_time/dt)

for i in range(num_process):
    x = ssa_step(x, t, t+dt)
    t+=dt
    mrna.append(x)
    time.append(t)


plt.figure(figsize=(20, 8))
plt.plot(time, mrna)
plt.xlabel('Time (s)', fontsize=25)
plt.xlim(0, end_time)
plt.xticks(np.arange(0, 1050, 50))
plt.ylabel('Number of mRNA molecules', fontsize=25)
plt.ylim(0, 45)
plt.tick_params(axis='both', labelsize=20)
plt.grid(True)
plt.show()



