#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 07:34:29 2023

@author: Elephes
"""

import numpy as np
from numpy.random import rand 
import matplotlib.pyplot as plt

k0 = 0.2
k1 = 0.01

end_time = int(200)
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

# Calculate mean in intervals of 10 seconds
interval_means = []
interval_times = range(0, end_time + 1, 10)
for start_time in interval_times:
    end_time_interval = start_time + 10
    indices = [i for i, val in enumerate(time) if start_time <= val < end_time_interval]
    interval_mean = np.mean([mrna[i] for i in indices])
    interval_means.append(interval_mean)


#plot the trajectory
plt.figure(figsize=(10, 8))
plt.plot(time, mrna)
plt.xlabel('Time (s)', fontsize=25)
plt.xlim(0, 100)
plt.xticks(np.arange(0, 110, 10))
plt.ylabel('Number of mRNA molecules', fontsize=25)
plt.ylim(0, 45)
plt.tick_params(axis='both', labelsize=20)
plt.grid(True)
plt.show()


# Plot the means
plt.figure(figsize=(10, 8))
plt.bar(interval_times, interval_means, color='Salmon', width=10, align='edge', edgecolor='black')  # Adjust the color and width as needed
plt.xlabel('Time (s)', fontsize=25)
plt.ylabel('Mean mRNA molecules', fontsize=25)
plt.xlim(0, 200)
plt.xticks(np.arange(0, end_time + 20, 20))
plt.ylim(10, 30)
plt.tick_params(axis='both', labelsize=20)
plt.grid(axis='y')
plt.show()

