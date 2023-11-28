#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:55:16 2023

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as plt

# Given parameters
k0 = 0.2
k1 = 0.01
initial_mRNA = 0
simulation_time = 800
time_step = 1  # Time step for Euler's method

# Function representing the differential equation
def mRNA_dynamics(mRNA, t):
    return k0 - k1 * mRNA

# Euler's method for numerical integration
def euler_method(func, initial_value, time_range, step_size):
    values = [initial_value]
    times = np.arange(0, time_range, step_size)
    for t in times[:-1]:
        next_value = values[-1] + step_size * func(values[-1], t)
        values.append(next_value)
    return times, values

# Perform simulation
time_points, mRNA_values = euler_method(mRNA_dynamics, initial_mRNA, simulation_time, time_step)



# Plot the results
plt.figure(figsize=(8, 6), dpi = 3000)
plt.plot(time_points, mRNA_values)
plt.xlabel('Time (s)', fontsize = 20)
plt.xlim(0, 800)
plt.ylabel('Number of mRNA', fontsize = 20)
plt.ylim(0, 25)
plt.tick_params(axis='both', labelsize=16)

t_label = 700
mRNA_at_t_label = np.interp(t_label, time_points, mRNA_values)

plt.axvline(x=t_label, color='grey', linestyle='--', linewidth=3, label='t=700')
plt.axhline(y=mRNA_at_t_label, color='grey', linestyle='--', linewidth=3, label=f'mRNA at t={t_label}')
plt.text(t_label-0.5, mRNA_at_t_label+0.5, f'Number of mRNA at t={t_label}: {mRNA_at_t_label:.2f}', fontsize=18, ha='right')

plt.grid(True)
plt.show()
