#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:46:54 2023

@author: SUNG, Fan-Hsin Elephes 
Research Postgraduate. MRes in Systems Biology
Department of Life Sciences, Imperial College of Science, Technology and Medicine
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
        
        if min_idx == 0 or min_idx == 1:
            mRNA += stoichiometry[min_idx]
            tau = rt[min_idx]
        elif min_idx == 2 or min_idx == 3:
            t_new = tIn
            rt_new = rt[0:2]
            min_idx_new = np.argmin(rt_new)
            tau = rt_new[min_idx_new]
            tOut_new = t + tau
            P += stoichiometry[min_idx]
            t_new += rt[min_idx]
            while t_new < tOut_new:
                rt_next = reaction_times(mRNA, P)
                rt_next_new = rt_next[2:4]
                midx = np.argmin(rt_next_new)
                t_new += rt_next_new[midx]
                P += stoichiometry[midx+2]
            mRNA += stoichiometry[min_idx_new]  
        t += tau 
    return mRNA, P

cycle = int(input("How many sets of figures do you want? "))
for i in range(cycle):
    dt = 0.05
    mRNA = 0
    P = 0
    t = 0.0
    end_time = 1000
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


    
    indices_to_delete = np.where(np.array(time) <= 50)
    time_SS = np.delete(time, indices_to_delete)
    mRNA_level_SS = np.delete(mRNA_level, indices_to_delete)
    Protein_level_SS = np.delete(Protein_level, indices_to_delete)
    
    mean_mRNA_expression_SS = np.mean(mRNA_level_SS)
    mean_Protein_expression_SS = np.mean(Protein_level_SS)

    V_mRNA_expression_SS = np.var(mRNA_level_SS)
    V_Protein_expression_SS = np.var(Protein_level_SS)


    fig, ax1 = plt.subplots(figsize=(17, 8))


    ax1.set_xlabel('Time (t)', fontsize=30)
    ax1.set_ylabel('mRNA Expression Level', color='b', fontsize=27)
    ax1.set_xlim(0, end_time)
    ax1.set_ylim(0, 75)
    ax1.plot(time, mRNA_level, color='b')
    ax1.tick_params(axis='x', labelsize=20) 
    ax1.tick_params(axis='y', labelcolor='b', labelsize=20)


    ax2 = ax1.twinx()
    ax2.set_ylabel('Protein Expression Level', color='C1', fontsize=27)
    ax2.plot(time, Protein_level, color='C1')
    ax2.tick_params(axis='y', labelcolor='C1', labelsize=20)
    ax2.set_ylim(0, 220)

    ax1.text(0.05, 0.95, f'Mean mRNA_SS: {mean_mRNA_expression_SS:.2f}\nVariance mRNA_SS: {V_mRNA_expression_SS:.2f}', 
             transform=ax1.transAxes, fontsize=15, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

    ax2.text(0.05, 0.85, f'Mean Protein_SS: {mean_Protein_expression_SS:.2f}\nVariance Protein_SS: {V_Protein_expression_SS:.2f}', 
             transform=ax2.transAxes, fontsize=15, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))


    ax1.grid(axis='x', linestyle='-', alpha=0.7)
    plt.grid(True)
    plt.legend()
    plt.xticks(np.arange(0, 1100, 100))
    dpi_value = 3000
    plt.show()

  

    fano = np.var(mRNA_level_SS)/np.mean(mRNA_level_SS)
    
    plt.figure(figsize=(10, 8))
    plt.hist(mRNA_level_SS, bins=range(0, max(mRNA_level_SS)), density=True, alpha=0.7, color='C0', edgecolor='black')
    plt.xlim(min(mRNA_level_SS), max(mRNA_level_SS))
    plt.xlabel("Number of mRNA", fontsize = 20)
    plt.ylabel("Probability Density (%)", fontsize = 20)
    plt.tick_params(axis='both', labelsize=16)

    mu = k0/k1
    ar = np.arange(0, max(mRNA_level_SS))
    plt.plot(ar + 0.5, poisson.pmf(ar, mu), color='C1', linewidth=3, label='Poisson Distribution')

    fano_label = f'Fano Factor: {fano:.2f}'
    plt.annotate(fano_label, xy=(0.85, 0.6), xycoords='axes fraction', fontsize=18, color='red', ha='center')

    plt.legend(fontsize = 18)
    plt.grid(True)
    dpi_value = 3000
    plt.show()


    indices_to_delete = np.where(np.array(time) <= 50)
    time_SS = np.delete(time, indices_to_delete)
    Protein_level_SS = np.delete(Protein_level, indices_to_delete)

    fano_P = np.var(Protein_level_SS)/np.mean(Protein_level_SS)


    plt.figure(figsize=(10, 8))
    plt.hist(Protein_level_SS, bins=range(0, max(Protein_level_SS)), density=True, alpha=0.7, color='C4', edgecolor='black')
    plt.xlim(min(Protein_level_SS), max(Protein_level_SS))
    plt.xlabel("Number of Protein", fontsize = 20)
    plt.ylabel("Probability Density (%)", fontsize = 20)
    plt.tick_params(axis='both', labelsize=16)

    fano_label = f'Fano Factor (Protein): {fano_P:.2f}'
    plt.annotate(fano_label, xy=(0.75, 0.6), xycoords='axes fraction', fontsize=18, color='red', ha='center')

    plt.legend(fontsize = 18)
    plt.grid(True)
    dpi_value = 3000
    plt.show()

