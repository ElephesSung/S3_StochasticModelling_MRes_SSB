#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 14:22:37 2023

@author: Elephes
"""

from numpy.random import rand 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd  


def generate_transcription_rates(num_TXNrates, nominal_value, min_order_of_magnitude, max_order_of_magnitude):
    random_factors = np.random.uniform(low=min_order_of_magnitude, high=max_order_of_magnitude, size=num_TXNrates)
    transcription_rates = nominal_value * 10 ** random_factors
    return transcription_rates

num_TXNrates = 500
nominal_value = 1
min_order_of_magnitude = -3  # Corresponds to 10^(-3) = 0.001
max_order_of_magnitude = 0   # Corresponds to 10^0 = 1

transcription_rates = generate_transcription_rates(num_TXNrates, nominal_value, min_order_of_magnitude, max_order_of_magnitude)

mean_mrna = []
mean_p = []

var_mrna = []
var_p =[]

std_mrna = []
std_p = []

noise_mrna = []
noise_p = []

fano_mrna = []
fano_p =[]


for TXN in transcription_rates:
    k0 = TXN
    k1 = 0.01
    k2 = 5
    k3 = 1

    stoichiometry = [1, -1, 1, -1]
    
    mRNA_SS=[]
    Protein_SS=[]
    

    
    Cell_No_Each_TXN = range(500)
    
    
    for bibi in Cell_No_Each_TXN:
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
        dt = 10
        mRNA = 0
        P = 0
        t = 0.0
        end_time = 60
        num_process = int(end_time/dt)

        mRNA_level = []
        Protein_level = []
        time = [t]

        for i in range(num_process):
            mRNA, P = ssa_step(mRNA, P, t, t+dt)
            t += dt
            mRNA_level.append(mRNA)
            Protein_level.append(P)
            time.append(t)
            
        mRNA_SS.append(mRNA_level[-1])
        Protein_SS.append(Protein_level[-1])
    
    mean_mRNA_SS = np.mean(mRNA_SS)
    var_mRNA_SS = np.var(mRNA_SS)
    fano_mRNA_SS = var_mRNA_SS / mean_mRNA_SS
    
    mean_Protein_SS = np.mean(Protein_SS)
    var_Protein_SS = np.var(Protein_SS)
    fano_Protein_SS = var_Protein_SS / mean_Protein_SS

    fano_mrna.append(fano_mRNA_SS)
    fano_p.append(fano_Protein_SS)

    std_mRNA_SS = np.sqrt(var_mRNA_SS)
    std_Protein_SS = np.sqrt(var_Protein_SS)
    
    cv_mRNA_SS = (std_mRNA_SS / mean_mRNA_SS) * 100
    cv_Protein_SS = (std_Protein_SS / mean_Protein_SS) * 100
    
    noise_mRNA_SS = var_mRNA_SS / (mean_mRNA_SS ** 2)
    noise_Protein_SS = var_Protein_SS / (mean_Protein_SS ** 2)
    
    mean_mrna.append(mean_mRNA_SS)
    mean_p.append(mean_Protein_SS)
    
    noise_mrna.append(noise_mRNA_SS)
    noise_p.append(noise_Protein_SS)
    
    var_mrna.append(var_mRNA_SS)
    var_p.append(var_Protein_SS)
    
    std_mrna.append(std_mRNA_SS)
    std_p.append(std_Protein_SS)

data = {
    'Transcription Rate': transcription_rates,
    'Mean mRNA': mean_mrna,
    'Var mRNA': var_mrna,
    'Std mRNA': std_mrna,
    'Noise mRNA': noise_mrna,
    'Fano mRNA': fano_mrna,
    'Mean Protein': mean_p,
    'Var Protein': var_p,
    'Std Protein': std_p,
    'Noise Protein': noise_p,
    'Fano Protein': fano_p,
}

# Create a DataFrame from the dictionary
df = pd.DataFrame(data)

# Save the DataFrame to an Excel file
excel_filename = 'simulation_results.xlsx'
df.to_excel(excel_filename, index=False)

    


# Scatter plot of mean_mrna and noise_mrna with hollow, larger circles
plt.figure(figsize=(10, 8), dpi = 3000)
plt.scatter(mean_mrna, noise_mrna, label='mRNA', color='blue', edgecolors='C4', s=100, facecolors='none')
plt.xscale('log') 
plt.yscale('log') 
plt.xlabel('Mean mRNA Number', fontsize=26)
plt.ylabel('mRNA Noise', fontsize=26)
plt.tick_params(axis='both', labelsize=18)
plt.grid(True)
plt.show()


# Scatter plot of mean_p and noise_p with hollow, larger circles
plt.figure(figsize=(10, 8), dpi = 3000)
plt.scatter(mean_p, noise_p, label='Protein', color='orange', edgecolors='C9', s=100, facecolors='none')
plt.xscale('log')  # Set x-axis to logarithmic scale
plt.yscale('log') 
plt.xlabel('Mean Protein Number', fontsize=26)
plt.ylabel('Protein Noise', fontsize=26)
plt.tick_params(axis='both', labelsize=18)
plt.grid(True)
plt.show()


# Histogram of transcription_rates
plt.figure(figsize=(10, 8), dpi = 3000)
log_bins = np.logspace(np.log10(min(transcription_rates)), np.log10(max(transcription_rates)), 30)
plt.hist(transcription_rates, bins=log_bins, alpha=0.5, label='Transcription Rates', color='DarkRed', edgecolor='black')
plt.xscale('log') 
plt.xlabel('Transcription Rates', fontsize=30)
plt.ylabel('Probability Density (%)', fontsize=30)
plt.tick_params(axis='both', labelsize=23)
plt.grid(True)
plt.show()

# Histogram of mean_mrna
plt.figure(figsize=(10, 8), dpi = 3000)
log_bins = np.logspace(np.log10(min(mean_mrna)), np.log10(max(mean_mrna)), 30)
plt.hist(mean_mrna, bins=log_bins, alpha=0.5, color='Crimson', edgecolor='black')
plt.xscale('log') 
plt.xlim(10 ** -1.3, 10 ** 2.3)
plt.xlabel('Mean mRNA Number', fontsize=30)
plt.ylabel('Probability Density (%)', fontsize=30)
plt.tick_params(axis='both', labelsize=23)
plt.grid(True)
plt.show()


# Histogram of mean_p
plt.figure(figsize=(10, 8), dpi = 3000)
plt.hist(mean_p, bins=30, alpha=0.5, color='Salmon', edgecolor='black')
plt.xlabel('Mean Protein Number', fontsize=30)
plt.ylabel('Probability Density (%)', fontsize=30)
plt.tick_params(axis='both', labelsize=23)
plt.axvline(np.mean(mean_p), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(mean_p):.4f}')
plt.grid(True)
plt.legend(fontsize=18)
plt.show()


# Histogram of fano_mrna
plt.figure(figsize=(10, 8), dpi = 3000)
plt.hist(fano_mrna, bins=30, alpha=0.5, color='C4', edgecolor='black')
plt.xlabel('mRNA Fano Factor', fontsize=30)
plt.ylabel('Probability Density (%)', fontsize=30)
plt.tick_params(axis='both', labelsize=23)
plt.axvline(np.mean(mean_p), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(fano_mrna):.4f}')
plt.grid(True)
plt.legend(fontsize=18)
plt.show()



# Histogram of fano_p
plt.figure(figsize=(10, 8), dpi = 3000)
plt.hist(fano_p, bins=30, alpha=0.5, color='C4', edgecolor='black')
plt.xlabel('Protein Fano Factor', fontsize=30)
plt.ylabel('Probability Density (%)', fontsize=30)
plt.tick_params(axis='both', labelsize=23)
plt.axvline(np.mean(mean_p), color='red', linestyle='dashed', linewidth=2, label=f'Mean: {np.mean(fano_p):.4f}')
plt.grid(True)
plt.legend(fontsize=18)
plt.show()
    


    
    
    
    
    
    
    
    
    