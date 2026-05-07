# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:08:16 2026

@author: laure
"""
import matplotlib.pyplot as plt
import pyspike as spk
from numpy import *

tmin, tmax = 0, 10
num_trains = 3
spikes = [
    [0, 1, 2, 4, 7, 10],
    [0, 3, 4, 6, 10],
    [0, 2, 5, 10]
]

# dataset = 8

# if dataset == 8:
#    tmax=10;
#    tmin=0
#    spike_trains = []
#    spike_trains.append(spk.SpikeTrain([0,1.9, 3.9, 6,10], [tmin, tmax]))
#    spike_trains.append(spk.SpikeTrain([0,2, 6.1,10], [tmin, tmax]))
#    spike_trains.append(spk.SpikeTrain([0,2.1, 4.1, 6.2,10], [tmin, tmax]))
#    spike_trains.append(spk.SpikeTrain([0,2.2, 4.2, 5.9,10], [tmin, tmax]))
#    num_trains = len(spike_trains)

def f_spike_conform (spikes, tmin, tmax):
    n= len(spikes)
    for i in range (0,n):
        if spikes[i][0]!=tmin:
            spikes[i].insert(0,tmin)
        if spikes[i][-1]!=tmax:
            spikes[i].append(tmax)
    return spikes

def f_ISI_calculate (spikes):
    ISI_distance = []
    for i in range (len(spikes)):
        ISI_distance_i = []
        for j in range (len(spikes[i])-1):
            ISI_distance_i.append(spikes[i][j+1]-spikes[i][j])
        ISI_distance.append(ISI_distance_i)
    return ISI_distance

def I_t (value1, value2):
    return abs (value1 - value2) / max(value1, value2)

def f_ISI_distance (spikes, tmin, tmax):
     
     n = len (spikes)
     isi_values = f_ISI_calculate(spikes)
     
     if n == 1:
        plt.figure()
        # On utilise les temps de spikes du premier train pour l'axe X
        # On ajoute la dernière valeur de ISI pour que la taille corresponde à spikes[0]
        y_plot = isi_values[0] + [isi_values[0][-1]]
        
        plt.step(spikes[0], y_plot, where='post')
        
        plt.title(f"Single Spike Train ISI Evolution")
        # dist_val n'existe pas ici, on calcule la moyenne des ISI directement
        avg_isi = sum(isi_values[0]) / (tmax - tmin)
        plt.suptitle(f"Average ISI: {avg_isi:.4f}")
        
        plt.ylim(0, max(y_plot) * 1.1) # On adapte l'échelle Y
        plt.show()
        return avg_isi
    
     ISI_matrix = zeros((n, n))
     pair_data = [] # Pour stocker les évolutions temporelles
     all_t_events = [tmin, tmax]
    
     # Parcours des paires
     for i in range(n):
         for j in range(i + 1, n):
            # Union des temps de spikes (équivalent de unique([spikes{i}, spikes{j}]))
             t_all = sorted(list(set(spikes[i]) | set(spikes[j])))
             all_t_events.extend(t_all)
            
             i_ij_acc = 0 # Intégrale
             it_list = [] # Évolution temporelle
        
             for k in range(len(t_all) - 1):
                t_mid = (t_all[k] + t_all[k+1]) / 2.0
                
                idx_x = len([s for s in spikes[i][:-1] if s <= t_mid]) - 1
                idx_y = len([s for s in spikes[j][:-1] if s <= t_mid]) - 1
                
                val_x = isi_values[i][idx_x]
                val_y = isi_values[j][idx_y]
                
                current_i_t = I_t(val_x, val_y)
            
                i_ij_acc += current_i_t * (t_all[k+1] - t_all[k])
                it_list.append(current_i_t)
        
             dist_val = i_ij_acc / (tmax - tmin)
             ISI_matrix[i, j] = dist_val
             ISI_matrix[j, i] = dist_val
        
             pair_data.append({'t': t_all, 'It': it_list})
            
             plt.figure()
             plt.step(t_all, it_list + [it_list[-1]], where='post')
             plt.title(f"Pair {i+1} vs {j+1}")
             plt.suptitle(f"Pairwise ISI-distance: {dist_val:.4f}")            
             plt.ylim(0, 1)
             plt.margins(x=0, y=0)
             plt.show()
             
     t_global = unique(all_t_events)
     I_matrix = zeros((len(pair_data), len(t_global) - 1))
   
     for p in range(len(pair_data)):
         tp = pair_data[p]['t']
         itp = pair_data[p]['It']
         for k in range(len(t_global) - 1):
             t_mid = (t_global[k] + t_global[k+1]) / 2.0
             idx = where(tp[:-1] <= t_mid)[0][-1]
             I_matrix[p, k] = itp[idx]
           
     I_pop_mean = mean(I_matrix, axis=0)
     dt = diff(t_global) 
     global_isi_val = sum(I_pop_mean * dt) / (tmax - tmin)


     # Plot global
     plt.figure()
     plt.step(t_global, append(I_pop_mean, I_pop_mean[-1]), where='post', color='red')
     plt.title("Evolution of Population Average ISI distance")
     plt.suptitle(f"Global ISI-distance: {global_isi_val:.4f}")   
     plt.xlim(0, 10)
     plt.ylim(0, 1)
     plt.show()
     
     plt.figure()
     plt.imshow(ISI_matrix, cmap='viridis', extent=[0.5, n+0.5, n+0.5, 0.5])
     plt.xticks(range(1, n + 1)); 
     plt.yticks(range(1, n + 1));
     plt.xlabel('Spike_trains');
     plt.ylabel('Spike_trains');
     plt.colorbar()
     plt.title("Matrix of the ISI-distance")
     plt.show()
    
     return ISI_matrix, pair_data, sorted(list(set(all_t_events)))
  
f_ISI_distance(spike_trains,0,10)   