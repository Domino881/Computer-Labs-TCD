# This code was written during the 2023 Theoretical Physics
# Student Association summer workshop.
#
# Its purpose is to simulate an Ising model using the 
# Markov Chain Monte Carlo method and estimate the energy
# and magnetisation of the system.

import numpy as np
from termcolor import colored
import matplotlib.pyplot as plt
import time

xdim, ydim = 10,10

class Lattice():
    def __init__(self, matrix):
        self.matrix = matrix
        self.expanded = self.matrix
        self.expanded = np.append(self.expanded, [self.matrix[0]], 0)
        self.expanded = np.append([self.matrix[-1]], self.expanded, 0)
        self.expanded = np.insert(self.expanded, [0], [[0]]+[[bb[-1]] for bb in self.matrix]+[[0]], 1)
        self.expanded = np.insert(self.expanded, [xdim+1], [[0]]+[[bb[0]] for bb in self.matrix]+[[0]], 1)
        
        
    def __getitem__(self, tup):
        ix,iy = tup
        ix = ix if ix<len(self.matrix) else ix-len(self.matrix)
        iy = iy if iy<len(self.matrix[0]) else iy-len(self.matrix[0])
        
        return self.matrix[ix][iy]
    
    def __setitem__(self, tup, newval):
        ix, iy = tup
        ix = ix if ix<len(self.matrix) else ix-len(self.matrix)
        iy = iy if iy<len(self.matrix[0]) else iy-len(self.matrix[0])
        
        self.matrix[ix][iy] = newval
        self.expanded[ix+1][iy+1] = newval
        
        
    def nn(self, tup):
        ix, iy = tup
        return [self[ix + dx, iy + dy] for dx,dy in [[-1,0],[1,0],[0,-1],[0,1]]]
    
    def __repr__(self):
        for row in self.matrix:
            print(''.join(f'{colored(f"{x:3}","magenta" if x>0 else "cyan")}' for x in row))
        return "E = " + str(self.H()) + ",   M = " + str(self.M())
    
    def H(self):
        E = 0 
        
 
        E += np.sum(self.expanded[0:-2, 1:-1] * self.matrix)
        E += np.sum(self.expanded[2:, 1:-1] * self.matrix)
        E += np.sum(self.expanded[1:-1, 0:-2] * self.matrix)
        E += np.sum(self.expanded[1:-1, 2:] * self.matrix)
        # for ix in range(xdim):
        #     for iy in range(ydim):
        #         E += 0.5 * np.sum(self[(ix,iy)] * self.nn((ix,iy)))

        return -E/2
    
    def M(self):
        return np.sum(self.matrix)

burnin = 100
N = 500


E_avg_s=[]
M_avg_s=[]
C_v_s=[]

T_s = np.linspace(0.8, 7, 30)

for T in T_s:
    
    beta = 1 / T
    l = Lattice(np.random.choice([-1,1], (xdim,ydim)))

    
    energies=[]
    magnetisations=[]
    
    start = time.time()

    for n in range(N):  
        for ix in range(xdim):
            for iy in range(ydim):
                idx = (ix,iy)
                deltaH = 2*l[idx] * np.sum(l.nn(idx))
                
                if deltaH <= 0:
                    l[idx] = -l[idx]
                else:
                    r = np.random.rand()
                    if r <= np.exp(-beta * deltaH):
                        l[idx] = -l[idx]
            
        En = l.H()
        Mn = l.M()
        energies.append(En)
        magnetisations.append(Mn)
        
        # print(l)
    print(f"Time per sweep: {(time.time()-start)/N*1e6:.2f}us")
    
    
    E_avg = np.mean(energies[burnin:])
    M_avg = np.mean(magnetisations[burnin:])
    
    C_v = beta**2 * np.sum((energies - E_avg)**2) / (N-burnin-1)
    
    E_avg_s.append(E_avg / (xdim*ydim))
    M_avg_s.append(abs(M_avg / (xdim*ydim)))
    C_v_s.append(C_v / (xdim*ydim))


F, axes = plt.subplots(1,2, figsize=(15,4),sharex=True)

axes[0].plot(energies, label="E")
axes[0].legend()
axes[1].plot(magnetisations, label="M")
axes[1].legend()

tttmp = plt.figure(1)
F, axes = plt.subplots(3,1, figsize=(8,8),sharex=True)

axes[0].plot(T_s, E_avg_s, 'o', label="normalized E_avg")
axes[0].legend()
axes[1].plot(T_s, M_avg_s, 'o', label="normalized M_avg")
axes[1].legend()
axes[2].plot(T_s, C_v_s, 'o', label="normalized C_v")
axes[2].legend()

list_loc=[0.0,2.27,5.0]
axes[0].xaxis.set_ticks(list_loc)

print(l)


print(f"E_avg = {E_avg:.2f},  M_avg = {M_avg:.2f},  C_v = {C_v:.2f}")
