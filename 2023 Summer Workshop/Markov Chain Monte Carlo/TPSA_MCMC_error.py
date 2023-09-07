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


N = 2**14
burnin = 128

beta = 1 / 5.0
l = Lattice(np.random.choice([-1,1], (xdim,ydim)))
e = l.H()
m = l.M()

energies=[]
magnetisations=[]

start = time.time()

for n in range(N+burnin):  
    for ix in range(xdim):
        for iy in range(ydim):
            idx = (ix,iy)
            deltaH = 2*l[idx] * np.sum(l.nn(idx))
            
            if deltaH <= 0:
                l[idx] = -l[idx]
                e += deltaH
                m += 2 * l[idx]
            else:
                r = np.random.rand()
                if r <= np.exp(-beta * deltaH):
                    l[idx] = -l[idx]
                    e += deltaH
                    m += 2 * l[idx]
        
    # En = l.H()
    # Mn = l.M()
    energies.append(e)
    magnetisations.append(m)
    print(f"{n}/{N+burnin}")
    
print(f"Time per sweep: {(time.time()-start)/N*1e6:.2f}us")


E_avg = np.mean(energies[burnin:])
M_avg = np.mean(magnetisations[burnin:])

tmptmp=plt.figure(1)

binsizes = np.logspace(0,13,base=2.0, num=25, dtype=int)
variances=[]

for binsize in binsizes:
    binned = [np.mean(energies[burnin+i:burnin+i+binsize]) for i in range(0,N-binsize+1,binsize)]
    variances.append(np.var(binned)/(N/binsize))
    
plt.plot(binsizes, variances)
plt.xscale('log', base=2.0)

C_v = beta**2 * np.sum((energies - E_avg)**2) / (N-burnin-1)