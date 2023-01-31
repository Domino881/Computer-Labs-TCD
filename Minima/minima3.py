import matplotlib.pyplot as plt
import math
import numpy as np

c = 1.44
A = 1090.0
p = 0.033

x = np.arange(0.01,1.0,0.005)
    
def V(x):
    return A*np.exp(-x/p)-c/x

def dV(x):
    return -A*(1/p)*np.exp(-x/p)+c/(x**2)

def ddV(x):
    return A*(1/p**2)*np.exp(-x/p)-2*c/(x**3)

f, (plt1,plt2) = plt.subplots(2,1, sharex=True, figsize=(8,4))
plt2.plot(x,-dV(x), label="$F(x)=\\frac{-dV(x)}{dx}$")
plt2.plot(x,np.zeros(len(x)))
plt2.legend()
plt2.grid()
plt2.set(ylim=(-30,50))


def run(plotting, x0, tol):
    x1 = x0
    nsteps=0
    while abs(dV(x1)) > tol:
        if(plotting):
            plt.plot(x1, -dV(x1), color=(abs(math.atan(nsteps))/1.6,0,0),
                 marker='s')
            #plt.plot(x1, -dV(x1), 'rs')
            
        x1 = x1 - dV(x1)/ddV(x1)
        nsteps += 1

    if(plotting):
        print(x1,-dV(x1),V(x1))
    
    return x1

minim = run(True,0.2,0.0001)

plt1.plot(x,V(x), label="$V(x)$")
plt1.plot(minim, V(minim), '')
plt1.plot(x,np.zeros(len(x)),c='C1')

plt1.set(xlim=(x[0],x[-1]), ylim=(-6,4), xticks=np.arange(x[0]-0.01,x[-1],0.1))
plt1.grid()
plt1.legend(loc='upper right')

plt1.vlines(x=0.236053, ymin=-6, ymax=4, clip_on=False, linewidth=1,
            linestyles='dashed', color='grey')
plt2.vlines(x=0.236053, ymin=-30, ymax=65, clip_on=False, linewidth=1,
            linestyles='dashed', color='grey')

plt1.plot(0.236053, V(0.236053), 'o', c='C0', label="minimum")
plt2.plot(0.236053, 0, 'o', c='C0')

plt.show()