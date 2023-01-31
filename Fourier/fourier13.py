import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import square
plt.rcParams['figure.dpi']=200


# maxk=25

def Simpson(f, a, b, n=500):
    h = float(b-a) / n
    hs = np.array([2*j*h for j in range(1,n//2+1)])
    # print(a+hs[:-1], a+hs-h)
    return h/3 * (f(a) + 2*sum(f(a+hs[:-1])) + 4*sum(f(a+hs-h)) + f(b))


def FourierSeries(f, T):
    omega = 2*np.pi/T
    
    a = [1/T * Simpson(f,0,T)]
    for k in range(1,maxk+1):
        g = lambda x: f(x)*np.cos(x*k*omega)
        a.append(2/T*Simpson(g, 0, T))
        
    b = [0]
    for k in range(1,maxk+1):
        h = lambda x: f(x)*np.sin(x*k*omega)
        b.append(2/T*Simpson(h, 0, T))
        
    return a, b

def Fourier(f, T):
    a,b = FourierSeries(f, T)
    omega = 2*np.pi/T
    g = lambda x: np.sum([a[k]*np.cos(k*x*omega)+b[k]*np.sin(k*x*omega)\
                          for k in range(0, maxk)], axis=0)
    return g

x = np.arange(0,100,0.01)
T = np.pi*2

F, plots = plt.subplots(3,3, figsize=(25,15))
i=0
F.suptitle("Reconstructed square wave for different numbers of Fourier terms",
           fontsize=37, fontweight=200)

def func(x):
    return square(x)

for maxk in [1,2,3,5,10,20,30,50,200]:
    t = np.arange(0,2*T,0.01)
    plots[i//3][i%3].set(title=r"$k_{max}=$"+str(maxk))
    plots[i//3][i%3].plot(t, func(t), "gray")
    plots[i//3][i%3].plot(t, Fourier(func, T)(t), "C1", linewidth=3)

    plots[i//3][i%3].set(xlim=[0,4*np.pi],ylim=[-1.5,1.5])
    i += 1
    
    
F, plots = plt.subplots(3,3, figsize=(25,15))
i=0
F.suptitle("Reconstructed rectangle wave for different numbers of Fourier terms",
           fontsize=37, fontweight=200)
    
def func(x):
    x = x-(x//(2*np.pi)*2*np.pi)
    return 1 - 2*(x < 2/5*np.pi)


for maxk in [1,2,3,5,10,20,30,50,200]:
    t = np.arange(0,2*T,0.01)
    plots[i//3][i%3].set(title=r"$k_{max}=$"+str(maxk))
    plots[i//3][i%3].plot(t, func(t), "gray")
    plots[i//3][i%3].plot(t, Fourier(func, T)(t), "C1", linewidth=3)

    plots[i//3][i%3].set(xlim=[0,4*np.pi],ylim=[-1.5,1.5])
    i += 1
