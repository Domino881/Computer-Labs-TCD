import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import square
plt.rcParams['figure.dpi']=200


maxk=25

def Simpson(f, a, b, n=500):
    h = float(b-a) / n
    hs = np.array([2*j*h for j in range(1,n//2+1)])
    # print(a+hs[:-1], a+hs-h)
    return h/3 * (f(a) + 2*sum(f(a+hs[:-1])) + 4*sum(f(a+hs-h)) + f(b))

def func(x):
    return square(x)

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

F, (plt1) = plt.subplots(1,1,figsize=(15,4))
plt1.set(title="Fourier coefficients of a square function")
plt1.plot(np.arange(0,26,0.01), 4/(np.pi*np.arange(0,26,0.01)), "--", c="grey")
plt1.text(s=r"$\frac{4}{\pi\cdot k}$", x=0,y=1.83, fontsize="x-large", color="grey")

plt1.plot(np.arange(maxk+1), 
         FourierSeries(func, 2*np.pi)[0],'ro', alpha=0.5, label=r"$a_k$")
plt1.plot(np.arange(maxk+1), 
         FourierSeries(func, 2*np.pi)[1],'bo', alpha=0.5, label=r"$b_k$")
plt1.legend()
plt1.grid()
plt1.set(xticks=range(26), xlabel="k", ylabel="value", ylim=[-0.2,2],xlim=[-0.5,25.5])
