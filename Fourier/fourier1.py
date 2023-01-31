import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import square
plt.rcParams['figure.dpi']=200


maxk=10

def Simpson(f, a, b, n=500):
    h = float(b-a) / n
    hs = np.array([2*j*h for j in range(1,n//2+1)])
    # print(a+hs[:-1], a+hs-h)
    return h/3 * (f(a) + 2*sum(f(a+hs[:-1])) + 4*sum(f(a+hs-h)) + f(b))

def func(x):
    # x = x-2*(x//2)
    # return np.exp(x)
    return square(5/np.pi*x)

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
T = 2/5*np.pi**2

print(Simpson(lambda x: np.exp(x), 0, 1, n=8))
print(np.exp(1)-1)

for (f,lab) in [(lambda x: np.sin(x), r"$sin(t)$"),
                (lambda x: np.cos(x)+3*np.cos(2*x)-4*np.cos(3*x), r"$cos(t)+3cos(2t)-4cos(3t)$"),
                (lambda x: np.sin(x)+3*np.sin(3*x)+5*np.sin(5*x), r"$sin(t)+3sin(3t)+5sin(5t)$"),
                (lambda x: np.sin(x)+2*np.cos(3*x)+3*np.sin(5*x), r"$sin(t)+2cos(3t)+3sin(5t)$")]:
    F, (plt1, plt2) = plt.subplots(1,2,figsize=(15,4))
    
    plt1.plot([-100,100],[0,0], 'grey')
    plt1.plot(np.arange(0,np.pi*4,0.01), f(np.arange(0,np.pi*4,0.01)),
              label=lab)
    plt1.set(xlim=[0,4*np.pi], xlabel="t", ylabel="f(t)")
    plt1.legend()
    
    plt2.plot(np.arange(maxk+1), 
             FourierSeries(f, 2*np.pi)[0],'ro', alpha=0.5, label=r"$a_k$")
    plt2.plot(np.arange(maxk+1), 
             FourierSeries(f, 2*np.pi)[1],'bo', alpha=0.5, label=r"$b_k$")
    plt2.legend()
    plt2.grid()
    plt2.set(xticks=range(11), xlabel="k", ylabel="value")