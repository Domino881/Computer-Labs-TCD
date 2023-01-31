import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import square
from datetime import datetime
from matplotlib.ticker import MultipleLocator

startTime = datetime.now()
plt.rcParams['figure.dpi']=200

def func(x):
    return np.sin(0.45*np.pi*x)#*np.exp(-0.2*x)
    # return (0 < x)&(x < 5)

def Nh(N, h=0, ran=0):
    if h>0:
        return N, h
    elif ran:
        return N, ran/N
    else:
        return False

def FourierSeries(f):
    N = len(f)
    omega = 2*np.pi/N

    mn = np.mgrid[0:N,0:N][0]
    mn = np.multiply(mn, mn.T)
    e_vector = np.exp(-1j*mn*omega)

    F = np.dot(e_vector,f)
    return F

def FourierSeriesFFT(f):
    N = len(f)
    if N%2 != 0:
        return FourierSeries(f)
    else: 
        F_even = FourierSeriesFFT(f[::2])
        F_odd = FourierSeriesFFT(f[1::2])
    
    terms = np.exp(-2j * np.pi * np.arange(N) / N)
    return np.concatenate([F_even + terms[:N//2] * F_odd,
                           F_even + terms[N//2:] * F_odd])


def Fourier(f, samples):
    F_n = FourierSeriesFFT(f)
    omega = 2*np.pi/samples
    mn= np.mgrid[0:samples,0:samples][0]
    mn = np.multiply(mn,mn.T)

    return 1/samples*np.dot(np.exp(1j*omega*mn),F_n)
    # return f
    
i=0
    
for (NN, hh, rann) in [(128, 0.1, 0), (128, 0, 2/0.45)]:
    N, h = Nh(NN , hh, rann)
    print(N,h)
    x = np.arange(0,N*h,h)
    y = func(x)
    
    plt.figure(i)
    plt.plot(x,y, 's', c='C1', markersize=2, label="sampled points")
    plt.plot(np.arange(-1,30,0.01), func(np.arange(-1,30,0.01)), 'C1', alpha=0.4)
    plt.xlim([-0.5,2/0.45+0.5])
    yy = np.real(Fourier(y, N))
    # plt.plot(x,yy,'o',c='blue', markersize=2, label="calculated points")
    plt.legend(loc="upper right",framealpha=1)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("f(t)")
    
    plt.figure(i+1)
    plt2 = plt.subplot()
    a = FourierSeries(y)
    
    ar = abs(np.real(a))<abs(np.imag(a))
    plt2.stem(np.arange(N)+0.1, ~ar * np.real(a), label=r"$F_{real}$",
              linefmt='r', markerfmt=" ",basefmt="grey")
    plt2.stem(np.arange(N), np.imag(a), label=r"$F_{imaginary}$", 
              linefmt="b", markerfmt=" ",basefmt="grey")
    plt2.stem(np.arange(N)+0.1, ar * np.real(a),
              linefmt='r', markerfmt=" ",basefmt="grey")
    plt2.set(xlim=[0,N], xticks=range(0,N,10), xlabel="n")
    plt2.xaxis.set_minor_locator(MultipleLocator(1))
    plt2.legend()
    i+=2
    
plt.figure(i)
N, h = Nh(128 , 0.1, 0)
x = np.arange(0,N*h,h)
plt.plot(np.arange(-1,N*h+1,0.01), func(np.arange(-1,N*h+1,0.01)), 'C1', alpha=0.4)
plt.plot(x, Fourier(func(x), N), 'C0o', markersize=2, label="reconstructed points")
plt.xlim([-0.5,N*h+0.5])
plt.grid()
plt.legend()

print("time of execution: " + str(datetime.now() - startTime))

# plt.plot(range(maxk+1), b, 'bo')