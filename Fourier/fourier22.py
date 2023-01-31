import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import square
from datetime import datetime
from matplotlib.ticker import MultipleLocator

startTime = datetime.now()
plt.rcParams['figure.dpi']=200

def func(x):
    return np.sin(6*np.pi*x)#*np.exp(-0.2*x)
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
    
def format_func(x, ja):
    if x-np.floor(x) > 0.01:
        return r"$\nu_n$"
    else:
        return int(x)

    
i=0
for hh in [0.6,0.2,0.1,0.04]:
    N, h = Nh(32, hh, 1)
    print(N,h)
    x = np.arange(0,N*h,h)
    y = func(x)
    
    plt.figure(i)
    F, plts = plt.subplots(2,2,figsize=(20,5))
    F.suptitle("DFT for $f(t)=cos(6\pi t)$ and h="+str(h), 
               fontsize=25, fontweight=200)

    plts[0][0].plot(x,y, 's', c='C1', markersize=2, label="sampled points")
    plts[0][0].plot(np.arange(-1,N*h+1,0.01), func(np.arange(-1,N*h+1,0.01)), 'C1', alpha=0.4)
    plts[0][0].set(xlim=[-0.5,N*h+0.5])
    yy = np.real(Fourier(y, N))
    # plt.plot(x,yy,'o',c='blue', markersize=2, label="calculated points")
    plts[0][0].legend(loc="upper right",framealpha=1)
    plts[0][0].set(xlabel="t")
    
    a = FourierSeriesFFT(y)
    
    ar = abs(np.real(a))<abs(np.imag(a))
    plts[1][0].stem(2*np.pi/(N*h)*np.arange(N)+0.1, ~ar * np.real(a),
                    label=r"$F_{real}$", linefmt='r', markerfmt=" ",basefmt="grey")
    plts[1][0].stem(2*np.pi/(N*h)*np.arange(N), np.imag(a),
                    label=r"$F_{imaginary}$", linefmt="b", markerfmt=" ",basefmt="grey")
    plts[1][0].stem(2*np.pi/(N*h)*np.arange(N)+0.1, ar * np.real(a),
                    linefmt='r', markerfmt=" ",basefmt="grey")
    plts[1][1].plot(2*np.pi/(N*h)*np.arange(N),np.abs(a)**2,
                    c="purple", alpha=0.5, label="power spectrum")
    plts[1][1].legend(framealpha=1)
    plts[1][0].set(xticks=list(plts[1][0].get_xticks()) + [np.pi/h])
    plts[1][0].set(xlim=[0,np.pi/h],xlabel=r"$\omega_n$")
    plts[1][0].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    plts[1][1].set(xticks=list(plts[1][0].get_xticks()) + [np.pi/h])
    plts[1][1].set(xlim=[0,np.pi/h],xlabel=r"$\omega_n$")
    plts[1][1].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    #xticks=np.arange(0,2*np.pi/h,10),
    # plts[1][0].xaxis.set_minor_locator(MultipleLocator(1))
    plts[1][0].legend(framealpha=1)
    gg = 0.97
    plts[1][0].set_facecolor((gg,gg+0.01,gg+0.02))
    plts[1][1].set_facecolor((gg,gg+0.01,gg+0.02))
    
    plts[0][1].plot(x,yy,c='blue', markersize=2, label="reconstructed function")
    plts[0][1].plot(x,yy,'bo', markersize=2)
    plts[0][1].legend(loc="upper right")
    plts[0][1].set(xlabel="t")

    plts[0][1].plot(np.arange(-1,N*h+1,0.01), func(np.arange(-1,N*h+1,0.01)), 'C1', alpha=0.4)
    plts[0][1].set(xlim=[-0.5,N*h+0.5])
    i+=4

