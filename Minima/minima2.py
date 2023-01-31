import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib import style

plt.style.use('default')

x = np.arange(0,11,0.2)

def f(x):
    return (x-5)**2-4
def fderiv(x):
    return 2*(x-5)

ff, (plt1,plt2) = plt.subplots(2,1,sharex=True,sharey=True)
plt1.plot(x,f(x), label="$f\ (x)$")
plt1.plot(x,np.zeros(len(x)),'black')
plt1.set(xticks=np.arange(0,11,1))
plt.xlim(0,10)
plt.ylim(-10,25)
plt1.grid()
plt1.legend(loc='upper right')

plt2 = plt.subplot(212)
plt2.plot(np.zeros(len(x)),'black')
plt2.plot(x, fderiv(x), 'orange', label="$\\frac{df(x)}{dx}$")
plt2.set(xticks=np.arange(0,11,1))
plt.xlim(0,10)
plt.grid()
plt2.legend()

def run(plotting, x0, tol):
    x1 = x0
    nsteps=0
    while abs(f(x1)) > tol:
        # if(plotting):
        #     plt1.plot(x1, f(x1), color=(abs(math.atan(nsteps))/1.6,0,0),
        #          marker='s')

        nsteps += 1
        x1 = x1 - f(x1)/fderiv(x1)

    if(plotting):
        print(x1,f(x1))
    
    return nsteps
        
run(True, 1, 0.0001)
#run(True, 5.7, 0.001)

x = np.arange(1e-4,5,1e-4)
#y = [run(False, 1, i) for i in x]

# plt.figure(2)
# plt.plot(x,y)
# plt.xscale('log')

# x = np.arange(-50,210,1)
# plt.figure(3)
# plt.plot(x,f(x), x,np.zeros(len(x)),'black')
# run(True,5.01,0.01)


# plt3 = plt.subplot(111)

# x = np.arange(1e-4,3,1e-4)
# y = [run(False, 5.01, i) for i in x]

# plt3.plot(x,y, label='$x_0=5.01$')


# x = np.arange(1e-4,3,1e-4)
# y = [run(False, 1, i) for i in x]

# plt3.plot(x,y, label='$x_0=1$')
# plt3.set(xscale='log', xlabel="tolerance", ylabel="n_steps", ylim=(0,12),
#          yticks=np.arange(0,12,1))

# plt3.grid()
# plt3.legend()


plt.show()