import matplotlib.pyplot as plt
import math
import numpy as np


x = np.arange(0,10,0.2)
np.set_printoptions(precision=3)

def f(x):
    return (x-5)**2-4

plot = plt.subplot(211)
plot.plot(x,f(x), x,np.zeros(50))
plot.set(xlim=(0,10), xticks=np.arange(0,11,1))
plot.grid(True)

def run(plotting, tol, x1, x3):
    nsteps=0
    if(f(x3) > f(x1)):
        x1, x3 = x3, x1
        
    x2 = (x1+x3)/2.
    
    while (np.abs(f(x2)) > tol):
        if(plotting):
            plt.plot(x2,f(x2),
                     color=(abs(math.atan(nsteps))/1.6,0,0),marker='s')
            #print(x2,f(x2))

        if (f(x2) > 0):
            x1 = x2
        else:
            x3 = x2
            
        x2 = (x1+x3)/2
        nsteps += 1
        print('%0.3f' % x2, '& %0.5f' % f(x2))

    return nsteps

plot = plt.subplot(211)
plot.plot(x,f(x), '#1f77b4', x,np.zeros(50), 'orange')
plot.set(xlim=(0,10), xticks=np.arange(0,11,1))
plot.grid(True)
run(True, 0.0001, 1., 6.)

# plot = plt.subplot(212)
# plot.plot(x,f(x), '#1f77b4', x,np.zeros(50), 'orange')
# plot.set(xlim=(0,10), xticks=np.arange(0,11,1))
# plot.grid(True)
run(True, 1e-50, 6., 9.)

# x = np.arange(1e-4,5,1e-4)
# y = [run(False, i, 1., 6.) for i in x]

# plt1 = plt.subplot()

# plt1.plot(x,y)
# plt1.set(xscale='log', xlabel='tolerance', ylabel='n_steps')
# plt.grid(True)

plt.show()
