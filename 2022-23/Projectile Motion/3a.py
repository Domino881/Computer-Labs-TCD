import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

g, b, dt, tmax, figi = 9.81, 0.5, 0.001, 15.0, 0
plt.rcParams['figure.dpi'] = 200

def gen_y(m, y0, vy0, t):
    vy = np.zeros(len(t))
    vy[0] = vy0
    
    for i in range(len(t)-1):
        vy[i+1] = vy[i] - g*dt - b/m * vy[i] * dt
   
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(len(t)-1):
        y[i+1] = y[i] + vy[i]*dt
    return y

def gen_x(m, y0, vy0, t):
    vx = np.zeros(len(t))
    vx[0] = vx0
    
    for i in range(len(t)-1):
        vx[i+1] = vx[i] - b/m * vx[i] * dt
   
    x = np.zeros(len(t))
    x[0] = x0
    for i in range(len(t)-1):
        x[i+1] = x[i] + vx[i]*dt
    return x

t = np.arange(0,tmax,dt)

x0, y0 = 0, 0
V = 5

masses = np.arange(0.01, 10, 0.5)
angles = np.arange(25.0, 45.0, 1)

best_dist = np.zeros(len(masses))
best_theta = np.zeros(len(masses))

for m in masses:
    for theta in angles:
        vx0 = V*np.cos(theta)
        vy0 = V*np.sin(theta)
        
        y = gen_y(m, y0, vy0, t)
        x = gen_x(m, x0, vx0, t) 
        
        max_i = 0
        for i in range(len(y)):
            if y[i]<0 and y[i-1]>0:
                max_i = i
                break
        distance = y[max_i-1]/(y[max_i-1]-y[max_i])*x[nax_i]-y[max_i]/(y[max_i-1]-y[max_i])*x[max_i-1]
        if distance > best_dist[m]:
            best_dist[m] = distance
            best_theta[m] = theta


