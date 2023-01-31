import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

g, b, dt, tmax, figi = 9.81, 0.5, 0.01, 15.0, 0
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

def gen_x(m, x0, vx0, t):
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

masses = np.arange(0.1, 5, 0.2)
angles = np.arange(25, 45.0, 0.2)*np.pi/180

best_dist = np.zeros(len(masses))
best_theta = np.zeros(len(masses))

for i in range(len(masses)):
    m = masses[i]
    print("calculating: m="+str(m))
    for theta in angles:
        
        vx0 = V*np.cos(theta)
        vy0 = V*np.sin(theta)
        
        y = gen_y(m, y0, vy0, t)
        x = gen_x(m, x0, vx0, t) 
        
        max_i = 0
        for j in range(len(y)):
            if y[j]<0 and y[j-1]>0:
                max_i = j
                break
            
        distance = y[max_i-1]/(y[max_i-1]-y[max_i])*x[max_i]-y[max_i]/(y[max_i-1]-y[max_i])*x[max_i-1]
        # print(m, theta*180/np.pi, distance)

        if distance > best_dist[i]:
            best_dist[i] = distance
            best_theta[i] = theta*180/np.pi

plt.plot(masses, best_theta)
plt.xlim([0,masses[-1]])
plt.plot([0,100],[45,45])
