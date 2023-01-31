import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

g, dt, tmax, figi = 9.81, 0.01, 15.0, 0
plt.rcParams['figure.dpi'] = 200

def gen_v_quadr(m, b, c, vx0, vy0, t):
    vx, vy = vx0, vy0
    Vx=[vx0]
    Vy=[vy0]
    for tt in t:
        vx += -c/m*np.sqrt(vx**2+vy**2)*vx*dt
        vy += -g*dt-c/m*np.sqrt(vx**2+vy**2)*vy*dt
        Vx.append(vx)
        Vy.append(vy)
    return Vx, Vy

def gen_xy(m, b, c, x0, y0, vx0, vy0, t):
    ############### y
    vy = np.zeros(len(t))
    vy[0] = vy0
    
    for i in range(len(t)-1):
        vy[i+1] = vy[i] - g*dt - b/m * vy[i] * dt
   
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(len(t)-1):
        y[i+1] = y[i] + vy[i]*dt
        
    ############### x
    vx = np.zeros(len(t))
    vx[0] = vx0
    
    for i in range(len(t)-1):
        vx[i+1] = vx[i] - b/m * vx[i] * dt
   
    x = np.zeros(len(t))
    x[0] = x0
    for i in range(len(t)-1):
        x[i+1] = x[i] + vx[i]*dt
        
    return x,y
        
def gen_xy_quadr(m, b, c, x0, y0, vx0, vy0, t):
    vx, vy = gen_v_quadr(m, b, c, vx0, vy0, t)
    
    y = np.zeros(len(t))
    y[0] = y0
    x = np.zeros(len(t))
    x[0] = x0
    
    for i in range(len(t)-1):
        y[i+1] = y[i] + vy[i]*dt
    for i in range(len(t)-1):
        x[i+1] = x[i] + vx[i]*dt
    return x,y

t = np.arange(0,tmax,dt)

x0, y0 = 0, 0
vx0, vy0 = 5,5
b = 0.5
c=0.5
m = 3

x,y = gen_xy(m,b,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="linear air resistance")
x,y = gen_xy_quadr(m,b,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="quadratic air resistance")
x,y = gen_xy(m,0,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="no air resistance")

plt.ylim(bottom=0, top=1.2*max(y))
plt.xlim([0, 6])
plt.legend()
#################################################

plt.figure(2)
x0, y0 = 0, 0
vx0, vy0 = 10,5
b = 0.5
c=0.5
m = 3

x,y = gen_xy(m,b,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="linear air resistance")
x,y = gen_xy_quadr(m,b,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="quadratic air resistance")
x,y = gen_xy(m,0,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="no air resistance")

plt.ylim(bottom=0, top=1.2*max(y))
plt.xlim([0, 12])
plt.legend()
#################################################

plt.figure(3)
x0, y0 = 0, 0
vx0, vy0 = 5,5
b = 0.5
c = 0.5
m = 10

x,y = gen_xy(m,b,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="linear air resistance")
x,y = gen_xy_quadr(m,b,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="quadratic air resistance")
x,y = gen_xy(m,0,c,x0,y0,vx0,vy0,t)
plt.plot(x,y, label="no air resistance")

plt.ylim(bottom=0, top=1.2*max(y))
plt.xlim([0, 6])
plt.legend()
