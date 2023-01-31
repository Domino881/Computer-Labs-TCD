import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

g, b, dt, tmax, figi = 9.81, 0.5, 0.001, 15.0, 0
plt.rcParams['figure.dpi'] = 200

def gen_y(m, y0, vy0, t, alg):
    vy = (vy0+m*g/b)*np.exp(-b*t/m) - m*g/b
    if alg=='1':
        y=np.zeros(len(t))
        y[0]=float(y0)
        
        for i in range(len(t)-1):
            y[i+1]=y[i]+dt*vy[i]
            
    elif alg=='2':
        y = -g*m*m*np.exp(-b/m*t)/b**2-g*m*t/b-m*vy0*np.exp(-b/m*t)/b+\
                 g*m*m/b**2+m*vy0/b+y0
    else:
        return False
    
   
    l, r=0, len(t)
    while l<r-1:
        mid = (l+r)//2
        if y[mid]>0:
            l=mid
        else:
            r=mid
    return y[:r], r

    
def gen_x(l, m, x0, vx0, t, alg):
    if vx0 == 0.0:
        return False
    
    vx = vx0*np.exp(-b/m*t)
    
    if alg=='1':
        x=np.zeros(l)
        x[0]=x0
        for i in range(l-1):
            x[i+1]=x[i]+dt*vx[i]
        return x[:l]
    elif alg=='2':
        x = -m*vx0*np.exp(-b/m*t)/b + m*vx0/b + x0
        return x[:l]
    else:
        return False

def anal_y_squared(t, theta, m, V):
    y0 = 0
    vy0 = V*np.sin(theta)
    return (-g*m*m*np.exp(-b/m*t)/b**2-g*m*t/b-m*vy0*np.exp(-b/m*t)/b+\
             g*m*m/b**2+m*vy0/b+y0)**2
def anal_x(x0, vx0, m, V, t):
    return -m*vx0*np.exp(-b/m*t)/b + m*vx0/b + x0

def distance(theta, m, V, t, alg):
    if alg in ['1', '2']:
        theta=theta[0]
        if theta <= 0.0 or theta > np.pi/2:
            return 1000
        y, l = gen_y(m, 0, V*np.sin(theta), t, alg)
        x = gen_x(l, m, 0, V*np.cos(theta), t, alg)
        
        res = y[-2]/(y[-2]-y[-1])*x[-1]-y[-1]/(y[-2]-y[-1])*x[-2]
        return -res
    elif alg=='3':
        tt = optimize.minimize(anal_y_squared, np.pi/4, args=(theta, m, V),
                               options={'maxiter':6})['x'][0]
        return -anal_x(0, V*np.cos(theta), m, V, tt)


def find_angle(m, t, alg):
    print("\rmass is "+str(m), end='', flush=True)
    V = 8.0
    bds = optimize.Bounds(0.1, np.pi/2-0.1)
    return 360/(2*np.pi)*optimize.minimize(distance, np.pi/4, 
                                           args=(m,V,t,alg),
                                           options={'maxiter':15})['x'][0]

t = np.arange(0,tmax,dt)

x0, y0 = 0, 0.0
vx0, vy0 = 5.0, 5.0

y, l = gen_y(3, y0, vy0, t, '1')
x = gen_x(l, 3, x0, vx0, t, '1') 

# print(find_angle3(0.3, t))
# print([distance(np.pi/4+e, 0.3, 5.0, t) for e in np.arange(-0.1,0.1,0.05)])
# plt.axvline(x=45)
# plt.legend()
        
plt.figure(2)
plt.plot(x,y, label="with air resistance")
plt.plot((x0+vx0*t),(y0+vy0*t-0.5*g*t**2), label="in vacuum")
plt.ylim([0,1.6])
plt.xlim([0,6])
plt.xlabel("Distance (m)")
plt.ylabel("Height (m)")
plt.title("Trajectory under air resistance \n as compared to vacuum")
plt.legend()


# plt.figure(2)
# print(find_angle3(10, t))

mms = np.logspace(-1,0.9,15)

plt.figure(3)

angles = [find_angle(mm, t, '2') for mm in mms]
# # print()

# mms2 = mms[::3]
# angles2 = angles[::3]
# plt.xlabel("Distance (m)")
# plt.ylabel("Height (m)")
# plt.title("Optimal trajectories for different mass particles")
# plt.xlim([-0.3,8])
# plt.ylim([0,1.6])
# for i in range(len(mms2)):
#     V=8.0
#     y, l = gen_y(mms2[i], 0.0, V*np.sin(np.pi/180*angles2[i]), t, '1')
#     x = gen_x(l, mms2[i], 0.0, V*np.cos(np.pi/180*angles2[i]), t, '1')
#     plt.plot(x,y, label="m="+"{:.1f}".format(mms2[i]))
# plt.legend()

# plt.figure(4, dpi=300)

plt.grid()
plt.plot(mms, angles)
plt.plot(mms, np.full(np.shape(mms), 45.0))
plt.xlim([0,6.1])
# plt.xlim([0,1])
plt.xlabel("Mass of particle (kg)")
plt.ylabel(r"Optimum angle ($^\circ$)")
plt.title("Optimal launching angle vs. mass of particle")
    