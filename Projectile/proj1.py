import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi']=200

B = 1.6e-4
C = 0.25
D = 0.07

def f(V):
    return B*D*V + C*D*D*V**2

Vrange = np.arange(0,10,0.1)

# plt.figure(1)
# plt.plot(Vrange, f(Vrange))
# plt.xlabel("V")
# plt.ylabel("f(V)")

# F, (plt1,plt2) = plt.subplots(1,2, figsize=(18,5))
# plt1.plot(D*Vrange,B*D*Vrange, label=r"$bV$")
# plt1.set(xlabel=r"$D\times V$")
# plt1.legend(loc="upper right")


# plt2.plot(D*Vrange,C*D*D*(Vrange**2), label=r"$cV^2$")
# plt2.set(xlabel=r"$D\times V$")
# plt2.legend(loc="upper right")

# Vrange=np.arange(0.004,0.5,0.002)
# plt.figure(4)
# plt.xlim([0,D*0.5])
# plt.plot([0,0.5], [1,1], '--', c='grey')
# plt.plot(D*Vrange, C/B*D*Vrange**2, label=r"$c/b\cdot V$")
# plt.xlabel(r"$D\times V$")
# plt.legend()
# # plt.xscale('log')
# plt.yscale('log')


def V_y(m, g, b, dt, tmax):
    t,v=0,0
    vy=[0.0]
    while (t<tmax):
        v += -g*dt - (b/m)*v*dt
        t += dt
        vy.append(v)
    return np.arange(0,tmax,dt), vy[:1+int(tmax//dt)]

def y(t, v_y):
    y=np.zeros(len(t))
    y[0]=5.0
    for i in range(len(t)-1):
        # if y[i] <= 0.0:
        #     y=y[:i+1]
        #     break
        y[i+1]=y[i]+dt*v_y[i]
    return t, y

def anal_y(t, m, b):
    res = 5.0 - m**2*g/(b**2)*np.exp(-b*t/m)-m*g/b*t + m**2*g/b**2
    l, r=0, len(t)
    while l<r-1:
        mid = (l+r)//2
        if res[mid]>0:
            l=mid
        else:
            r=mid
    return r, res

m, g, b, dt, tmax = 0.3 ,9.81, 0.5, 0.1, 10.0
t, Vs = V_y(m, g, b, dt, tmax)
anal_v = m*g/b*(np.exp(-b*t/m)-1)

F, (plt3, plt4) = plt.subplots(1,2, figsize=(18,5))
plt3.plot([-1,15], [-m*g/b,-m*g/b], '--', c='gray')
plt3.annotate(r"$V_T$", (0,-m*g/b+0.1), c='gray')

plt3.set(title=r"Comparison of numerical and analytical solution of $V_y(t)$")
plt3.set(xlabel="Time (s)")
plt3.set(ylabel=r"$V_y$ (m/s)")
# for mm in [0.1, 0.2, 0.3]:
    # t, Vs = V_y(mm, g, b, dt, tmax)
plt3.plot(t ,Vs, label="numerical solution")
plt3.plot(t, anal_v, label="analytical solution")
plt3.set(xlim=[-0.0,3])
plt3.set(ylim=[-m*g/b-1,1])
plt3.legend(loc='upper right')

# plt.figure(6)
plt4.set(title=r"Error in $V_y$ vs. time")
plt4.set(xlabel="Time (s)")
plt4.set(ylabel="Absolute Error")
plt4.set(xlim=[-0.0,3])

plt4.plot(t,np.abs(Vs-anal_v))

# tfalls=[]
# plt.figure(7)
# plt.title("Height vs. time for different masses")
# plt.ylabel("Height (m)")
# plt.xlabel("Time (s)")
# for mm in [0.07,0.1,0.2,1.0,20.0]:
#     tfall, ya = anal_y(t, mm, b)
#     plt.plot(t ,ya, label=r"$m=$"+"{:10.2f} kg".format(mm))
# plt.ylim([0,6])
# plt.xlim([0,5])
# plt.legend()

# plt.figure(8)
# plt.title("Time to reach y=0 vs. mass of the object")
# plt.xlabel("Mass (kg)")
# plt.ylabel("Time (s)")
# for mm in np.logspace(-1.4,0.1,101):
#     tfall, ya = anal_y(t,mm,b)
#     tfalls.append(tfall*dt)
# plt.plot(np.logspace(-1.4,0.1,101),tfalls)
# plt.xscale('log')


# for i in range(20):
#     if i not in [4]:
#         plt.close(i)