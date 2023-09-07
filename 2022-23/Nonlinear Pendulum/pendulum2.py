import numpy as np
import matplotlib.pyplot as plt

k = 0.5
A = 0.9
phi = 0.66667
iteration_number = 0
a = 0

def f(theta, omega, t, linear=False):
    if not linear:
        return -np.sin(theta) - k*omega + A*np.cos(phi*t)
    else:
        return -theta - k*omega + A*np.cos(phi*t)

def step(theta, omega, t, dt, alg_type):
    if alg_type == "linear trap":
        k1a = dt * omega
        k1b = dt * f(theta, omega, t, linear=True)
        k2a = dt * (omega + k1b)
        k2b = dt * f(theta + k1a, omega + k1b, t + dt, linear=True)
        theta += (k1a + k2a) / 2
        omega += (k1b + k2b ) / 2
    
    elif alg_type == "non-linear trap":
        k1a = dt * omega
        k1b = dt * f(theta, omega, t, linear=False)
        k2a = dt * (omega + k1b)
        k2b = dt * f(theta + k1a, omega + k1b, t + dt, linear=False)
        theta += (k1a + k2a) / 2
        omega += (k1b + k2b ) / 2
        
    elif alg_type == "runge-kutta":
        k1a = dt * omega
        k1b = dt * f(theta, omega, t)
        k2a = dt * (omega + k1b/2)
        k2b = dt * f(theta + k1a/2, omega + k1b/2, t + dt/2)
        k3a = dt * (omega + k2b/2)
        k3b = dt * f(theta + k2a/2, omega + k2b/2, t + dt/2)
        k4a = dt * (omega + k3b)
        k4b = dt * f(theta + k3a, omega + k3b, t + dt)
        theta = theta + (k1a + 2*k2a + 2*k3a + k4a) / 6
        omega = omega + (k1b + 2*k2b + 2*k3b + k4b) / 6
    
    else:
        return "no"

    t = t + dt
    if theta > np.pi+a:
        theta -= 2 * (np.pi)
    if theta < -np.pi+a:
        theta += 2* (np.pi)
    return theta, omega, t

dt = 0.01
nsteps = 0
xmax = 50

k = 0.5
phi = 0.66667

F, plts = plt.subplots(nrows=2, ncols=3,figsize=(15,8))
thetas, omegas = [],[]

F.suptitle("Phase portraits for different values of A",
           fontweight='bold')

def phase_portrait(A1, nsteps, transient, a1):
    global A
    A = A1
    global a
    a = a1
    theta, omega, t = 3.0, 0.0, 0.0

    for i in range(nsteps):
        theta, omega, t = step(theta,omega,t,dt,"runge-kutta")
        if i>transient:
            thetas.append(theta)
            omegas.append(omega)


phase_portrait(0.9, 10000, 5000, 0)
plts[0][0].set(title=r"$A=0.9, transient=10\ 000$")
plts[0][0].plot(omegas, thetas, 'bo', markersize=0.3)
thetas, omegas = [],[]

phase_portrait(1.07, 20000, 15000, 0.1)
plts[0][1].set(title=r"$A=1.07, transient=20\ 000$")
plts[0][1].plot(omegas, thetas, 'bo', markersize=0.3)
thetas, omegas = [],[]

phase_portrait(1.35, 150000, 120000, 0)
plts[0][2].set(title=r"$A=1.35, transient=120\ 000$")
plts[0][2].plot(omegas, thetas, 'bo', markersize=0.3)
thetas, omegas = [],[]

phase_portrait(1.47, 30000,20000, 0.2)
plts[1][0].set(title=r"$A=1.47, transient=30\ 000$")
plts[1][0].plot(omegas, thetas, 'bo', markersize=0.3)
thetas, omegas = [],[]

phase_portrait(1.5, 45000, 35000, 0.2)
plts[1][1].set(title=r"$A=1.5, transient=45\ 000$")
plts[1][1].plot(omegas, thetas, 'bo', markersize=0.3)
thetas, omegas = [],[]

phase_portrait(0, 35000, 25000, 0.4)
plts[1][2].set(title=r"$A=0.0, transient=35\ 000$")
plts[1][2].plot(omegas, thetas, 'bo', markersize=0.3)
thetas, omegas = [],[]


plt.show()