import numpy as np
import matplotlib.pyplot as plt

k = 0.5
A = 0.0
phi = 0.66667

def graphs_init():
    F, (plt1, plt2) = plt.subplots(1,2,sharey=True, figsize=(16,5))
    F2, plt3 = plt.subplots(figsize=(14,5))
    #F.tight_layout()
    
    plt1.set(xlim=[0,xmax], ylim=[-np.pi-0.1,np.pi+0.1], title="Theta vs. t")
    plt2.set(xlim=[0,xmax], ylim=[-np.pi-0.1,np.pi+0.1], title="Omega vs. t")
    plt3.set(xlim=[0,xmax], ylim=[-np.pi-0.1,np.pi+0.1], 
              title="Trapezoid vs. Runge-Kutta")
    plt3.set(yticks=[-np.pi,-np.pi/2, 0, np.pi/2, np.pi], xlabel=r"$t$")
    plt3.set_yticklabels([r"$-\pi$",r"$-\frac{\pi}{2}$",
                           "0",r"$\frac{\pi}{2}$",r"$\pi$"], fontsize=14)

    plt2.grid()
    plt1.set(yticks=[-np.pi,-np.pi/2, 0, np.pi/2, np.pi], xlabel=r"$t$")
    plt1.set_yticklabels([r"$-\pi$",r"$-\frac{\pi}{2}$",
                          "0",r"$\frac{\pi}{2}$",r"$\pi$"], fontsize=14)
    plt1.grid()
    plt3.grid()
    plt2.set(xlabel=r"$t$")
    plt3.set(xlabel=r"$t$")
    return F, plt1, plt2, plt3

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
    if (np.abs(theta) > np.pi):
        theta -= 2 * np.pi * np.abs(theta) / theta
    return theta, omega, t


t = 0.0
dt = 0.01
nsteps = 0
xmax = 50

for theta0, omega0 in [[0.2, 0.0], [1.0, 0.0], [3.1415, 0.0], [0.0, 1.0],
                      [3.0,0.0]]:
    t=0
    
    F, plt1, plt2, plt3 = graphs_init()
    F.suptitle(r"$\theta_0=$"+str(theta0)+r",    $\omega_0=$" + str(omega0))
    plt3.set(title=r"Trapezoid vs. Runge-Kutta ($\theta_0=$"+str(theta0)+
              r",  $\omega_0=$"+str(omega0)+")")
    plt3.set(title="Theta vs. time for two different methods"+
             r" ($\theta_0=3.1415,\ \omega_0=0.0$)")
    
    theta, omega, thetas, omegas = theta0, omega0, [theta0], [omega0]
    for i in range(int(xmax//dt)):
        theta, omega, t = step(theta,omega,t,dt, "linear trap")
        thetas.append(theta)
        omegas.append(omega)
        
    # plt1.plot(np.arange(0,xmax,dt), thetas, c=[1,0.4,0], label=r"linear")
    # plt2.plot(np.arange(0,xmax,dt), omegas, c=[1,0.7,0.1], label=r"linear")
    
    # t=0    
    # theta, omega, thetas, omegas = theta0, omega0, [theta0], [omega0]
    # for i in range(int(xmax//dt)):
    #     theta, omega, t = step(theta,omega,t,dt, "non-linear trap")
    #     thetas.append(theta)
    #     omegas.append(omega)
        
    # plt1.plot(np.arange(0,xmax,dt), thetas, c=[0,0.4,1], label=r"nonlinear")
    # plt2.plot(np.arange(0,xmax,dt), omegas, c=[0.1,0.7,1], label=r"nonlinear")
    # plt3.plot(np.arange(0,xmax,dt), thetas, c=[1,0.7,0.1], label=r"$\theta$"+
    #           "(Trapezoid)")
    
    t=0    
    theta, omega, thetas, omegas = theta0, omega0, [theta0], [omega0]
    for i in range(int(xmax//dt)):
        theta, omega, t = step(theta,omega,t,dt, "runge-kutta")
        thetas.append(theta)
        omegas.append(omega)
        
    plt1.plot(np.arange(0,xmax,dt),thetas,c=[0,0.4,1], label=r"$\theta$")
    plt2.plot(np.arange(0,xmax,dt),omegas,c=[1,0.7,0.1], label=r"$\omega$")

    
    plt1.legend(loc="upper right")
    plt2.legend(loc="upper right")
    plt3.legend(loc="upper right")

plt.show()