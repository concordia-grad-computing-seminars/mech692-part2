import numpy as np

def grid(xL,xR,npts,nghost):
    dx  = (xR-xL)/npts
    xi  = np.linspace(xL-nghost*dx+dx/2,xR+nghost*dx-dx/2,npts+2*nghost)
    ind = np.arange(-nghost,npts+nghost,1) 
    return dx, xi, ind

def BC_constant(Pin, nghost, value):
    Pout = np.copy(Pin)
    for i in np.arange(0,nghost,1):
        Pout[i]    = value
        Pout[-i-1] = value
    return Pout

def BC_periodic(Pin, nghost, value):
    Pout = np.copy(Pin)
    for i in np.arange(0,nghost,1):
        Pout[i]    = Pout[-2*nghost + i]
        Pout[-i-1] = Pout[ 2*nghost -1 - i]
    return Pout

def dt_march_11(Pin, t, xi, dt, dx, c, npts, nghost, nsteps, BC_update):
    P2 = 0*Pin
    P1 = np.copy(Pin)
    for j in np.arange(1,nsteps+1,1):
        P2[nghost:npts+nghost:1] = P1[nghost:npts+nghost:1] - (c*dt/dx)*(P1[nghost:npts+nghost:1] - P1[nghost-1:npts+nghost-1:1])
        t = t + dt
        P2 = BC_update(P2,nghost,0)
        P1 = P2
    return P2

def dt_march_12(Pin, t, xi, dt, dx, c, npts, nghost, nsteps, BC_update):
    P2 = 0*Pin
    P1 = np.copy(Pin)
    for j in np.arange(1,nsteps+1,1):
        P2[nghost:npts+nghost:1] = P1[nghost:npts+nghost:1] - 0.5*(c*dt/dx)*(P1[nghost+1:npts+nghost+1:1] - P1[nghost-1:npts+nghost-1:1])
        t = t + dt
        P2 = BC_update(P2,nghost,0)
        P1 = P2
    return P2


def dt_march_22(Pin, t, xi, dt, dx, c, npts, nghost, nsteps, BC_update):
    P2 = 0*Pin
    P1 = np.copy(Pin)
    for j in np.arange(1,nsteps+1,1):
        P2[nghost:npts+nghost:1] = P1[nghost:npts+nghost:1]-0.5*(c*dt/dx)*(P1[nghost+1:npts+nghost+1:1]-P1[nghost-1:npts+nghost-1:1])+0.5*(c*dt/dx)**2*(P1[nghost+1:npts+nghost+1:1]-2*P1[nghost:npts+nghost:1]+P1[nghost-1:npts+nghost-1:1])
        t = t + dt;
        P2 = BC_update(P2,nghost,0);
        P1 = P2;
    return P2

def CFDApp(integrator, CFL, npts):
    xL = 0
    xR = 5
    nghost = 1
    c = 1
    L = xR-xL
    BC_update = BC_periodic   #this is how we choose the boundary condition
    ##BC_update = BC_constant  #this would select the constant boundary condition

    dx,xi, ind = grid(xL,xR,npts,nghost)  #Generate the discretization a.k.a. the grid
    P0 = np.sin(2*np.pi*xi/(xR-xL))                #Calculate the initial condition
    P0 = BC_update(P0,nghost,0)
    t  = 0

    dt = CFL*dx/c

    nsteps = L/(c*dt)

    Pend = integrator(P0,t,xi,dt,dx,c,npts,nghost,nsteps,BC_update);

    AbsErr = max(abs(Pend - P0))
    
    return dx, dt, xi, P0, Pend, AbsErr