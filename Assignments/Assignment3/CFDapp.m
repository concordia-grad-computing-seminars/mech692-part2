function [dx,dt,AbsErr] = CFDapp(integrator, CFL, npts)
    xL = 0;
    xR = 5;
    %npts = 100
    nghost = 1;
    c = 1;
    L = xR-xL;
    BC_update = @BC_periodic;  %this is how we choose the boundary condition
    %%BC_update = @BC_constant;  %this would select the constant boundary condition

    [dx,ind,xi] = grid(xL,xR,npts,nghost);  %Generate the discretization a.k.a. the grid
    P0 = sin(2*pi*xi/(xR-xL));              %Calculate the initial condition
    P0 = BC_update(P0,nghost,0);
    t  = 0;
    %CFL = 0.8;
    dt = CFL*dx/c;
    %plot(xi,P0,'o--')              %Plot the initial condition as a check

    nsteps = L/(c*dt);

    Pend = integrator(P0,t,xi,dt,dx,c,npts,nghost,nsteps,BC_update);

    AbsErr = max(abs(Pend - P0));
end
    