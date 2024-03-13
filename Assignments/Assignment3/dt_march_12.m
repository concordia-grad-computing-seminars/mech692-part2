function [P2] = dt_march_12(Pin,t,xi,dt,dx,c,npts,nghost,nsteps,BC_update)
  P2 = 0*Pin;
  P1 = Pin;
  for j = 1:nsteps
    P2(nghost+1:npts+nghost) = P1(nghost+1:npts+nghost) - 0.5*(c*dt/dx)*(P1(nghost+2:npts+nghost+1) - P1(nghost:npts+nghost-1));
    t = t + dt;
    P2 = BC_update(P2,nghost,0);
    P1 = P2;
  end
end
  