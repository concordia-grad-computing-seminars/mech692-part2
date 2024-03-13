function [dx,ind,xi] = grid(xL,xR,npts,nghost)
  dx  = (xR-xL)/npts;
  ind = 1-nghost:npts+nghost;
  xi  = (ind-1)*dx + dx/2;
end