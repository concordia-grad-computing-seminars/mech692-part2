function [Pout] = BC_periodic(Pin,nghost, value)
  Pout = Pin;
  Pout(1:nghost)         = Pin(end-2*nghost+1:end-nghost);
  Pout(end-nghost+1:end) = Pin(nghost+1:2*nghost);
end