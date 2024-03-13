function [Pout] = BC_periodic(Pin,nghost, value)
  Pout = Pin;
  Pout(1:nghost)         = value;
  Pout(end-nghost+1:end) = value;
end