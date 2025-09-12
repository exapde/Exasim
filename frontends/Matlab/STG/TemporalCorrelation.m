function [uPrime, vPrime, wPrime] = TemporalCorrelation(uPrime, vPrime, wPrime, us, vs, ws, dt, T)

c1 = exp(-pi * dt / (2.0 * T));
c2 = sqrt(1.0 - exp(-pi * dt / T));

uPrime = uPrime * c1 + us * c2;
vPrime = vPrime * c1 + vs * c2;
wPrime = wPrime * c1 + ws * c2;                
