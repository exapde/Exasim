function c = nasa9_Gcoeff(a, b)
% nasa9eval_G = a3 - b2 - (T*a4)/2 + (a2 + b1)/T - a1/(2*T^2) - (T^2*a5)/6 - (T^3*a6)/12 - (T^4*a7)/20 - a3*log(T) + (a2*log(T))/T
% nasa9eval_G = c1 + c2*T  + c3*T2 + c4*T3 + c5*T4 + c6*Tinv + c7*T2inv + c8*logT + c9*logTTinv; 

c1 = a(3) - b(2);
c2 = -a(4)/2;
c3 = -a(5)/6;
c4 = -a(6)/12;
c5 = -a(7)/20;
c6 = a(2) + b(1);
c7 = -a(1)/2;
c8 = -a(3);
c9 = a(2);
c = [c1 c2 c3 c4 c5 c6 c7 c8 c9];



