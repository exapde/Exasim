function [fT, dfT] = elementaryfunctions(T)

logT = log(T);
Tinv = 1/T;
T2 = T * T;
T3 = T2 * T;
T4 = T3 * T;
T2inv = 1/T2;
logTTinv = logT*Tinv;

fT = [sym(1) T T2 T3 T4 Tinv T2inv logT logTTinv];
dfT = [sym(0), sym(1), 2*T, 3*T2, 4*T3, -T2inv, -2*Tinv*T2inv, Tinv, (1 - logT)*T2inv];
 

% nana9 = zeros(length(fT), length(sw));
% for i = 1:length(fT)
%   for j = 1:length(sw)
%     nana9(i,j) = fT(i)*sw(j);
%   end  
% end
% nasa9eval_G = a3 - b2 - (T*a4)/2 + (a2 + b1)/T - a1/(2*T^2) - (T^2*a5)/6 - (T^3*a6)/12 - (T^4*a7)/20 - a3*log(T) + (a2*log(T))/T
% nasa9eval_G = c1 + c2*T  + c3*T2 + c4*T3 + c5*T4 + c6*Tinv + c7*T2inv + c8*logT + c9*logTTinv; 

% lnkf_r_sym = logForwardRateCoefficients(A_r, beta_r, theta_r, nr, T);
% function lnkf_r = logForwardRateCoefficients(A_r, beta_r, theta_r, nr, Tf_r)
%     lnkf_r = zeros(nr,1,class(Tf_r));
%     for ir = 1:nr %TODO: Tf_r should be a vector for 2T
%         lnkf_r(ir) = log(A_r(ir)) + beta_r(ir) * log(Tf_r) - theta_r(ir)/Tf_r;
%     end
% end


