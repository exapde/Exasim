syms T rho1 rho2 rho3 rho4 rho5

kinetics_params = kinetics;
[species_thermo_structs, Mw, RU] = thermodynamicsModels;
[blottner_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport();
Mw = Mw';

Tsw1 = species_thermo_structs{1}.T1;
Tsw2 = species_thermo_structs{1}.T2;

A_r = kinetics_params.A_r;
beta_r = kinetics_params.beta_r;
theta_r = kinetics_params.theta_r;
nu_f_kj = kinetics_params.nu_f_kj;
nu_b_kj = kinetics_params.nu_b_kj;
alpha_jr = kinetics_params.alpha_jr;
T_ref = kinetics_params.T_ref;
P_atm = kinetics_params.P_atm;
nr = kinetics_params.nr;
ns = kinetics_params.ns;

alpha = 1e8;
[fT, dfT] = elementaryfunctions(T);
% fT = [sym(1) T T2 T3 T4 Tinv T2inv logT logTTinv];
% dfT = [sym(0), sym(1), 2*T, 3*T2, 4*T3, -T2inv, -2*Tinv*T2inv, Tinv, (1 - logT)*T2inv];
[fsw, dfsw] = switchfunctions(T, Tsw1, Tsw2, alpha);

rhos = [rho1 rho2 rho3 rho4 rho5]; % rho_i_dim
rho_tilde = rhos ./ Mw;
P = T * sum(rhos ./ Mw) * RU; % p_dim

rho_tsum = sum(rho_tilde);
Xs = rho_tilde/rho_tsum; % X
dXs = sym(zeros(ns, ns));
for i = 1:ns
  for j = 1:ns
    if (i==j)
      t1 = 1/Mw(i);
    else
      t1 = 0;
    end
     dXs(i, j) = (t1*rho_tsum - (1/Mw(j))*rho_tilde(i))/(rho_tsum*rho_tsum);
  end
end

Dij = sym(zeros(ns, ns));
dDijdT = sym(zeros(ns, ns));
dDijdP = sym(zeros(ns, ns));
lnT = fT(8);
dlnT = dfT(8);
Tinv = fT(6);
for i = 1:ns
  for j = 1:ns    
    params = gupta_structs{i,j};
    expD = exp(params.D);
    t1 = params.A * lnT*lnT + params.B * lnT + params.C;    
    dt1 = 2*params.A * lnT * dlnT + params.B * dlnT;    
    Tterm = T^t1;                            % Tterm in evalCurveFit_guptaDij
    dTterm = Tterm*(t1 * Tinv + dt1 * lnT ); %t1*Tterm/T + dt1 * Tterm * lnT
    Dij(i,j) = (expD * Tterm / P) * 10.1325; % Dij in averageDiffusionCoeffs   
    dDijdT(i,j) = (expD * dTterm / P) * 10.1325;    
    dDijdP(i,j) = -(expD * Tterm / (P*P)) * 10.1325;    
  end
end

Ds = sym(zeros(1, ns));
DsdT = sym(zeros(1, ns));
DsdP = sym(zeros(1, ns));
Dsdrho = sym(zeros(ns, ns));
for i = 1:ns
  denom = sym(0); % denom in averageDiffusionCoeffs
  denomdT = sym(0);
  denomdP = sym(0);
  denomdrho = sym(zeros(1,ns));
  for j = 1:ns
    if i ~= j
      denom = denom + Xs(j) / Dij(i,j);
      denomdT = denomdT - Xs(j) * dDijdT(i,j) / (Dij(i,j)*Dij(i,j));
      denomdP = denomdP - Xs(j) * dDijdP(i,j) / (Dij(i,j)*Dij(i,j));
      for k = 1:ns
        denomdrho(k) = denomdrho(k) + dXs(j,k) / Dij(i,j);
      end
    end
  end
  Ds(i) = (1.0 - Xs(i)) / denom; % D_vec
  DsdT(i) = -denomdT*(1.0 - Xs(i)) / (denom*denom);
  DsdP(i) = -denomdP*(1.0 - Xs(i)) / (denom*denom);
  for k = 1:ns
    Dsdrho(i,k) = (-dXs(i,k)*denom - (1.0 - Xs(i))*denomdrho(k)) / (denom*denom);
  end
end    

mus = sym(zeros(1, ns));
musT = sym(zeros(1, ns));
for i = 1:ns
  params = gupta_mu_structs{i};
  t1 = params.A * lnT + params.B;
  dt1 = params.A * dlnT;    
  Tterm = T^t1;
  dTterm = Tterm*(t1 * Tinv + dt1 * lnT ); 
  expC = exp(params.C);
  mus(i) = 0.1 * expC * Tterm; % mu_i
  musT(i) = 0.1 * expC * dTterm;
end

phi = sym(zeros(1,ns));
phiT = sym(zeros(1, ns));
phirho = sym(zeros(ns, ns));
for i = 1:ns
  for j = 1:ns
    if i == j
      phi(i) = phi(i) + Xs(j);      
      for k = 1:ns
        phirho(i,k) = phirho(i,k) + dXs(j,k); 
      end
    else
      mu_ratio = mus(i) / mus(j);
      muT_ratio = (musT(i)*mus(j) - musT(j)*musT(i))/ (mus(j)*mus(j));
      M_ratio = Mw(i) / Mw(j);
      tmp = 1.0 + sqrt(mu_ratio / sqrt(M_ratio));  % is this correct?    
      tmpT = 0.5*muT_ratio/sqrt(mu_ratio / sqrt(M_ratio));
      phi(i) = phi(i) + Xs(j) * tmp^2 / sqrt(8.0 *(1.0 + M_ratio)); % phi_i
      phiT(i) = phiT(i) + Xs(j) * 2 * tmp * tmpT / sqrt(8.0 *(1.0 + M_ratio));
      for k = 1:ns
        phirho(i,k) = phirho(i,k) + dXs(j,k) * tmp^2 / sqrt(8.0 *(1.0 + M_ratio));
      end
    end
  end
end

mud =  sum( mus .* Xs ./ phi); % mu_d_dim
mudT =  sum( (musT .* phi - mus .* phiT) .* Xs ./ (phi.*phi));
mudrho = sym(zeros(1,ns));
for i = 1:ns
  tm = mus(i)/phi(i);
  for j = 1:ns
    mudrho(j) = mudrho(j) + tm * dXs(i,j);
  end
end

cp = sym(zeros(1,ns));
cpT = sym(zeros(1,ns));
for i = 1:ns
  c1 = nasa9_cpcoeff(species_thermo_structs{i}.a1);
  c2 = nasa9_cpcoeff(species_thermo_structs{i}.a2);
  c3 = nasa9_cpcoeff(species_thermo_structs{i}.a3);
  h1 = sum(c1.*fT(1:7));
  h2 = sum(c2.*fT(1:7));
  h3 = sum(c3.*fT(1:7));  
  cp(i) = fsw(1)*h1 + fsw(2)*h2 + fsw(3)*h3;
  cp(i) = cp(i) * RU / Mw(i); % getCpsMass
  
  dh1 = sum(c1.*dfT(1:7));
  dh2 = sum(c2.*dfT(1:7));
  dh3 = sum(c3.*dfT(1:7));  
  cpT(i) = fsw(1)*dh1 + fsw(2)*dh2 + fsw(3)*dh3 + dfsw(1)*h1 + dfsw(2)*h2 + dfsw(3)*h3;
  cpT(i) = cpT(i) * RU / Mw(i);
end

lambda = mud *  (cp + 5/4 * RU./Mw); % lambda_i
lambdaT = mudT * (cp + 5/4 * RU./Mw) + mud * cpT;

kappa =  sum( lambda .* Xs ./ phi); % kappa_dim
kappaT =  sum( (lambdaT .* phi - lambda .* phiT) .* Xs ./ (phi.*phi));
kapparho = sym(zeros(1,ns));
for i = 1:ns
  tm = lambda(i)/phi(i);
  for j = 1:ns
    kapparho(j) = kapparho(j) + tm * dXs(i,j);
  end
end

h = sym(zeros(1,ns));
hT = sym(zeros(1,ns));
for i = 1:ns
  c1 = nasa9_hcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
  c2 = nasa9_hcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
  c3 = nasa9_hcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);
  h1 = sum(c1.*fT);
  h2 = sum(c2.*fT);
  h3 = sum(c3.*fT);  
  h(i) = fsw(1)*h1 + fsw(2)*h2 + fsw(3)*h3;
  h(i) = h(i) * T * RU / Mw(i); % h_vec
  
  dh1 = sum(c1.*dfT);
  dh2 = sum(c2.*dfT);
  dh3 = sum(c3.*dfT);  
  hT(i) = fsw(1)*dh1 + fsw(2)*dh2 + fsw(3)*dh3 + dfsw(1)*h1 + dfsw(2)*h2 + dfsw(3)*h3;
  hT(i) = (hT(i) * T + hT(i)) * RU / Mw(i);
end

rho_sum = sum(rhos);
Ys = rhos/rho_sum; % Y
dYs = sym(zeros(ns, ns));
for i = 1:ns
  for j = 1:ns
    if (i==j)
      t1 = 1;
    else
      t1 = 0;
    end
     dYs(i, j) = (t1*rho_sum - rhos(i))/(rho_sum*rho_sum);
  end
end


cv = (cp - 1.0) * RU ./ Mw; % getCvsMass
cvT = cpT * RU ./ Mw;
cvY = sum(cv.*Ys);          % mixtureFrozenCvMass
cvYT = sum(cvT.*Ys); 
denom = sum(rhos) * cvY;    % denom
denomT = sum(rhos) * cvYT;
denomrho = sym(zeros(1, ns));
for i = 1:ns
  denomrho(i) = cvY;
end

es = (h - 1.0)*T*RU ./ Mw; % e_i
esT = hT*T*RU ./ Mw + (h - 1.0)*RU ./ Mw;

dT_drhos = -es/denom; % dT_drho_i_dim
dT_drhoe = 1.0/denom; % dT_drhoe_dim

