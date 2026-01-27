function omega = kineticsource(T, rhos)
%syms T rho1 rho2 rho3 rho4 rho5
%rhos = [rho1; rho2; rho3; rho4; rho5];

kinetics_params = kinetics;
[species_thermo_structs, Mw, RU] = thermodynamicsModels;
rho_tilde = rhos ./ Mw;

Tsw1 = species_thermo_structs{1}.T1;
Tsw2 = species_thermo_structs{1}.T2;

A_r = kinetics_params.A_r;
beta_r = kinetics_params.beta_r;
theta_r = kinetics_params.theta_r;
nu_f_kj = kinetics_params.nu_f_kj;
nu_b_kj = kinetics_params.nu_b_kj;
alpha_jr = kinetics_params.alpha_jr;
%T_ref = kinetics_params.T_ref;
P_atm = kinetics_params.P_atm;
nr = kinetics_params.nr;
ns = kinetics_params.ns;

alpha = 1e4;
fT = elementaryfunctions(T);
fsw = switchfunctions(T, Tsw1, Tsw2, alpha);

pressureTerm = log(P_atm/RU) - fT(8);
Gformation = sym(zeros(1,ns));
for i = 1:ns
  c1 = nasa9_Gcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
  c2 = nasa9_Gcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
  c3 = nasa9_Gcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);
  h1 = sum(c1.*fT);
  h2 = sum(c2.*fT);
  h3 = sum(c3.*fT);  
  Gformation(i) = fsw(1)*h1 + fsw(2)*h2 + fsw(3)*h3 - pressureTerm;
end

lnkf_r = sym(zeros(1,nr));
for ir = 1:nr  
  lnkf_r(ir) = log(A_r(ir)) + beta_r(ir) * fT(8) - theta_r(ir)*fT(6);
end

nu_bf = zeros(ns, nr);
for ir = 1:nr
  for k = 1:ns
    nu_bf(k, ir) = (nu_b_kj(k,ir) - nu_f_kj(k,ir)) ;    
  end
end

lnkb_r = sym(zeros(1,nr));
for ir = 1:nr
  lnkb_r(ir) = lnkf_r(ir);
  for k = 1:ns
    lnkb_r(ir) = lnkb_r(ir) + nu_bf(k,ir) * Gformation(k);
  end
end

kf_r = exp(lnkf_r);
kb_r = exp(lnkb_r);

Alpha_rs = zeros(nr,ns+1);
Alpha_rs(1:3,2:(ns+1)) = alpha_jr';
Alpha_rs(4:nr,1) = 1;
thirdbody_r = (Alpha_rs)*[sym(1); rho_tilde(:)];

Rr = sym(zeros(1,nr));
for r = 1:nr
  Cf = 1;
  Cb = 1;        
  for s = 1:ns
    Cf = Cf * rho_tilde(s)^nu_f_kj(s,r);
    Cb = Cb * rho_tilde(s)^nu_b_kj(s,r);
  end
  Rr(r) = (kf_r(r) * Cf - kb_r(r) * Cb) * thirdbody_r(r);
end

nuMw = zeros(ns,nr);
for k = 1:ns
  for r = 1:nr
    nuMw(k,r) = Mw(k) * (nu_b_kj(k,r) - nu_f_kj(k,r)) ;
  end    
end

omega = sym(zeros(ns,1));
for k = 1:ns
    tmp = 0;    
    for r = 1:nr
        tmp = tmp + nuMw(k,r) * Rr(r);
    end        
    omega(k) = tmp;
end






