syms T rho1 rho2 rho3 rho4 rho5

kinetics_params = kinetics;
[species_thermo_structs, Mw, RU] = thermodynamicsModels;
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

rhos = [rho1 rho2 rho3 rho4 rho5];
rho_tilde = rhos ./ Mw;

pressureTerm = log(P_atm/RU) - fT(8);
dpressureTerm = -dfT(8);
Gformation = sym(zeros(1,ns));
dGformation = sym(zeros(1,ns));
for i = 1:ns
  c1 = nasa9_Gcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
  c2 = nasa9_Gcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
  c3 = nasa9_Gcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);
  h1 = sum(c1.*fT);
  h2 = sum(c2.*fT);
  h3 = sum(c3.*fT);  
  Gformation(i) = fsw(1)*h1 + fsw(2)*h2 + fsw(3)*h3 - pressureTerm;
  
  dh1 = sum(c1.*dfT);
  dh2 = sum(c2.*dfT);
  dh3 = sum(c3.*dfT);  
  dGformation(i) = fsw(1)*dh1 + fsw(2)*dh2 + fsw(3)*dh3 + dfsw(1)*h1 + dfsw(2)*h2 + dfsw(3)*h3 - dpressureTerm;
end

lnkf_r = sym(zeros(1,nr));
dlnkf_r = sym(zeros(1,nr));
for ir = 1:nr  
  lnkf_r(ir) = log(A_r(ir)) + beta_r(ir) * fT(8) - theta_r(ir)*fT(6);
  dlnkf_r(ir) = beta_r(ir) * dfT(8) - theta_r(ir)*dfT(6);
end

nu_bf = zeros(ns, nr);
for ir = 1:nr
  for k = 1:ns
    nu_bf(k, ir) = (nu_b_kj(k,ir) - nu_f_kj(k,ir)) ;    
  end
end

lnkb_r = sym(zeros(1,nr));
dlnkb_r = sym(zeros(1,nr));
for ir = 1:nr
  lnkb_r(ir) = lnkf_r(ir);
  dlnkb_r(ir) = dlnkf_r(ir);
  for k = 1:ns
    lnkb_r(ir) = lnkb_r(ir) + nu_bf(k,ir) * Gformation(k);
    dlnkb_r(ir) = dlnkb_r(ir) + nu_bf(k,ir) * dGformation(k);
  end
end
% lnkb_r = Gformation * nu_bf = (ng x ns) * (ns x nr) = ng x nr 
% dlnkb_r = dGformation * nu_bf = (ng x ns) * (ns x nr) = ng x nr 

kf_r = exp(lnkf_r);
dkf_r = kf_r.*dlnkf_r;

kb_r = exp(lnkb_r);
dkb_r = kb_r.*dlnkb_r;

Alpha_rs = zeros(nr,ns+1);
Alpha_rs(1:3,2:(ns+1)) = alpha_jr';
Alpha_rs(4:nr,1) = 1;
thirdbody_r = (Alpha_rs)*[sym(1); rho_tilde(:)];
dthirdbody_r = Alpha_rs(:,2:end);
for k = 1:ns
  dthirdbody_r(:,k) = dthirdbody_r(:,k)/Mw(k);
end

Rr = sym(zeros(1,nr));
dRrdT = sym(zeros(1,nr));
dRrdrho = sym(zeros(nr,ns));
for r = 1:nr
  Cf = 1;
  Cb = 1;        
  for s = 1:ns
    Cf = Cf * rho_tilde(s)^nu_f_kj(s,r);
    Cb = Cb * rho_tilde(s)^nu_b_kj(s,r);
  end
  dCf = sym(zeros(1,ns));
  dCb = sym(zeros(1,ns));
  for s = 1:ns
    dCf(s) = (nu_f_kj(s,r)/Mw(s))*(Cf/rho_tilde(s));
    dCb(s) = (nu_b_kj(s,r)/Mw(s))*(Cb/rho_tilde(s));
  end

  Rr(r) = (kf_r(r) * Cf - kb_r(r) * Cb) * thirdbody_r(r);
  dRrdT(r) = (dkf_r(r) * Cf - dkb_r(r) * Cb) * thirdbody_r(r);
  for s = 1:ns
    dRrdrho(r,s) = (kf_r(r) * dCf(s) - kb_r(r) * dCb(s)) * thirdbody_r(r) + (kf_r(r) * Cf - kb_r(r) * Cb) * dthirdbody_r(r,s);
  end        
end
% Rr = (ng x nr) 
% dRrdT * nuWw =  (ng x nr) 
% dRrdrho * nuWw =  (ng x ns * nr) 

nuMw = zeros(ns,nr);
for k = 1:ns
  for r = 1:nr
    nuMw(k,r) = Mw(k) * (nu_b_kj(k,r) - nu_f_kj(k,r)) ;
  end    
end

omega = sym(zeros(1,ns));
domegadT = sym(zeros(1,ns));
domegadrho = sym(zeros(ns,ns));
for k = 1:ns
    tmp = 0;
    dtmp = 0;    
    drho = sym(zeros(1,ns));
    for r = 1:nr
        tmp = tmp + nuMw(k,r) * Rr(r);
        dtmp = dtmp + nuMw(k,r) * dRrdT(r);        
        for j = 1:ns
          drho(j) = drho(j) + nuMw(k,r) * dRrdrho(r,j);        
        end
    end    
    
    omega(k) = tmp;
    domegadT(k) = dtmp;
    domegadrho(k,:) = drho;
end

% omega = Rr * nuWw =  (ng x nr) * (nr x ns)  = ng x ns 
% domegadT = dRrdT * nuWw =  (ng x nr) * (nr x ns)  = ng x ns 
% domegadrho = dRrdrho * nuWw =  (ng x ns * nr) * (nr x ns)  = ng x ns x ns






