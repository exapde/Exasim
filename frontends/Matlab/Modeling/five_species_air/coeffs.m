kinetics_params = kinetics;
[species_thermo_structs, Mw, RU] = thermodynamicsModels;
[blottner_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport();

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

Tsw1 = species_thermo_structs{i}.T1;
Tsw2 = species_thermo_structs{i}.T2;
fsw = switchfunctions(T, Tsw1, Tsw2, alpha);      

for i = 1:ns
    hc1(i,:) = nasa9_hcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
    hc2(i,:) = nasa9_hcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
    hc3(i,:) = nasa9_hcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);

    Gc1(i,:) = nasa9_Gcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
    Gc2(i,:) = nasa9_Gcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
    Gc3(i,:) = nasa9_Gcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);
    
    cpc1(i,:) = nasa9_cpcoeff(species_thermo_structs{i}.a1);
    cpc2(i,:) = nasa9_cpcoeff(species_thermo_structs{i}.a2);
    cpc3(i,:) = nasa9_cpcoeff(species_thermo_structs{i}.a3);    
end

nu_bf = zeros(ns, nr);
for ir = 1:nr
  for k = 1:ns
    nu_bf(k, ir) = (nu_b_kj(k,ir) - nu_f_kj(k,ir)) ;    
  end
end

for i = 1:ns
    c1 = nasa9_hcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
    c2 = nasa9_hcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
    c3 = nasa9_hcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);
    h1(i) = sum(c1.*fT);
    h2(i) = sum(c2.*fT);
    h3(i) = sum(c3.*fT);  
    %H(i) = fsw(1)*h1(i) + fsw(2)*h2(i) + fsw(3)*h3(i);      
end

