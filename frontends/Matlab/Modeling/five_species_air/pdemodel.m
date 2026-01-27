function pde = pdemodel
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbouhdg = @fbouhdg;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    % pde.avfield = @avfield;
    pde.sourcew = @sourcew;
    pde.initw = @initw;
    pde.eos = @eos;
    pde.monitor = @monitor;
end

function w0 = initw(x, mu, eta)
    % should not be used, overwrite by mesh.wdg 
    w0 = sym(zeros(1,1));
    w0(1) = 1;
end


function m = mass(u, q, w, v, x, t, mu, eta)
    ns = 5;
    ndim = 2;
    m = sym(ones(ns + ndim + 1, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)
    nch = 8;
    ns = 5;

    % ndim = 2;
    [species_thermo_structs, Mw, ~] = thermodynamicsModels();

    % Nondimensional params
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    mu_scale    = eta(5);
    kappa_scale = eta(6);

    rho_i = sym(zeros(ns,1));
    rho = sym(0);

    % Conservative Variables
    for ispecies = 1:ns
        % rho_i(ispecies) = 0 + lmax(u(ispecies)-0,alphaClip); %subspecies density
        rho_i(ispecies) = u(ispecies);
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(ns+1);
    rhov = u(ns+2);
    rhoE = u(ns+3);

    drho_dx_i = -q(1:ns);
    drhou_dx  = -q(ns+1);
    drhov_dx  = -q(ns+2);
    drhoE_dx  = -q(ns+2+1);
    drho_dy_i = -q((nch+1:nch+ns));
    drhou_dy  = -q(nch+ns+1);
    drhov_dy  = -q(nch+ns+2);
    drhoE_dy  = -q(nch+ns+2+1);
    av = v(1);

    rhoinv = 1.0 / rho;
    uv = rhou * rhoinv; %velocity
    vv = rhov * rhoinv;
    E = rhoE * rhoinv; %energy

    % Mutation outputs
    rho_i_dim = rho_i * rho_scale;
    T = w(1);
    T_dim = T * T_scale;
    p_dim = pressure(T_dim, rho_i_dim, Mw);
    p = p_dim / rhoe_scale;

    H = E + p*rhoinv; %enthalpy

    % Inviscid + AV Fluxes
    for ispecies = 1:ns
        fi(ispecies,1) = rho_i(ispecies) * uv - av.*drho_dx_i(ispecies);
    end
    fi(ns + 1,1) = rhou * uv + p     - av.*drhou_dx;
    fi(ns + 2,1) = rhov * uv         - av.*drhov_dx;
    fi(ns + 3,1) = rhou * H          - av.*drhoE_dx;

    for ispecies = 1:ns
        fi(ispecies,2) = rho_i(ispecies) * vv - av.*drho_dy_i(ispecies);
    end
    fi(ns + 1,2) = rhou * vv      - av.*drhou_dy;
    fi(ns + 2,2) = rhov * vv + p  - av.*drhov_dy;
    fi(ns + 3,2) = rhov * H       - av.*drhoE_dy;

    % Viscous fluxes
    beta = 0;
    Ec          = eta(9);
    Pr          = eta(10);
    Re          = eta(11);

    [dT_drho_i_dim, dT_drhoe_dim, D_vec, h_vec, mu_d_dim, kappa_dim] = transportcoefficients(T_dim, rho_i_dim);

    drho_dx = sum(drho_dx_i);
    drho_dy = sum(drho_dy_i);
    %X = X_i(rho_i_dim,Mw);
    uv = rhou * rhoinv; %velocity
    vv = rhov * rhoinv;
    %E = rhoE .* rhoinv; 
    du_dx = (drhou_dx - drho_dx*uv)*rhoinv;
    dv_dx = (drhov_dx - drho_dx*vv)*rhoinv;
    du_dy = (drhou_dy - drho_dy*uv)*rhoinv;
    dv_dy = (drhov_dy - drho_dy*vv)*rhoinv;
    uTu2      = 0.5 * (uv * uv + vv * vv);
    duTu2_dx  = uv * du_dx + vv * dv_dx; 
    duTu2_dy  = uv * du_dy + vv * dv_dy;
    
    dT_drho_i = dT_drho_i_dim / T_scale * rho_scale;
    dT_drhoe = dT_drhoe_dim / T_scale * rhoe_scale;

    dre_drho  = Ec*-uTu2;
    dre_duTu2 = Ec*-rho;
    dre_drhoE = Ec*1.0;
    dre_dx    = dre_drho * drho_dx + dre_duTu2 * duTu2_dx + dre_drhoE * drhoE_dx;
    dre_dy    = dre_drho * drho_dy + dre_duTu2 * duTu2_dy + dre_drhoE * drhoE_dy;
    dT_dx     = sum(dT_drho_i .* drho_dx_i) +  dT_drhoe * dre_dx;
    dT_dy     = sum(dT_drho_i .* drho_dy_i) +  dT_drhoe * dre_dy;

    h_scale = u_scale^2;
    D_scale = u_scale;

    mu_d = mu_d_dim / mu_scale;
    kappa = kappa_dim / kappa_scale;
    D_vec = D_vec / D_scale;
    h_vec = h_vec / h_scale;

    %%%%%%%% Calculation of J_i
    dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rhoinv * rhoinv;
    dY_dy_i = (drho_dy_i * rho - rho_i * drho_dy) * rhoinv * rhoinv;

    J_i_x = -rho * D_vec .* dY_dx_i + rho_i .* sum(D_vec .* dY_dx_i);
    J_i_y = -rho * D_vec .* dY_dy_i + rho_i .* sum(D_vec .* dY_dy_i);

    %%%%%%%% Stress tensor tau
    txx = mu_d * 2.0/3.0 * (2 * du_dx - dv_dy) / Re + beta * (du_dx + dv_dy);
    txy = mu_d * (du_dy + dv_dx) / Re;
    tyy = mu_d * 2.0/3.0 * (2 * dv_dy - du_dx) / Re + beta * (du_dx + dv_dy);

    % VISCOUS FLUX
    
    for i = 1:ns
        fv(i,1) = -J_i_x(i); 
        fv(i,2) = -J_i_y(i);
    end

    fv(ns + 1, 1) = txx;
    fv(ns + 2, 1) = txy;
    fv(ns + 3,1) = uv * txx + vv * txy - (sum(h_vec.*J_i_x) - kappa*dT_dx / (Re*Pr*Ec));
    
    fv(ns + 1, 2) = txy;
    fv(ns + 2, 2) = tyy;
    fv(ns + 3,2) = uv * txy + vv * tyy - (sum(h_vec.*J_i_y) - kappa*dT_dy / (Re*Pr*Ec));

    f = fi - fv;
end

function s = source(u, q, w, v, x, t, mu, eta)
    
    ns = 5;
    % Nondimensional params
    rho_scale   = eta(1);
    u_scale     = eta(2);
    T_scale     = eta(4);
    L_scale     = eta(8);
    omega_scale = rho_scale * u_scale / L_scale;

    % Mutation outputs
    rho_i_dim = u(1:ns) * rho_scale;
    % rho_i_dim = 0 + lmax(rho_i_dim,alphaClip); %subspecies density
    T = w(1) * T_scale;
    omega_i = kineticsource(T, rho_i_dim);

    s(1:ns) = omega_i / omega_scale;
    s(ns+1) = 0.0;
    s(ns+2) = 0.0;
    s(ns+3) = 0.0;
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    % uinf = initu(x, mu, eta);

    uinf = [0.000000650000000   0.228300980000000   0.010260100000000   0.754307040000000   0.007131230000000   0.996653279785810 0   0.607507866832124];
    uinf = uinf(:);

    f_out = (u - uhat);
    f_in = (uinf - uhat);

    % wall boundary condition    
    un = u(ns+1).*n(1) + u(ns+2).*n(2);  
    ui = u;
    ui(ns+1) = ui(ns+1) - n(1).*un;
    ui(ns+2) = ui(ns+2) - n(2).*un;
    fh = ui - uhat;

    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    ub = 0*[f_in f_in fh fh f_in f_out];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    % uinf = initu(x, mu, eta);

    uinf = [0.000000650000000   0.228300980000000   0.010260100000000   0.754307040000000   0.007131230000000   0.996653279785810 0   0.607507866832124];
    uinf = uinf(:);

    f_out = (u - uhat);
    f_in = (uinf - uhat);

    % wall boundary condition    
    un = u(ns+1).*n(1) + u(ns+2).*n(2);  
    ui = u;
    ui(ns+1) = ui(ns+1) - n(1).*un;
    ui(ns+2) = ui(ns+2) - n(2).*un;
    fh = ui - uhat;

    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    fb = 0*[f_in f_in fh fh f_in f_out];
end

function u0 = initu(x, mu, eta)
    ns = 5;
    rho_scale = eta(1);
    u_scale = eta(2);
    rhoe_scale = eta(3);

    u0(1:ns) = mu(1:ns) / rho_scale;
    u0(ns+1:ns+2) = mu(ns+1:ns+2) / (rho_scale * u_scale);
    u0(ns+2+1) = mu(ns + 2 + 1) / rhoe_scale;
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    mu_scale    = eta(5);
    kappa_scale = eta(6);
    cp_scale    = eta(7);
    L_scale     = eta(8);

    uinf = initu(x, mu, eta);

    % uinf = [0.000000650000000   0.228300980000000   0.010260100000000   0.754307040000000   0.007131230000000   0.996653279785810 0   0.607507866832124];
    uinf = uinf(:);

    f_out = (u - uhat);
    f_in = (uinf - uhat);

    % wall boundary condition    
    un = u(ns+1).*n(1) + u(ns+2).*n(2);  
    ui = u;
    ui(ns+1) = ui(ns+1) - n(1).*un;
    ui(ns+2) = ui(ns+2) - n(2).*un;
    fh = ui - uhat;


    %%% Isothermal wall
    [species_thermo_structs, Mw, RU] = thermodynamicsModels();
    
    rho_i_wall = u(1:ns);
    rho_i_wall_dim = rho_i_wall * rho_scale;
    T_wall = eta(12);

    p_wall = pressure(T_wall, rho_i_wall_dim, Mw);
    X_wall = X_i(rho_i_wall_dim, Mw);
    e_dim = mixtureEnergyMass(T_wall, p_wall, X_wall, Mw, species_thermo_structs);
    rhoE_dim = (sum(rho_i_wall_dim) * e_dim);
    rhoE_wall = rhoE_dim / rhoe_scale;

    uf = u;
    uf(6:7) = 0;
    f_iso = uf - uhat;
    f_iso(8) = rhoE_wall - uhat(8);

    %%% Noncatalytic wall
    % set dY/dn = 0
    % Step 1: get pressure from flow (dP/dn = 0)
    T_flow     = w(1) * T_scale;
    rho_i_flow = u(1:5)*rho_scale;
    p_flow = pressure(T_flow, rho_i_flow, Mw);

    % Step 2: get mass fractions from flow (dY/dn = 0)
    Y_i_flow = Y_i(rho_i_flow);

    % [rho_i_noncat, e_noncat] = density_energy_from_YiTP(Y_i_flow, T_wall, p_flow);
    % rhoE_noncat = sum(rho_i_noncat) * e_noncat;
    X_i_flow = (Y_i_flow .* Mw) / sum(Y_i_flow .* Mw);

    % Step 3: set state with Y, P, T_wall (double check how this is done) - I think really just need density
    rho_noncat = density(T_wall, p_flow, X_i_flow, Mw);
    rho_i_noncat = rho_noncat * Y_i_flow;
    e_dim_noncat = mixtureEnergyMass(T_wall, p_flow, X_i_flow, Mw, species_thermo_structs);
    rhoE_dim_noncat = (sum(rho_i_noncat) * e_dim_noncat);
    rhoE_noncat = rhoE_dim_noncat;
    u_noncat = 0*uhat;
    u_noncat(1:5) = rho_i_noncat / rho_scale;
    u_noncat(6:7) = 0;
    u_noncat(8) = rhoE_noncat / rhoe_scale;
    f_noncat = u_noncat - uhat;

    %%% Noncatalytic wall use uh
    rho_i_flow_uh = uhat(1:5)*rho_scale;
    p_flow_uh = pressure(T_flow, rho_i_flow_uh, Mw);

    % Step 2: get mass fractions from flow (dY/dn = 0)
    Y_i_flow_uh = Y_i(rho_i_flow_uh);
    [rho_i_noncat_uh, e_noncat_uh] = density_energy_from_YiTP(Y_i_flow_uh, T_wall, p_flow_uh);
    rhoE_noncat_uh = sum(rho_i_noncat_uh) * e_noncat_uh;

    u_noncat_uh = 0*uhat;
    u_noncat_uh(1:5) = rho_i_noncat_uh / rho_scale;
    u_noncat_uh(6:7) = 0;
    u_noncat_uh(8) = rhoE_noncat_uh / rhoe_scale;
    f_noncat_uh = u_noncat_uh - uhat;

    %%% Supercatalytic wall: u
    % specify Y_i to Y_eq
    % Step 1: get pressure from flow (dP/dn = 0): computed above
    % Step 2: specify mass fractions
    Y_i_cat = [0; 0; 0; 0.7624; 1.0-0.7624];
    X_i_cat = (Y_i_cat .* Mw) / sum(Y_i_cat .* Mw);
    rho_flow = sum(rho_i_flow);
    rho_i_supercat = rho_flow * Y_i_cat;
    e_dim_supercat = mixtureEnergyMass(T_wall, p_flow, X_i_cat, Mw, species_thermo_structs);
    % rhoE_dim_cat = (sum(rho_i_cat) * e_dim_cat);
    % Step 4: set state with Y, P, T_wall (double check how this is done)
    % rho_cat = density(T_wall, p_flow, X_i_cat, Mw);
    % rho_i_cat = rho_cat * Y_i_cat;
    % e_dim_cat = mixtureEnergyMass(T_wall, p_flow, X_i_cat, Mw, species_thermo_structs);
    % rhoE_dim_cat = (sum(rho_i_cat) * e_dim_cat);
    % [rho_i_supercat, e_supercat] = density_energy_from_YiTP(Y_i_cat, T_wall, p_flow);
    rhoE_cat = sum(rho_i_supercat)*e_dim_supercat;
    u_cat = 0*uhat;
    u_cat(1:5) = rho_i_supercat / rho_scale;
    u_cat(6:7) = 0;
    u_cat(8) = rhoE_cat / rhoe_scale;
    f_cat = u_cat - uhat;

    %%% Supercatalytic wall: uh
    % specify Y_i to Y_eq
    % Step 1: get pressure from flow (dP/dn = 0): computed above
    % Step 2: specify mass fractions
    Y_i_cat = [0; 0; 0; 0.7624; 1.0-0.7624];
    % X_i_cat = (Y_i_cat .* Mw) / sum(Y_i_cat .* Mw);

    % Step 4: set state with Y, P, T_wall (double check how this is done)
    % rho_cat = density(T_wall, p_flow, X_i_cat, Mw);
    % rho_i_cat = rho_cat * Y_i_cat;
    % e_dim_cat = mixtureEnergyMass(T_wall, p_flow, X_i_cat, Mw, species_thermo_structs);
    % rhoE_dim_cat = (sum(rho_i_cat) * e_dim_cat);
    [rho_i_supercat_uh, e_supercat_uh] = density_energy_from_YiTP(Y_i_cat, T_wall, p_flow_uh);
    rhoE_cat_uh = sum(rho_i_supercat_uh)*e_supercat_uh;
    u_cat = 0*uhat;
    u_cat(1:5) = rho_i_supercat_uh / rho_scale;
    u_cat(6:7) = 0;
    u_cat(8) = rhoE_cat_uh / rhoe_scale;
    f_cat_uh = u_cat - uhat;



        %%% Partially catalytic wall
    % J_i n = w_cat
    % Step 1: evaluate catalytic source term
    % Step 2: enforce viscous flux of continuity to 0
    % Step 3: grab wall temperature
    gam_i = eta(13:17);
    w_1 = T_wall / T_scale;
    uf = u;
    uf(8) = rhoE_wall;
    uf(6:7) = 0;

    f_species   = flux_visc_species(u, q, w, v, x, t, mu, eta);  
    f_species_iso   = flux_visc_species(uf, q, w_1, v, x, t, mu, eta);  

    fn_species = f_species(:,1)*n(1)+f_species(:,2)*n(2);
    fn_species_iso = f_species_iso(:,1)*n(1)+f_species_iso(:,2)*n(2);
    J_cat = 0*fn_species;
    for is = 1:5
        J_cat(is) = gam_i(is) * sqrt(T_wall / (2*pi)) * sqrt(RU ./ Mw(is)) .* u(is) * rho_scale; %TODO: some ambiguity here...
    end
    fn_species_iso(1:5) = fn_species_iso(1:5) + J_cat / (rho_scale*u_scale);

    f_cat_gam      = f_iso;
    % f_cat_gam(1:5) = f_cat_gam(1:5) + fn_species(1:5); %+ tau * (u(1:5) - uhat(1:5));
    f_cat_gam(1:5) = fn_species - fn_species_iso + tau *  (u(1:5) - uhat(1:5));

    % f_cat_gam      = f_noncat;
    % f_cat_gam(1:2) = f_cat_gam(1:2) - J_cat(1:2) / (rho_scale*u_scale);

    f_noncat_noflux       = f_noncat;
    f_noncat_noflux(1:5)  = f_species(:,1)*n(1)+f_species(:,2)*n(2) + tau * (u(1:5)-uhat(1:5));

    q = q(:);
    f_grad = q(1:8)*n(1) + q(9:16)*n(2) + tau*(u(:) - uhat(:));


    % supsersonic inflow, supersonic outflow, isothermal, noncat, supercat, partial cat  %inv. wall  no flux
    % fb = [f_in               f_out              f_iso    f_noncat f_cat      f_cat_gam     fh         f_grad];
    fb = [f_in               f_out              f_iso    f_noncat f_cat      f_cat_gam       fh         f_grad      f_noncat_uh   f_noncat_noflux    f_cat_uh];
end

function f = eos(u, q, w, v, x, t, mu, eta)
% Nondimensional params
    kinetics_params = kinetics();
    [species_thermo_structs, Mw, ~] = thermodynamicsModels();
    ns = kinetics_params.ns;
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    mu_scale    = eta(5);
    kappa_scale = eta(6);
    cp_scale    = eta(7);
    L_scale     = eta(8);
    Ec          = eta(9);

    rho_i = u(1:ns) * rho_scale;
    rhou = u(ns+1) * (rho_scale * u_scale);
    rhov = u(ns+2) * (rho_scale * u_scale);
    rhoE = u(ns+3) * rhoe_scale;

    rhoe = Ec * (rhoE - 0.5 * (rhou*rhou + rhov*rhov) / sum(rho_i));
    rho_tilde = rho_i ./ Mw;
    alpha = -sum(rho_tilde);
    f = f_T(w(1)*T_scale, rho_tilde, rhoe, alpha, species_thermo_structs);
end

function f = sourcew(u, q, w, v, x, t, mu, eta)
% Nondimensional params
    kinetics_params = kinetics();
    ns = kinetics_params.ns;
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    Ec          = eta(9);
    
    rho_i = u(1:ns) * rho_scale;
    rhou = u(ns+1) * (rho_scale * u_scale);
    rhov = u(ns+2) * (rho_scale * u_scale);
    rhoE = u(ns+3) * rhoe_scale;

    rhoe = Ec * (rhoE - 0.5 * (rhou*rhou + rhov*rhov) / sum(rho_i));
    f = equationofstate(w(1)*T_scale, rho_i, rhoe);
end

function m = monitor(u, q, w, v, x, t, mu, eta)
m(1) = w(1);
end
