 
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
end

% function sw = sourcew(u, q, w, v, x, t, mu, eta)
%     ns = 5;
%     %should not be used, should be overwritten mutationOutputs_2d_dim.cpp
%     sw = sym(zeros(26,1));
% 
%     rmin = 0.0;
%     alphaClip = mu(19);
% 
%     rho_i = rmin + lmax(u(1:ns) - rmin, alphaClip);
%     sw(1:ns) = rho_i;
%     sw(6) = 0;
%     sw(7) = 0;
% end

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
    nd = 2;
    ns = 5;

    nenergy = 1;
    % ndim = 2;
    [~, Mw, ~] = thermodynamicsModels();

    % Nondimensional params
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    mu_scale    = eta(5);
    kappa_scale = eta(6);
    cp_scale    = eta(7);
    L_scale     = eta(8);

    rho_i = sym(zeros(ns,1));
    rho = sym(0);

    rmin = 0.0;
    alphaClip = mu(19);

    % Conservative Variables
    for ispecies = 1:ns
%         rho_i(ispecies) = rmin + lmax(u(ispecies)-rmin,alphaClip); %subspecies density
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
    T = w(1) *T_scale;
    p = pressure(T, rho_i_dim, Mw) / rhoe_scale;

    H = E + p*rhoinv; %enthalpy

    % Fluxes
    for ispecies = 1:ns
        f(ispecies,1) = rho_i(ispecies) * uv - av.*drho_dx_i(ispecies);
    end
    f(ns + 1,1) = rhou * uv + p     - av.*drhou_dx;
    f(ns + 2,1) = rhov * uv         - av.*drhov_dx;
    f(ns + 3,1) = rhou * H          - av.*drhoE_dx;

    for ispecies = 1:ns
        f(ispecies,2) = rho_i(ispecies) * vv - av.*drho_dy_i(ispecies);
    end
    f(ns + 1,2) = rhou * vv      - av.*drhou_dy;
    f(ns + 2,2) = rhov * vv + p  - av.*drhov_dy;
    f(ns + 3,2) = rhov * H       - av.*drhoE_dy;

    % fi = [ru, ru*uv+p, rv*uv, ru*h, ...
            % rv, ru*vv, rv*vv+p, rv*h];
    % fi = reshape(fi,[4,2]);    
end



function s = source(u, q, w, v, x, t, mu, eta)
    
    ns = 5;
    rmin = 0.0;

    [species_thermo_structs, Mw, ~] = thermodynamicsModels();
    kinetics_params = kinetics();

    % Nondimensional params
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    mu_scale    = eta(5);
    kappa_scale = eta(6);
    cp_scale    = eta(7);
    L_scale     = eta(8);
    omega_scale = rho_scale * u_scale / L_scale;

    % Mutation outputs
    rho_i_dim = u(1:ns) * rho_scale;
    T = w(1) * T_scale;
    omega_i = netProductionRatesTotal(rho_i_dim, T, Mw, kinetics_params, species_thermo_structs);

    s(1:ns) = omega_i / omega_scale;
    s(ns+1) = 0.0;
    s(ns+2) = 0.0;
    s(ns+3) = 0.0;
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    ndim = 2;
    [species_thermo_structs, Mw, ~] = thermodynamicsModels();

    ub = sym(zeros(ns+ndim+1, 3));
    rho_scale   = eta(1);
    u_scale     = eta(2);
    rhoe_scale  = eta(3);
    T_scale     = eta(4);
    mu_scale    = eta(5);
    kappa_scale = eta(6);
    cp_scale    = eta(7);
    L_scale     = eta(8);

    Ec = mu(23);
%%%%%% C-CODE MANUALLY WRITTEN
%     ub(:,1) = u(:);
%     ub(:,2) = u(:); 

    uinflow = initu(x, mu, eta);

    uoutflow = u;

    uadiabatic = u;
    uadiabatic(ns+1:ns+ndim) = 0.0;


    ub(:,1) = uinflow;
    ub(:,2) = uoutflow;
    ub(:,3) = uadiabatic;
    ub(:,4) = uinflow;
    ub(:,5) = uoutflow;
    ub(:,6) = uadiabatic;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    ns = 5;
    ndim = 2;
    fb = sym(zeros(ns+ndim+1, 6));

%     fiso(1:ns) = 0.0;
    % fb(:,1) = 0;
    % fb(:,2) = 0;
    % fb(:,3) = 0;
    % fb(:,1) = 0;
    % fb(:,2) = 0;
    % fb(:,3) = 0;
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
    fb = [f_in f_in fh fh f_in f_out];
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

    Ec = 1;

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

    Ec = mu(23);

    rho_i = u(1:ns) * rho_scale;
    rhou = u(ns+1) * (rho_scale * u_scale);
    rhov = u(ns+2) * (rho_scale * u_scale);
    rhoE = u(ns+3) * rhoe_scale;

    rhoe = Ec * (rhoE - 0.5 * (rhou*rhou + rhov*rhov) / sum(rho_i));
    rho_tilde = rho_i ./ Mw;
    alpha = -sum(rho_tilde);
    f = f_T(w(1)*T_scale, rho_tilde, rhoe, alpha, species_thermo_structs);
end