function f = flux(u, q, w, v, x, t, mu, eta)
    nch = 8;
    nd = 2;
    ns = 5;

    nenergy = 1;
    % ndim = 2;
    [species_thermo_structs, Mw, ~] = thermodynamicsModels();

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

    % Conservative Variables
    for ispecies = 1:ns
        rho_i(ispecies) = u(ispecies);
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(ns+1);
    rhov = u(ns+2);
    rhoE = u(ns+3);

    drho_dx_i = -q(1:ns);
    drho_dy_i = -q((nch+1:nch+ns));

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

    % Viscous fluxes
    Ec          = eta(9);
    Pr          = eta(10);
    Re          = eta(11);

    [blottner_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport();
    drho_dx = sum(drho_dx_i);
    drho_dy = sum(drho_dy_i);
    X = X_i(rho_i_dim,Mw);
    
    Y = Y_i(rho_i_dim);

    D_vec = averageDiffusionCoeffs(T_dim, X, Y, Mw, p_dim, gupta_structs);

    D_scale = u_scale;

    D_vec = D_vec / D_scale;

    %%%%%%%% Calculation of J_i
    dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rhoinv * rhoinv;
    dY_dy_i = (drho_dy_i * rho - rho_i * drho_dy) * rhoinv * rhoinv;

    J_i_x = -rho * D_vec .* dY_dx_i + 0*rho_i .* sum(D_vec .* dY_dx_i);
    J_i_y = -rho * D_vec .* dY_dy_i + 0*rho_i .* sum(D_vec .* dY_dy_i);
    
    for i = 1:ns
        fv(i,1) = -J_i_x(i); 
        fv(i,2) = -J_i_y(i);
    end

    f = -fv;
end
