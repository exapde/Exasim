function pde = pdemodel_full_model
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    r = x(1);
    m = r*[1.0 1.0 1.0 0];        % Multiply by r for axisymmetric
end

function f = flux(u, q, w, v, x, t, mu, eta)
    %%% Setup %%%
    % Unpacking the input vectors into variables for readability
    % u = num2cell(u);
    % q = num2cell(q);
    % x = num2cell(x);
    % mu = num2cell(mu);

    % DON'T FORGET that q=-grad(u), so for the gradients of eqns 1-3, need to negate the derivative to get the actual gradient.
    % [ne, np, nn, phi] = u{:};
    ne = u(1);
    np = u(2);
    nn = u(3);
    phi = u(4);

    dne_dr = q(1);
    dnp_dr = q(2);
    dnn_dr = q(3);
    Er = q(4);
    dne_dz = q(5);
    dnp_dz = q(6);
    dnn_dz = q(7);
    Ez = q(8);

    r = x(1);
    z = x(2);

    r0 = mu(1);
    z0 = mu(2);
    s0 = mu(3);
    Nmax = mu(4);
    e = mu(5);
    epsilon0 = mu(6);
    Ua = mu(7);
    gamma = mu(8);
    E_bd = mu(9);
    r_tip = mu(10);
    n_ref = mu(11);
    N = mu(12);

    normE = (Er^2+Ez^2)^0.5;
    % swarm = get_swarm_params(normE, N);
    % swarm = num2cell(swarm);
    % [alpha, eta2, beta, D, mue, mup, mun] = swarm{:};
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = -.0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%%
    D_star = D/(mue*E_bd*r_tip);   % Nondimensionalized diffusion coefficient

    % Diffusive/viscous flux is (+) because of point above - gradient q is passed to the function as -grad(u).
    f1 = r*(-[Er, Ez].*ne           + D_star*[dne_dr dne_dz]);
    f2 = r*(mup/mue*[Er, Ez].*np  + D_star*[dnp_dr dnp_dz]);
    f3 = r*(-mun/mue*[Er, Ez].*nn + D_star*[dnn_dr dnn_dz]);
    f4 = r*[Er Ez];

    f = [f1; f2; f3; f4];
end

function s = source(u, q, w, v, x, t, mu, eta)
    % Unpacking the input vectors into variables for readability
    % u = num2cell(u);
    % q = num2cell(q);
    % x = num2cell(x);
    % mu = num2cell(mu);

    % DON'T FORGET that q=-grad(u), so for the gradients of eqns 1-3, need to negate the derivative to get the actual gradient.
    % [ne, np, nn, phi] = u{:};
    % [dne_dr, dnp_dr, dnn_dr, Er, dne_dz, dnp_dz, dnn_dz, Ez] = q{:};
    % [r, z] = x{:};
    % [r0, z0, s0, Nmax, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N] = mu{:};

    ne = u(1);
    np = u(2);
    nn = u(3);
    phi = u(4);

    dne_dr = q(1);
    dnp_dr = q(2);
    dnn_dr = q(3);
    Er = q(4);
    dne_dz = q(5);
    dnp_dz = q(6);
    dnn_dz = q(7);
    Ez = q(8);

    r = x(1);
    z = x(2);

    r0 = mu(1);
    z0 = mu(2);
    s0 = mu(3);
    Nmax = mu(4);
    e = mu(5);
    epsilon0 = mu(6);
    Ua = mu(7);
    gamma = mu(8);
    E_bd = mu(9);
    r_tip = mu(10);
    n_ref = mu(11);
    N = mu(12);

    normE = (Er^2+Ez^2)^0.5;
    % swarm = get_swarm_params(normE, N);
    % swarm = num2cell(swarm);
    % [alpha, eta2, beta, D, mue, mup, mun] = swarm{:};
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = -.0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%%

    s1 = (alpha-eta)*r_tip*ne*normE - beta*epsilon0/(e*mue)*ne*np;
    s2 = alpha*r_tip*ne*normE - beta*epsilon0*np/(e*mue) *(nn + ne);
    s3 = eta*r_tip*ne*normE - beta*epsilon0/(e*mue)*nn*np;
    s4 = np-ne-nn;

    s = [s1, s2, s3, s4];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    % Unpacking the input vectors into variables for readability
    % u = num2cell(u);
    % q = num2cell(q);
    % x = num2cell(x);
    % mu = num2cell(mu);

    % DON'T FORGET that q=-grad(u), so for the gradients of eqns 1-3, need to negate the derivative to get the actual gradient.
    % [ne, np, nn, phi] = u{:};
    % [dne_dr, dnp_dr, dnn_dr, Er, dne_dz, dnp_dz, dnn_dz, Ez] = q{:};
    % [r, z] = x{:};
    % [r0, z0, s0, Nmax, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N] = mu{:};
    ne = u(1);
    np = u(2);
    nn = u(3);
    phi = u(4);

    dne_dr = q(1);
    dnp_dr = q(2);
    dnn_dr = q(3);
    Er = q(4);
    dne_dz = q(5);
    dnp_dz = q(6);
    dnn_dz = q(7);
    Ez = q(8);

    r = x(1);
    z = x(2);

    r0 = mu(1);
    z0 = mu(2);
    s0 = mu(3);
    Nmax = mu(4);
    e = mu(5);
    epsilon0 = mu(6);
    Ua = mu(7);
    gamma = mu(8);
    E_bd = mu(9);
    r_tip = mu(10);
    n_ref = mu(11);
    N = mu(12);

    normE = (Er^2+Ez^2)^0.5;
    % swarm = get_swarm_params(normE, N);
    % swarm = num2cell(swarm);
    % [alpha, eta2, beta, D, mue, mup, mun] = swarm{:};
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = -.0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%%

    f = flux(u, q, w, v, x, t, mu, eta);
    ndotE = n(1)*Er + n(2)*Ez;
    alpha = 1000*tanh(ndotE);
    % alpha is 0 when ndotE < 0 => neumann
    % alpha is 1 when ndotE > 0 => homogeneous dirichlet
    % Undefined/0.5 when ndotE = 0
    % Note that the BCs are nondimensionalized
    D_star = D/(mue*E_bd*r_tip);   % Nondimensionalized diffusion coefficient

    % Total flux
    fb_e = f(1,1)*n(1) + f(1,2)*n(2) + tau*(u(1)-uhat(1));
    fb_p = f(2,1)*n(1) + f(2,2)*n(2) + tau*(u(2)-uhat(2));
    fb_n = f(3,1)*n(1) + f(3,2)*n(2) + tau*(u(3)-uhat(3));
    fb_poi = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4));

    % Diffusive flux component removed - "convective flux only"
    fb_conv_e = -ne*             (n(1)*Er + n(2)*Ez) + tau*(u(1)-uhat(1));    % Bringing (-) out front
    fb_conv_p =  np*(mup/mue)* (n(1)*Er + n(2)*Ez) + tau*(u(2)-uhat(2));
    fb_conv_n = -nn*(mun/mue)* (n(1)*Er + n(2)*Ez) + tau*(u(3)-uhat(3));

    % Read "flux for equation x, boundary y". Note it is the flux, not the derivative of the function. The two are different with the addition of the convective term.
    feq1b1 = (gamma*np*normE - ne*(n(1)*Er + n(2)*Ez))/D_star;
    feq2b1 = fb_conv_p;
    feq3b1 = fb_n;
    feq4b1 = fb_poi;
    
    % BC 2 is an axisymmetric BC, so inhomogeneous neumann conditions are applied to each equation.
    feq1b2 = 0;
    feq2b2 = 0;
    feq3b2 = 0;
    feq4b2 = 0;
    
    feq1b3 = (1-alpha)*fb_conv_e + alpha*fb_e;
    feq2b3 = (1-alpha)*fb_conv_p + alpha*fb_p;
    feq3b3 = (1-alpha)*fb_conv_n + alpha*fb_n;
    feq4b3 = fb_poi;
    
    feq1b4 = fb_conv_e;
    feq2b4 = fb_p;
    feq3b4 = fb_conv_n;
    feq4b4 = fb_poi;
    
    feq1b5 = (1-alpha)*fb_conv_e + alpha*fb_e;
    feq2b5 = (1-alpha)*fb_conv_p + alpha*fb_p;
    feq3b5 = (1-alpha)*fb_conv_n + alpha*fb_n;
    feq4b5 = 0;

    feq1 = [feq1b1; feq1b2; feq1b3; feq1b4; feq1b5].';
    feq2 = [feq2b1; feq2b2; feq2b3; feq2b4; feq2b5].';
    feq3 = [feq3b1; feq3b2; feq3b3; feq3b4; feq3b5].';
    feq4 = [feq4b1; feq4b2; feq4b3; feq4b4; feq4b5].';
    
    fb = [feq1; feq2; feq3; feq4];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    % Unpacking the input vectors into variables for readability
    % u = num2cell(u);
    % q = num2cell(q);
    % x = num2cell(x);
    % mu = num2cell(mu);

    % DON'T FORGET that q=-grad(u), so for the gradients of eqns 1-3, need to negate the derivative to get the actual gradient.
    % [ne, np, nn, phi] = u{:};
    % [dne_dr, dnp_dr, dnn_dr, Er, dne_dz, dnp_dz, dnn_dz, Ez] = q{:};
    % [r, z] = x{:};
    % [r0, z0, s0, Nmax, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N] = mu{:};

    ne = u(1);
    np = u(2);
    nn = u(3);
    phi = u(4);

    dne_dr = q(1);
    dnp_dr = q(2);
    dnn_dr = q(3);
    Er = q(4);
    dne_dz = q(5);
    dnp_dz = q(6);
    dnn_dz = q(7);
    Ez = q(8);

    r = x(1);
    z = x(2);

    r0 = mu(1);
    z0 = mu(2);
    s0 = mu(3);
    Nmax = mu(4);
    e = mu(5);
    epsilon0 = mu(6);
    Ua = mu(7);
    gamma = mu(8);
    E_bd = mu(9);
    r_tip = mu(10);
    n_ref = mu(11);
    N = mu(12);

    normE = (Er^2+Ez^2)^0.5;
    % swarm = get_swarm_params(normE, N);
    % swarm = num2cell(swarm);
    % [alpha, eta2, beta, D, mue, mup, mun] = swarm{:};
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = -.0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%%

    f = flux(u, q, w, v, x, t, mu, eta);
    ndotE = n(1)*Er + n(2)*Ez;
    alpha = 1000*tanh(ndotE);
    % alpha is 0 when ndotE < 0
    % alpha is 1 when ndotE > 0
    % Undefined/0.5 when ndotE = 0
    % Note that the BCs are nondimensionalized

    ueq1b1 = ne;
    ueq2b1 = np;
    ueq3b1 = 0;
    ueq4b1 = -Ua/(E_bd*r_tip);
    
    % BC 2 is an axisymmetric BC, so inhomogeneous neumann conditions are applied to each equation.
    ueq1b2 = ne;
    ueq2b2 = np;
    ueq3b2 = nn;
    ueq4b2 = phi;
    
    ueq1b3 = (1-alpha)*ne + alpha*0;
    ueq2b3 = (1-alpha)*np + alpha*0;
    ueq3b3 = (1-alpha)*nn + alpha*0;
    ueq4b3 = 0;
    
    ueq1b4 = ne;
    ueq2b4 = 0;
    ueq3b4 = nn;
    ueq4b4 = 0;
    
    ueq1b5 = (1-alpha)*ne + alpha*0;
    ueq2b5 = (1-alpha)*np + alpha*0;
    ueq3b5 = (1-alpha)*nn + alpha*0;
    ueq4b5 = phi;

    ueq1 = [ueq1b1; ueq1b2; ueq1b3; ueq1b4; ueq1b5].';
    ueq2 = [ueq2b1; ueq2b2; ueq2b3; ueq2b4; ueq2b5].';
    ueq3 = [ueq3b1; ueq3b2; ueq3b3; ueq3b4; ueq3b5].';
    ueq4 = [ueq4b1; ueq4b2; ueq4b3; ueq4b4; ueq4b5].';

    ub = [ueq1; ueq2; ueq3; ueq4];
end

function u0 = initu(x, mu, eta)
    mu = num2cell(mu);
    [r0, z0, s0, Nmax, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N] = mu{:};

    % Nondimensionalization
    r0_star = r0/r_tip;
    z0_star = z0/r_tip;
    sig0_star = s0/r_tip;
    n_max_star = Nmax/n_ref;

    r = x(1);
    z = x(2);

    % "star" indicating it is a nondimensional quantity
    n0_star = n_max_star* exp(-1/(2*sig0_star^2) * ((r-r0_star)^2 + (z-z0_star)^2));

    u0 = [n0_star; n0_star; 0; 0];    % Only set the seed particles for the electrons and positives
end

function params = get_swarm_params(normE, N)
    % Using the values in the appendix of this paper: https://iopscience.iop.org/article/10.1088/0022-3727/30/4/017/pdf
    % Note that the reduced E field is converted to V*cm2 for performing the calculations, then the resultant quantity is converted back from cm->m.
    % Expects the reduced E field in V*m2. Should be computed by dividing the L2 norm of the E field by the number density of neutrals.

    cm2m = 0.01;
    m2cm = 100;
    E_N = normE/N*m2cm^2;   % Reduced electric field in V*cm2. Expects E to be the L2 norm of the E field.

    % First ionization coefficient
    if E_N <= 1.5e-15
        alpha = 6.619e-17*exp(-5.593e-15/E_N)*cm2m^2;
    elseif E_N > 1.5e-15
        alpha = 2e-16*exp(-7.248e-15/E_N)*cm2m^2;
    end

    % Two-body attachment coefficient
    if E_N <= 1.05e-15
        eta2 = 6.089e-4*E_N - 2.893e-19*cm2m^2;
    elseif E_N > 1.05e-15
        eta2 = 8.889e-5*E_N + 2.567e-19*cm2m^2;
    end

    % Recombination coefficient
    beta = 2e-7*cm2m^3;

    % Electron mobility
    % Removing the sign term because these formulas were designed for computing the velocity, not the mobility. The mobility always has the same sign.
    if E_N <= 2.6e-17
        mue = -(3.38e4 + 6.87e22*E_N) * cm2m / E;   % Formula is for the velocity, so need to divide by the E field strength to get the mobility
    elseif ((E_N > 2.6e-17) && (E_N <= 1e-16))
        mue = -(1.63e6 + 7.2973e21*E_N) * cm2m / E;
    elseif ((E_N > 1e-16) && (E_N <= 2e-15))
        mue = -(1.3e6 + 1.03e22*E_N) * cm2m / E;
    elseif E_N > 2e-15
        mue = -(7.1e6 + 7.4e21*E_N) * cm2m / E;
    end
    
    % Negative ion mobility
    if E_N <= 5e-16
        mun = -1.86*cm2m^2;
    elseif E_N > 5e-16
        mun = -2.7*cm2m^2;
    end

    % Positive ion mobility
    % Assuming that P0/P is approximately 1 and that the units for the E field are V/cm
    mup = 2.34*cm2m^2;

    % Electron diffusion coefficient
    % Unit conversion (cm->m) is "baked in" through the mue calculated previously
    D = mue*0.3341e9*E_N^0.54069;

    params = [alpha, eta2, beta, D, mue, mup, mun];
end
% Notes: changed indexing into q, changed the way we handle the phi solution so that mue isn't singular

% Next step: make it nondimensional, check to get the poisson source term working (observe the potential changing at each time step)

% Need to show phi changing
% Need to plot the superposition of phi solutions together

% Is the time solution snapshot saved at the beginning of the timestep or at the end? i.e. is the IC saved or not?
% Sometimes it's difficult to tell if the sim has advanced at the first snapshot, or if it's the IC.

% Separate out the two components of phi when plotting
% Run the electrostatics case separately with the ne forcing term (even a small background charge should have an effect)
% Can I compute quantities like the reduced E field and swarm parameters without haveing to recompute the same values in each function?
% Change model to nondimensionalize by the mue computed at the Ebd strength
% and at STP (-> N). This will require a tweak to the mue in the model,
% probably by a ratio of mue_ref/mue or something like that.