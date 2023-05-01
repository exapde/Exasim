function pde = pdemodel_dd_orig
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    r = x(1);
    m = r*[1.0 0];        % Multiply by r for axisymmetric
end

function f = flux(u, q, w, v, x, t, mu, eta)
    %%%%%% SETUP %%%%%%
    % DON'T FORGET TO CHANGE THE INDEXING BACK TO THE FULL MODEL.
    ne = u(1);
    phi = u(2);

    dne_dr = q(1);
    Er = q(2) + v(2);
    dne_dz = q(3);
    Ez = q(4) + v(3);

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
    mue_ref = mu(13);

    normE = (Er^2+Ez^2)^0.5;
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = .0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%% END SETUP %%%%%%

    D_star = D/(mue_ref*E_bd*r_tip);   % Nondimensionalized diffusion coefficient

    % Diffusive/viscous flux is (+) because gradient q is passed to the function as -grad(u).
    f1 = r*(-(mue/mue_ref)*[Er, Ez].*ne + D_star*[dne_dr dne_dz]);
    % f2 = r*(mup/mue_ref*[Er, Ez].*np  + D_star*[dnp_dr dnp_dz]);
    % f3 = r*(-(mun/mue_ref)*[Er, Ez].*nn + D_star*[dnn_dr dnn_dz]);
    f4 = r*[Er Ez];
    f = [f1; f4];
end

function s = source(u, q, w, v, x, t, mu, eta)
    %%%%%% SETUP %%%%%%
    % DON'T FORGET TO CHANGE THE INDEXING BACK TO THE FULL MODEL.
    ne = u(1);
    phi = u(2);

    dne_dr = q(1);
    Er = q(2) + v(2);
    dne_dz = q(3);
    Ez = q(4) + v(3);

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
    mue_ref = mu(13);

    normE = (Er^2+Ez^2)^0.5;
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = .0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%% END SETUP %%%%%%

    % s1 = (alpha-eta)*(mue/mue_ref)*r_tip*ne*normE - beta*epsilon0/(e*mue_ref)*ne*np;
    % s2 = alpha*(mue/mue_ref)*r_tip*ne*normE - beta*epsilon0*np/(e*mue_ref) *(nn + ne);
    % s3 = eta*(mue/mue_ref)*r_tip*ne*normE - beta*epsilon0/(e*mue_ref)*nn*np;
    % s4 = ne+nn-np;  % Added the minus sign

    % s = [s1, s2, s3, s4];
    s = [sym(0.0), sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    %%%%%% SETUP %%%%%%
    % DON'T FORGET TO CHANGE THE INDEXING BACK TO THE FULL MODEL.
    ne = u(1);
    phi = u(2);

    dne_dr = q(1);
    Er = q(2) + v(2);
    dne_dz = q(3);
    Ez = q(4) + v(3);

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
    mue_ref = mu(13);

    normE = (Er^2+Ez^2)^0.5;
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = .0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%% END SETUP %%%%%%
    
    f = flux(u, q, w, v, x, t, mu, eta);
    fb_e = f(1,1)*n(1) + f(1,2)*n(2) + tau*(u(1)-uhat(1));
    % fb_p = f(2,1)*n(1) + f(2,2)*n(2) + tau*(u(2)-uhat(2));
    % fb_n = f(3,1)*n(1) + f(3,2)*n(2) + tau*(u(3)-uhat(3));
    fb_poi = f(2,1)*n(1) + f(2,2)*n(2) + tau*(u(2)-uhat(2)); % CHANGE INDEXING

    ndotE = n(1)*Er + n(2)*Ez;
    alpha = 0.5*(tanh(1000*ndotE)+1);

    fb_conv_e = -r*ne*(mue/mue_ref)*(n(1)*Er + n(2)*Ez) + tau*(u(1)-uhat(1));    % Bringing (-) out front
    % fb_conv_p =  r*np*(mup/mue_ref)*(n(1)*Er + n(2)*Ez) + tau*(u(2)-uhat(2));    % Bringing (-) out front
    % fb_conv_n = -r*nn*(mun/mue_ref)*(n(1)*Er + n(2)*Ez) + tau*(u(3)-uhat(3));    % Bringing (-) out front

    % 1. Needle
    % 2. Axisymmetric
    % 3. Outflow boundary at bottom
    % 4. Grounded cylinder and bottom plane
    % 5. Farfield (top and right)

    feq1b1 = 0;   % CHANGE ME  % Note it is the flux, not the derivative of the function. The two are different with the addition of the convective term.
    feq1b2 = 0;
    feq1b3 = (1-alpha)*fb_conv_e + alpha*fb_e;        % BC option 1
    feq1b4 = fb_conv_e;
    feq1b5 = 0;

    feq4b1 = fb_poi;
    feq4b2 = 0;
    feq4b3 = fb_poi;
    feq4b4 = fb_poi;
    feq4b5 = 0;

    % NOTE: WILL NEED TO NEGATE THE ALPHA FCN FOR THE POSITIVES

    feq1 = [feq1b1; feq1b2; feq1b3; feq1b4; feq1b5].';
    feq4 = [feq4b1; feq4b2; feq4b3; feq4b4; feq4b5].';
    
    fb = [feq1; feq4];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    %%%%%% SETUP %%%%%%
    % DON'T FORGET TO CHANGE THE INDEXING BACK TO THE FULL MODEL.
    ne = u(1);
    phi = u(2);

    dne_dr = q(1);
    Er = q(2) + v(2);
    dne_dz = q(3);
    Ez = q(4) + v(3);

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
    mue_ref = mu(13);

    normE = (Er^2+Ez^2)^0.5;
    alpha = 1.19e-21;
    eta = 2.28e-19;
    D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
    beta = 2e-13;
    mue = .0378;
    mup = 2.34e-4;
    mun = -2.7e-4;
    %%%%%% END SETUP %%%%%%
    
    ndotE = n(1)*Er + n(2)*Ez;
    alpha = 0.5*(tanh(1000*ndotE)+1);
    
    ueq1b1 = ne;
    ueq1b2 = ne;
    ueq1b3 = (1-alpha)*ne + alpha*0;
    ueq1b4 = ne;
    ueq1b5 = ne;

    ueq4b1 = 0;
    ueq4b2 = phi;
    ueq4b3 = 0;
    ueq4b4 = 0;
    ueq4b5 = phi;

    ueq1 = [ueq1b1; ueq1b2; ueq1b3; ueq1b4; ueq1b5].';
    ueq4 = [ueq4b1; ueq4b2; ueq4b3; ueq4b4; ueq4b5].';

    ub = [ueq1; ueq4];
end

function u0 = initu(x, mu, eta)
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
    mue_ref = mu(13);

    % Nondimensionalization
    r0_star = r0/r_tip;
    z0_star = z0/r_tip;
    sig0_star = s0/r_tip;
    n_max_star = Nmax/n_ref;

    % "star" indicating it is a nondimensional quantity
    n0_star = n_max_star* exp(-1/(2*sig0_star^2) * ((r-r0_star)^2 + (z-z0_star)^2));

    u0 = [n0_star; 0];    % Only set the seed particles for the electrons and positives
end


% Notes: changed indexing into q, changed the way we handle the phi solution so that mu_e isn't singular

% Next step: make it nondimensional, check to get the poisson source term working (observe the potential changing at each time step)

% Need to show phi changing
% Need to plot the superposition of phi solutions together

% Separate out the two components of phi when plotting
% Run the electrostatics case separately with the ne forcing term (even a small background charge should have an effect)

% Notes:
% With a negative sign for the convective velocity, the boundary seems to be anchored at 0 and the cloud moves away. Otherwise the cloud seems to move toward the tip with a positive sign for the convective velocity.


% Ask cuong why the diffusive flux should be 0 - probably because electrons don't diffuse across the surface, they have to get transported into it by the drift velocity


% Debugging steps
% sol contains the -grad(u)
% Plot everything: does the switch function  behave the way you want it to?
% Does the electric field look the way you want it to?
% Is at least one portion of the domain held at a dirichlet BC?
% Did you plot the alpha as a function of tanh(x) beforehand to verify that it was correct? Did you check it either in Desmos (or preferably both) a volume plot of the quantity in Erasim?
