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
    m = r*[1.0 1.0 1.0 0];        % Multiply by r for axisymmetric
end

function f = flux(u, q, w, v, x, t, mu, eta)
    r = x(1);
    Ex = v(2) + q(4);
    Ey = v(3) + q(8);
    % mu_e = 1.9163*((Ex^2 + Ey^2)^0.5)^(-0.25);     % Ionization coefficient [1/m]
    mu_e = 1;

    f1 = r*(mu(5)*[q(1) q(5)] + mu_e*[Ex, Ey].*u(1));
    f2 = r*(mu(5)*[q(2) q(6)] + mu_e*[Ex, Ey].*u(2));
    f3 = r*(mu(5)*[q(3) q(7)] + mu_e*[Ex, Ey].*u(3));
    f4 = r*[Ex Ey];
    f = [f1; f2; f3; f4];
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0), sym(0.0), sym(0.0), mu(12)/mu(13)*u(1)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);

    fb_e = f(1,1)*n(1) + f(1,2)*n(2) + tau*(u(1)-uhat(1));
    fb_p = f(2,1)*n(1) + f(2,2)*n(2) + tau*(u(2)-uhat(2));
    fb_n = f(3,1)*n(1) + f(3,2)*n(2) + tau*(u(3)-uhat(3));
    fb_poi = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4));

    % fb1 = [0; 0; 0; fb_e; 0].';
    % fb2 = [0; 0; 0; fb_p; 0].';
    % fb3 = [0; 0; 0; fb_n; 0].';
    % fb4 = [fb_poi; 0; 0; fb_poi; 0].';
    % fb = [fb1; fb2; fb3; fb4];

    % Read "flux for equation x, boundary y"
    feq1b1 = 0;
    feq2b1 = 0;
    feq3b1 = 0;
    feq4b1 = fb_poi;
    
    % BC 2 is an axisymmetric BC, so inhomogeneous neumann conditions are applied to each equation.
    feq1b2 = 0;
    feq2b2 = 0;
    feq3b2 = 0;
    feq4b2 = 0;
    
    feq1b3 = 0;
    feq2b3 = 0;
    feq3b3 = 0;
    feq4b3 = 0;
    
    feq1b4 = fb_e;
    feq2b4 = fb_p;
    feq3b4 = fb_n;
    feq4b4 = fb_poi;
    
    feq1b5 = 0;
    feq2b5 = 0;
    feq3b5 = 0;
    feq4b5 = 0;

    feq1 = [feq1b1; feq1b2; feq1b3; feq1b4; feq1b5].';
    feq2 = [feq2b1; feq2b2; feq2b3; feq2b4; feq2b5].';
    feq3 = [feq3b1; feq3b2; feq3b3; feq3b4; feq3b5].';
    feq4 = [feq4b1; feq4b2; feq4b3; feq4b4; feq4b5].';
    
    fb = [feq1; feq2; feq3; feq4];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    % ub1 = [u(1); u(1); u(1); 0; u(1)].';
    % ub2 = [u(2); u(2); u(2); 0; u(2)].';
    % ub3 = [u(3); u(3); u(3); 0; u(3)].';
    % ub4 = [0;    u(4); u(4); 0; u(4)].';
    % ub = [ub1; ub2; ub3; ub4];


    ne = u(1);
    np = u(2);
    nn = u(3);
    phi = u(4);

    ueq1b1 = ne;
    ueq2b1 = np;
    ueq3b1 = nn;
    ueq4b1 = 0;
    
    % BC 2 is an axisymmetric BC, so inhomogeneous neumann conditions are applied to each equation.
    ueq1b2 = ne;
    ueq2b2 = np;
    ueq3b2 = nn;
    ueq4b2 = phi;
    
    ueq1b3 = ne;
    ueq2b3 = np;
    ueq3b3 = nn;
    ueq4b3 = phi;
    
    ueq1b4 = 0;
    ueq2b4 = 0;
    ueq3b4 = 0;
    ueq4b4 = 0;
    
    ueq1b5 = ne;
    ueq2b5 = np;
    ueq3b5 = nn;
    ueq4b5 = phi;

    ueq1 = [ueq1b1; ueq1b2; ueq1b3; ueq1b4; ueq1b5].';
    ueq2 = [ueq2b1; ueq2b2; ueq2b3; ueq2b4; ueq2b5].';
    ueq3 = [ueq3b1; ueq3b2; ueq3b3; ueq3b4; ueq3b5].';
    ueq4 = [ueq4b1; ueq4b2; ueq4b3; ueq4b4; ueq4b5].';

    ub = [ueq1; ueq2; ueq3; ueq4];
end

function u0 = initu(x, mu, eta)
    x1 = x(1);
    y1 = x(2);
    sige = 0.01;
    sigp = 0.013;
    sign = 0.007;
    x0 = 0.0;
    y0 = 0.0;

    u0_ne = exp(-0.5*( (x1-x0)^2/sige^2 + (y1-y0)^2/sige^2)); 
    u0_np = exp(-0.5*( (x1-x0)^2/sigp^2 + (y1-y0)^2/sigp^2)); 
    u0_nn = exp(-0.5*( (x1-x0)^2/sign^2 + (y1-y0)^2/sign^2)); 

    u0 = [u0_ne; u0_np; u0_nn; 0];
end

% Notes: changed indexing into q, changed the way we handle the phi solution so that mu_e isn't singular

% Next step: make it nondimensional, check to get the poisson source term working (observe the potential changing at each time step)

% Need to show phi changing
% Need to plot the superposition of phi solutions together

% Separate out the two components of phi when plotting
% Run the electrostatics case separately with the ne forcing term (even a small background charge should have an effect)