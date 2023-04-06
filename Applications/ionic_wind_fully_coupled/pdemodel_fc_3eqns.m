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
    fb_conv_e = f(1,1)*n(1) + f(1,2)*n(2) + tau*(u(1)-uhat(1));
    fb1 = [0; 0; 0; fb_conv_e; 0].';

    fb_conv_p = f(2,1)*n(1) + f(2,2)*n(2) + tau*(u(2)-uhat(2));
    fb2 = [0; 0; 0; fb_conv_p; 0].';

    fb_conv_n = f(3,1)*n(1) + f(3,2)*n(2) + tau*(u(3)-uhat(3));
    fb3 = [0; 0; 0; fb_conv_n; 0].';

    fb_poi = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4));
    fb4 = [fb_poi; 0; 0; fb_poi; 0].';

    fb = [fb1; fb2; fb3; fb4];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub1 = [u(1); u(1); u(1); 0; u(1)].';
    ub2 = [u(2); u(2); u(2); 0; u(2)].';
    ub3 = [u(3); u(3); u(3); 0; u(3)].';
    ub4 = [0; u(4); u(4); 0; u(4)].';

    ub = [ub1; ub2; ub3; ub4];
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