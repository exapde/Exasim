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
    m = r*[1.0, 0];        % Multiply by r for axisymmetric
end

function f = flux(u, q, w, v, x, t, mu, eta)
    r = x(1);
    Ex = q(3);
    Ey = q(4);
    mu_e = 1.9163*((Ex^2 + Ey^2)^0.5)^(-0.25);     % Ionization coefficient [1/m]

    f1 = r*(mu(5)*[q(1) q(2)] + mu_e*[Ex, Ey].*u(1));      % kappa=0.01, c=(1, 0)
    f2 = r*[Ex Ey];
    f = [f1; f2];
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb_pre = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
    fb = [0 fb_pre(2); 0 0; 0 0; fb_pre(1) fb_pre(2); 0 0].';
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [u(1) 1; u(1) u(2); u(1) u(2); 0 0; u(1) u(2)].';
end

function u0 = initu(x, mu, eta)
    x1 = x(1);
    y1 = x(2);
    sigx = 0.01;
    sigy = 0.01;
    x0 = 0.0;
    y0 = 0.0;

    u0 = [exp(-0.5*( (x1-x0)^2/sigx^2 + (y1-y0)^2/sigy^2)); 1];
end

