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
    m = r*[1.0];        % Multiply by r for axisymmetric
end

function f = flux(u, q, w, v, x, t, mu, eta)
    Ex = v(2);
    Ey = v(3);
    mu_e = 1.9163*((Ex^2 + Ey^2)^0.5)^(-0.25);     % Ionization coefficient [1/m]

    f = mu(5)*q + mu_e*[Ex, Ey].*u(1);      % kappa=0.01, c=(1, 0)
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = sym(0.0);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0);  % uhat(1) changed to 0
    fb = [fb 0];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [0, u(1)];
end

function u0 = initu(x, mu, eta)
    x1 = x(1);
    y1 = x(2);
    sigx = 0.01;
    sigy = 0.01;
    x0 = 0.0;
    y0 = 0.0;

    u0 = [exp(-0.5*( (x1-x0)^2/sigx^2 + (y1-y0)^2/sigy^2) )];
end
