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

function m = mass(u, q, w, v, x, t, mu, eta)
    ns = 5;
    ndim = numel(x);
    m = sym(ones(ns + ndim + 1, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)
    f = fluxaxial2d(u, q, w, v, x, t, mu, eta);
end

function s = source(u, q, w, v, x, t, mu, eta)        
    s = sourceaxial2d(u, q, w, v, x, t, mu, eta);    
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    ndim = numel(x);
    ub = sym(ones(ns + ndim + 1, 1));
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    ndim = numel(x);
    fb = sym(ones(ns + ndim + 1, 1));
end

function u0 = initu(x, mu, eta)
    ns = 5;
    ndim = numel(x);
    u0 = sym(ones(ns + ndim + 1, 1));    
end

function w0 = initw(x, mu, eta)
    w0 = sym(ones(1,1));
end

function f = eos(u, q, w, v, x, t, mu, eta)
    f = eosnd(u, q, w, v, x, t, mu, eta);    
end

function f = sourcew(u, q, w, v, x, t, mu, eta)
    f = eosnd(u, q, w, v, x, t, mu, eta);        
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
  fb = fbouhdgaxialnd(u, q, w, v, x, t, mu, eta, uhat, n, tau);
end

