function pde = pdemodel_poisson
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    end
    
    function f = flux(u, q, w, v, x, t, mu, eta)
        r = x(1);
        f = r*[q(1) q(2)];
    end
    
    function s = source(u, q, w, v, x, t, mu, eta)
    s = sym(0.0);
    end
    
    function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fh = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));        % Change to "uhat(1))" from "0.0" fixed the problem
    fb = [fh 0 fh];
    end
    
    function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau) 
    ub = [0 u(1) 1];
    end
    
    function u0 = initu(x, mu, eta)
    u0 = sym(0.0);
    end
    