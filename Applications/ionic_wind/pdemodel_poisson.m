function pde = pdemodel
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
        % s = v(1);
    s = sym(0.0);
    end
    
    function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fh = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
    fb = [fh 0 fh]
    end
    
    function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    % ub = sym(0.0); 
    ub = [0 u(1) 1]
    end
    
    function u0 = initu(x, mu, eta)
    u0 = sym(0.0);
    end
    