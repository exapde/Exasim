function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
% pde.avfield = @avfield;
pde.fbouhdg = @fbouhdg;
pde.fhathdg = @fhathdg;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)

f = fluxaxialns(u, q, w, v, x, t, mu, eta);

end

% function f = avfield(u, q, w, v, x, t, mu, eta)
%     f = getavfield2d(u,q,v,mu);
% end

function s = source(u, q, w, v, x, t, mu, eta)

s = sourceaxialns(u, q, w, v, x, t, mu, eta);

end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
   
    fb = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

fb = fbouhdgaxialns(u, q, w, v, x, t, mu, eta, uhat, n, tau);

end

function fh = fhathdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

f = flux(uhat, q, w, v, x, t, mu, eta);

r = uhat(1);
ru = uhat(2);
rv = uhat(3);
nx = n(1);
ny = n(2);    
run = ru*nx + rv*ny;
r1 = 1/r;
un = run*r1;

s = (u - uhat);
for i = 1:4
    s(i) = r*tau(i)*(u(i) - uhat(i));
end
% s(2) = sqrt(un*un + 1e-10)*(u(2) - uhat(2)); 
% s(3) = sqrt(un*un + 1e-10)*(u(3) - uhat(3)); 

fh = f(:,1)*n(1) + f(:,2)*n(2) + s;

end

function u0 = initu(x, mu, eta)
    u0 = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end


