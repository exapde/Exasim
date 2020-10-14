function pde = pdemodel
pde.extsource1 = @extsource1;
pde.extsource2 = @extsource2;
end

%function f = flux(u, q, w, v, x, t, mu, eta)
% sources

function fi = finterface1(u, q, w, v, x, t, mu, eta, uhat, vhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fi = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end
function ub = uinterface1(u, q, w, v, x, t, mu, eta, uhat, vhat, n, tau)
ub = sym(0.0); 
end
function s = externalsource1(U, Q, W, V, x, t, mu, eta)
s = U{2};
end

function s = extsource2(U, Q, W, V, x, t, mu, eta)
s = U{1};
end

