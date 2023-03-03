function pde = pdemodel2
% pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
% pde.initw = @initw;
end

% function m = mass(u, q, w, v, x, t, mu, eta)
%     m = sym(0.0);        % Multiply by r for axisymmetric
% end

function f = flux(u, q, w, v, x, t, mu, eta)
    r = x(1);
    f = r*[q(1) q(2)];      % Poisson eqn has a (-) on the left, so it is fine to leave -grad here
end

function s = source(u, q, w, v, x, t, mu, eta)
    x1 = x(1);
    x2 = x(2);
    s = sym(1.0);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fh = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);

    fb = [fh 0];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    E_bd = mu(16);
    r_tip = mu(17);
    Ua = mu(14);
    % ub = [-Ua/(E_bd*r_tip), u(1), 0, 0, u(1)];    % Need to check that the solution variable for the current equation is w here and not u
    ub = [0, u(1)];    % Need to check that the solution variable for the current equation is w here and not u
end

function u0 = initu(x, mu, eta)
    u0 = sym(0.0);
end

% function w0 = initw(x, mu, eta)
%     w0 = sym([0.0, 0.0, 0.0]);
% end


% TODO add neumann BC to axisymmetric
% TODO set E field to 'mesh.vdg' in the pdeapp file for the convdiff model
% 'mesh.vdg=poi_sol'
% sol = npoly per elem x num components of the solution (which will be 3, potential + 2 gradients)
% for now, just assume the E field is constant throughout all time (gaussian dist). When we get that working we'll upgrade to wrapping it all in a for loop and coupling the poisson problem to convdiff
% Just use one equation for now - convect electrons in the presence of an E field

% todo urgent
% setup a github branch- notes.md can be the readme
% algorithms
% SWE concentration -> or aerodynamics
% openfoam num elements
% python import module or package