function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    r = x(1);
    m = r*[1.0, 1.0, 1.0];        % Multiply by r for axisymmetric
end

function f = flux(u, q, w, v, x, t, mu, eta)
    Ex = v(2);
    Ey = v(3);
    mu_e = 1.9163*((Ex^2 + Ey^2)^0.5)^(-0.25);     % Ionization coefficient [1/m], eventually replace this with a lookup table

    f1 = mu(5)*q(1:2) + mu_e*[Ex, Ey].*u(1);    % First term is negative because q=-grad(u)
    f2 = mu(6)*q(3:4) + mu(3)*[Ex, Ey].*u(2);
    f3 = mu(7)*q(5:6) + mu(4)*[Ex, Ey].*u(3);

    f = r*[f1(:) f2(:) f3(:)];
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = sym([0, 0, 0]);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    % pde fluxes
    f = flux(u, q, w, v, x, t, mu, eta);    
    
    % numerical flux
    fh = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u - uhat);

    ne = u(1);
    np = u(2);
    nn = u(3);
    Ex = w(2);       % CHECK THIS! Check the dimension of w (should be 3x1?) and make sure the sign is correct
    Ey = w(3);
    
    mu_e = 1.9163*((Ex^2 + Ey^2)^0.5)^(-0.25);     % Ionization coefficient [1/m]
    mu_p = mu(3);
    mu_n = mu(4);
    gamma = mu(15);
    
    % inviscid fluxes
    fix = [-ne*Ex (mu_p/mu_e)*np*Ex -(mu_n/mu_e)*nn*Ex]; fix = fix(:);     % Nondimensionalized
    fiy = [-ne*Ey (mu_p/mu_e)*np*Ey -(mu_n/mu_e)*nn*Ey]; fiy = fiy(:);
    fih = fix*n(1) + fiy*n(2); %+ tau*(u - uhat); check
    
    % boundary flux on the needle tip - Chen boundary 1
    fb1 = fh;
    fb1(1) = -gamma*np*(((mu_p/mu_e)*Ex)^2 + ((mu_p/mu_e)*Ey)^2)^0.5;            % Nondimensionalized   check sign here (-n)
    fb1(2) = 0;

    % axis symmetric boundary condition   - axisymmetric BC means the diffusive flux=0
    fb2 = [0; 0; 0]; 
    
    % Chen bdry 3
    En = Ex*n(1) + Ey*n(2);
    signEn = tanh(1e3*(-En));    
    alpha = 0.5 + 0.5*signEn;   % Alpha=1 when (-En) is positive: 
    fb3 = alpha*fh + (1-alpha)*0;        % Check for this sign being flipped in the switch
    
    % Grounded boundary -> 4 in paper
    fb4 = fh;
    fb4(1) = 0;
    fb4(3) = 0;
        
    fb = [fb1 fb2 fb3 fb4 fb3];     % Note: For eqns 1-3, BCs 3,5&6 are the same (open boundary)
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    E_bd = mu(16);
    r_tip = mu(17);
    Ua = mu(14);

    % needle tip - boundary 1
    ub1 = u;
    ub1(3) = 0;
    
    % axis symmetric boundary 2
    ub2 = u;
    
    % grounded boundary
    Ex = w(2);       % CHECK THIS! Check the dimension of w (should be 3x1?) and make sure the sign is correct
    Ey = w(3);
    En = Ex*n(1) + Ey*n(2);
    signEn = tanh(1e3*(-En));    
    alpha = 0.5 + 0.5*signEn;
    ub3 = alpha*[0;0;0] + (1-alpha)*0;
   
    ub4 = u;
    ub4(2) = 0;    
    
    ub = [ub1 ub2 ub3 ub4 ub3];     % Note: For eqns 1-3, BCs 3,5&6 are the same (open boundary)
end

function u0 = initu(x, mu, eta)
    x1 = x(1);
    y1 = x(2);
    sigx = 0.01;
    sigy = 0.01;
    x0 = 0.0;
    y0 = 0.0;

    g1 = exp(-0.5*( (x1-x0)^2/sigx^2 + (y1-y0)^2/sigy^2) );
    % g2 = exp(-0.5*( (x1-x0)^2/sigx^2 + (y1-y0)^2/sigy^2) );
    % g3 = exp(-0.5*( (x1-x0)^2/sigx^2 + (y1-y0)^2/sigy^2) );

    u0 = [g1]; % Initialize all to the same profile
end

