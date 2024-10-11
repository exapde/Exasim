function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
% pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
% pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    m = x(1)*sym([1.0; 1.0; 0]);    % CHECK
end

function f = flux(u, q, w, v, x, t, mu, eta)

    % Read in values from the p vector
    r_tilde = x(1);

    % Read in values from the u vector
    ne_tilde = u(1);
    ni_tilde = u(2);
    phi_tilde = u(3);
    dne_dr_tilde = q(1);
    dni_dr_tilde = q(2);
    Er_tilde = q(3);
    dne_dz_tilde = q(4);
    dni_dz_tilde = q(5);
    Ez_tilde = q(6);

    % Load physics param
    l_ref = mu(1);
    mu_ref = mu(2);
    E_ref = mu(3);
    e_eps0 = mu(4);

    % Compute transport coefficients
    normE_tilde = sqrt(Er_tilde^2 + Ez_tilde^2);
    normE = normE_tilde*E_ref;

    mue_tilde = (2.3987*normE^(-.26))/mu_ref;
    De_tilde = 4.3628e-3*normE^.22 / (l_ref*mu_ref*E_ref);

    fv = [De_tilde.*dne_dr_tilde,0,Er_tilde,...
        De_tilde.*dne_dz_tilde,0,Ez_tilde];       % The negative sign is included in eqns 1 and 3 bc q=-grad(u)

    fi = [-Er_tilde*mue_tilde*ne_tilde,0,0,...
        -Ez_tilde*mue_tilde*ne_tilde,0,0];

    % Multiply flux by r for axisymmetry
    f = r_tilde*(fi + fv);
    f = reshape(f,[3,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)

    % Read in values from the p vector
    r_tilde = x(1);

    % Read in values from the u vector
    ne_tilde = u(1);
    ni_tilde = u(2);
    phi_tilde = u(3);
    dne_dr_tilde = q(1);
    dni_dr_tilde = q(2);
    Er_tilde = q(3);
    dne_dz_tilde = q(4);
    dni_dz_tilde = q(5);
    Ez_tilde = q(6);

    % Load physics param
    l_ref = mu(1);
    mu_ref = mu(2);
    E_ref = mu(3);
    e_eps0 = mu(4);
    ne_star = mu(10);

    % Compute transport coefficients
    normE_tilde = sqrt(Er_tilde^2 + Ez_tilde^2);
    normE = normE_tilde*E_ref;
    
    mue_tilde = (2.3987*normE^(-.26))/mu_ref;
    De_tilde = 4.3628e-3*normE^.22 / (l_ref*mu_ref*E_ref);
    alpha = (1.1944e6+ 4.3666e26/normE^3)*exp(-2.73e7/normE);
    alpha_tile = alpha*l_ref;
    eta_tilde = 340.75*l_ref;
    alpha_bar_tilde = alpha_tile-eta_tilde;
    
    charge_prefactor_tilde = e_eps0/(E_ref*l_ref^2);
    
    se = alpha_bar_tilde*mue_tilde*normE_tilde*ne_tilde;
    sphi = ne_star*charge_prefactor_tilde*(ni_tilde - ne_tilde);    % Note sign flip in source term because the negative sign is absorbed in the diffusive flux in the LHS
    
    s = r_tilde*[se se sphi];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    r = 1;

    % Read in values from the u vector
    ne_tilde = u(1);
    ni_tilde = u(2);
    phi_tilde = u(3);
    dne_dr_tilde = q(1);
    dni_dr_tilde = q(2);
    Er_tilde = q(3);
    dne_dz_tilde = q(4);
    dni_dz_tilde = q(5);
    Ez_tilde = q(6);

    % Load physics param
    l_ref = mu(1);
    E_ref = mu(3);
    phi0 = mu(5);

    fh = sym([0 0 0])';     % Initialize fbou vector
    % BDRY 1: Bottom electrode: species set to homogeneous nuemann and potential set to dirichlet

    fh1 = fh;
    % Inflow for electrons
    fh1(1) = r*((dne_dr_tilde*n(1)+ dne_dz_tilde*n(2)) + tau*(ne_tilde-uhat(1)));

    % Outflow for ions -- extrapolat u = uhat
    fh1(2) = r*tau*(ni_tilde-uhat(2));

    % Potential: grounded dirichlet
    phi_bottom_electrode = 0;
    fh1(3) = r.*tau .*(phi_bottom_electrode-uhat(3));

    % BDRY 2: Right "farfield"
    % Species + potential all have homogeneous neumann
    f = flux(uhat, q, w, v, x, t, mu, eta);
    fh2 = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat); % zero flux for all species
    
    % BDRY 3: Bottom electrode: species set to homogeneous nuemann and potential set to dirichlet

    fh3 = fh;
    % Outflow for electrons
    fh3(1) = r*tau*(ne_tilde-uhat(1));

    % Inflow for ions
    fh3(2) = r*((dni_dr_tilde*n(1)+ dni_dz_tilde*n(2)) + tau*(ni_tilde-uhat(2)));

    % Potential: grounded dirichlet
    phi0_tilde = phi0/(E_ref*l_ref);
    fh3(3) = r*tau*(phi0_tilde-uhat(3));
    
    % BDRY 4: Symmetry boundary: extrapolate m=u or u=uhat
    fh4 = fh;
    fh4(1) = r*tau*(u(1)-uhat(1));
    fh4(2) = r*tau*(u(2)-uhat(2));
    fh4(3) = r*tau*(u(3)-uhat(3));

    fb = [fh1 fh2 fh3 fh4];
end

function u0 = initu(x, mu, eta)
    u0 = sym([0.0, 0.0, 0.0]);     % Will be overriden by initial condition specified in pdeapp
end


