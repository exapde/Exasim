 
function pde = pdemodel
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    pde.sourcew = @sourcew;
    pde.initw = @initw;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    m = sym(ones(nspecies + 2, 1));
    m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)

    nspecies = 5;

    f = sym(zeros(nspecies + 2,1));
    rho_i = sym(zeros(5,1));
    rho = sym(0);

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = u(ispecies); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    rhoE = u(nspecies+2);

    rhoinv = 1.0 / rho;
    u = rhou * rhoinv; %velocity
    E = rhoE*rhoinv; %energy

    % Mutation outputs

    p = w(nspecies + 1); %pressure

    H = E+p*rhoinv; %enthalpy

    % Fluxes
    for ispecies = 1:nspecies
        f(ispecies) = rho_i(ispecies) * u;
    end
    f(nspecies + 1) = rhou * u + p;
    f(nspecies + 2) = rhou * H;
end

function s = source(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;

    s = sym(zeros(nspecies + 2,1));
    omega_i = sym(zeros(nspecies,1));
    rho_i = sym(zeros(5,1));
    rho = 0;

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = u(ispecies); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    rhoE = u(nspecies+2);

    rhoinv = 1/ rho;
    u = rhou * rhoinv; %velocity
    E = rhoE*rhoinv; %energy

    % Mutation outputs
    for ispecies = 1:nspecies
        omega_i(ispecies) = w(ispecies);
    end
    p = w(nspecies + 1); %pressure

    H = E+p*rhoinv; %enthalpy

    % Sources
    for ispecies = 1:nspecies
        s(ispecies) = omega_i(ispecies);
    end
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ub = sym(zeros(nspecies+2, 2));

    rhoL = mu(1:nspecies);
    rhouL = mu(nspecies+1);
    rhoEL = mu(nspecies+2);

    rhoR = mu(nspecies+3:2*nspecies+2);
    rhouR = mu(2*nspecies+3);
    rhoER = mu(2*nspecies+4);

    ub(:,1) = [rhoL(:); rhouL; rhoEL];
    ub(:,2) = [rhoR(:); rhouR; rhoER];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    fb = sym(zeros(nspecies+2, 2));

    f = flux(u, q, w, v, x, t, mu, eta);

    fb(:,1) = f*n(1) + tau(1)*(u-uhat);
    fb(:,2) = f*n(1) + tau(1)*(u-uhat);
end

function u0 = initu(x, mu, eta)
    % To be manually modified
    nspecies = 5;
    u0 = sym(zeros(nspecies+2));

    rhoL = mu(1:nspecies);
    rhouL = mu(nspecies+1);
    rhoEL = mu(nspecies+2);

    rhoR = mu(nspecies+3:2*nspecies+2);
    rhouR = mu(2*nspecies+3);
    rhoER = mu(2*nspecies+4);

    u0 = [rhoL(:) + rhoR(:); rhouL + rhouR; rhoEL + rhoER];
end

function sw = sourcew(u, q, w, v, x, t, mu, eta)
    nspecies = 5;

    sw = sym(zeros(nspecies + 2,1));
    rho_i = sym(zeros(5,1));
    rho = sym(0);

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = u(ispecies); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    rhoE = u(nspecies+2);

    rhoinv = 1.0 / rho;
    u = rhou * rhoinv; %velocity
    E = rhoE*rhoinv; %energy

    % Mutation inputs
    rhoe = rhoE - 0.5 * rho * u^2;
    for i = 1:nspecies
        sw(i) = rho_i(i);
    end
    sw(nspecies + 1) = rhoe;
end

function w0 = initw(x, mu, eta)
    nspecies = 5;
    w0 = sym(zeros(nspecies + 1,1));
end