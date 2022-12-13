function avField = avfieldLaplacian2d(u, q, w, v, x, t, mu, eta, nspecies, nd)
    % Mutation outputs
    p = eta(1);
    gam = eta(2);
    gam1 = gam - 1;
    nenergy = 1;
    ncu = nspecies + nd + nenergy;
    % artificial viscosity
    porder = mu(15);

    % regularization parameters for the bulk viscosity
    kb = mu(16);
    sb0   = mu(17); %0.02 
    sbmin = 0.0;
    sbmax = mu(18);% / sqrt(gam*gam - 1.0); %2.5 TODO: unsure that this should be changing 

    % regularization parameters
    alpha = mu(22);
    sigm = mu(23);
    rmin = 0.0;
    Hmin = 1.0e-4;

    % mesh size
    hm = v(2);

    % Get base variables
    rho_i = u(1:nspecies);
    rhou = u(nspecies+1);
    rhov = u(nspecies+nd);
    rhoE = u(nspecies+nd+1);

%     drho_dx_i = -q(1:nspecies) .* dlmax(u(1:nspecies)-rmin,alpha);
    drho_dx_i = -q(1:nspecies);
    drhou_dx  = -q(nspecies+1);
    drhov_dx  = -q(nspecies+2);
%     drhoE_dx  = -q(nspecies+ndim+1);
    drho_dy_i = -q(ncu+1:ncu+nspecies);
    drhou_dy  = -q(ncu+nspecies+1);
    drhov_dy  = -q(ncu+nspecies+2);
%     drhoE_dy  = -q(ncu+nspecies+ndim+1);

    % Regularization 
    rho_i = rmin + lmax(rho_i-rmin,alpha); % need to double check regularizatino here
%     drho_dx_i = drho_dx_i .*  dlmax(u(1:nspecies)-rmin,alpha);
%     drho_dy_i = drho_dy_i .*  dlmax(u(1:nspecies)-rmin,alpha);
    drho_dx = sum(drho_dx_i);
    drho_dy = sum(drho_dy_i);
    % Simple derived working quantities
    rho = sum(rho_i);
    rho_inv = 1./rho;
    uv = rhou.*rho_inv;
    vv = rhov.*rho_inv;
    E = rhoE * rho_inv; %energy
    H = E + p*rho_inv; %enthalpy
    H = Hmin + lmax(H - Hmin, alpha);

    % Critical speed of Sound
    c_star = sqrt((2.*gam1.*H) ./ (gam+1)); %TODO: this is fine for 1 temp but not 2

    % Computing derivatives for the sensors
    du_dx = (drhou_dx - drho_dx.*uv).*rho_inv;
    dv_dy = (drhov_dy - drho_dy.*vv).*rho_inv;
    du_dy = (drhou_dy - drho_dy*uv)*rho_inv;
    dv_dx = (drhov_dx - drho_dx*vv)*rho_inv;

    div_v = (du_dx + dv_dy); %TODO: check signs. Actually this probably is important, we want to make sure it's applied in negative v. pos dilitation
                  %      pretty sure my signs are okay. There is something strange about the code I was given. Trhouly don't understand the signs     
    % limit  divergence and vorticity
    div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);

    DucrosRatio = 1;
    sb = - (hm./porder) .* (div_v./c_star) .* DucrosRatio;
    sb = limiting(sb,sbmin,sbmax,alpha,sb0); % TODO: what should sbmin, sbmax, alpha, and sb0 be 

    % Artificial Bulk viscosity
    avb = (kb.*hm./(porder)) .* (sqrt(uv.^2 + vv.^2) + c_star) .* sb;
    
    % Assign artificial viscosities
    avField(1) = avb;  %  bulk
end