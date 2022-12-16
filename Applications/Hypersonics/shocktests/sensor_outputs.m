function o = sensor_outputs(u, q, w, v, x, t, mu, eta, hm)
    % Mutation outputs
    nspecies = 5;
    nenergy = 1;
    ndim = 1;
    p = eta(1);
    gam = eta(2);
    gam1 = gam - 1;
    
    % artificial viscosity
    porder = mu(15);

    % regularization parameters for the bulk viscosity
    kb = mu(16);
    sb0   = mu(17); %0.02 
    sbmin = 0.0;
    sbmax = mu(18);% / sqrt(gam*gam - 1.0); %2.5 TODO: unsure that this should be changing 

    % regularization parameters
    alpha = 1.0e12;
    rmin = 0.0;
    Hmin = 1.0e-4;

    % mesh size
%     hm = v(2);

    % % Get base variables
    % rvec = u(1:nspecies);
    % ru = u(nspecies+1);
    % % rv = udg(3);
    % rE = u(nspecies+nd+1);

    % rx = sum(-q(1:nspecies));
    % rux = -q(nspecies+1);

    % % Regularization 
    % rvec = rmin + lmax(rvec-rmin,alpha); % need to double check regularizatino here

    dT_drho_i = eta(2*nspecies+3:3*nspecies+2);
    dT_drhoe = eta(3*nspecies+3);

    rho_i = rmin + lmax(u(1:nspecies) - rmin, alpha);
    drho_dx_i = -q(1:nspecies) .* dlmax(u(1:nspecies)-rmin,alpha)
    rho = sum(rho_i);
    drho_dx = sum(drho_dx_i);

    rhou = u(nspecies+1);
    drhou_dx = -q(nspecies+1);
    rhoE = u(nspecies+ndim+1);
    drhoE_dx = -q(nspecies+ndim+1);

    % Some useful derived quantities 
    rho_inv = 1.0 / rho;
    drho_dx_inv = 1.0 / drho_dx;
    uv = rhou * rho_inv; %velocity
    dudx = (drhou_dx - drho_dx*uv)*rho_inv;
    
    % Simple derived working quantities
    uv = rhou.*rho_inv;
    E = rhoE * rho_inv; %energy
    H = E + p*rho_inv; %enthalpy
    H = Hmin + lmax(H - Hmin, alpha);

    % Critical speed of Sound
    c_star = sqrt((2.*gam1.*H) ./ (gam+1)); %TODO: this is fine for 1 temp but not 2

    % Computing derivatives for the sensors
    duv_dx = (drhou_dx - drho_dx.*uv).*rho_inv;
    div_v = (duv_dx); %TODO: check signs. Actually this probably is important, we want to make sure it's applied in negative v. pos dilitation
                  %      pretty sure my signs are okay. There is something strange about the code I was given. Truly don't understand the signs 
    % limit  divergence and vorticity
    sigm = 1e4;
    div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);
    sb = - (hm./porder) .* (div_v./c_star);

    % % Dilatation Sensor sb
    % re = rhoE - rho * (0.5 * uTu); 
    uTu2 = 0.5 * uv * uv;
    duTu2_dx = uv .* duv_dx; 
    dre_drho = -uTu2;
    dre_duTu2 = -rho;
    dre_drhoE = 1.0;
    dre_dx = dre_drho * drho_dx + dre_duTu2 * duTu2_dx + dre_drhoE * drhoE_dx;
    Tx = sum(dT_drho_i .* drho_dx_i) +  dT_drhoe * dre_dx;
    sT = hm./porder .* (Tx)
    
    %%%%%%%% Modify appropriate coefficients 
    % kappa = kappa + avk;
    % viscshear = viscshear + avs;
    % viscbulk = 0 + avb;
    % D_vec = D_vec + avDvec;

    %%%%%%%% Calculation of J_i
    % for ispecies = 1:nspecies
    %     Y_i(ispecies) = rho_i(ispecies) / rho;
    %     dYdx_i(ispecies) = (q(ispecies) * rho - rho_i(ispecies) * drhodx) * drhodxinv^2;
    % end
    Y_i = rho_i ./ rho;
%     dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rho_inv * rho_inv; 
    dY_dx_i = (drho_dx_i - Y_i * drho_dx) * rho_inv; 

    sY = hm./porder .* dY_dx_i
    o(1) = sb;
    o(2) = sT;
    o(2+1:2+nspecies) = sY;

end
