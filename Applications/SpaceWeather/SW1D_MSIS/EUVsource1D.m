function s_EUV = EUVsource1D(u, x, t, mu, eta)
    %Parameters that need to be defined
    nspecies = 4;
    nWaves = 37;
    
    %Computation itself
    r = x(1);
    
    gam = mu(1);
    gam1 = gam - 1.0;
    
    Fr = mu(4);
    omega = Fr/sqrt(gam);
    K0 = mu(5);
    M0 = mu(6);
    
    R0 = mu(16);
    H0 = mu(18);
    
    latitude = mu(23)*pi/180;
    longitude = mu(24)*pi/180;
    declinationSun0 = mu(8)*pi/180;
    doy = mu(11);
    t0 = mu(21);
    seconds = t*t0;
    Ndays = doy + seconds/86400;
    
    declinationSun = asin(-sin(declinationSun0)*cos(2*pi*(Ndays+9)/365.24 + pi*0.0167*2*pi*(Ndays-3)/365.24));
    
    %% computation of angles
    %define local time
    long_offset = omega*t + 2*pi*doy - 3*pi/4;
    localTime = longitude + long_offset;
    cosChi = sin(declinationSun)*sin(latitude) + cos(declinationSun)*cos(latitude)*cos(localTime);
    
    absSinChi = sqrt(1-cosChi^2);
    
    %Computation F10.7 (let's assume it constant at first, the variation is at another scale)
    F10p7 = mu(9);
    F10p7_81 = mu(10);
    F10p7_mean = 0.5*(F10p7 + F10p7_81);
    
    rtilde = u(1);
    rho = exp(rtilde);
    T = u(3)/sqrt(rho);

    % Quantities
    gravity = (R0^2/(r^2))/gam;
    H = T/(gam*gravity);
    
    %Chapman integral 
    Rp = rho*H;
    Xp = r/H;
    y = sqrt(Xp/2)*abs(cosChi);
    
    Ierf = 0.5*(1+tanh(100*(8-y)));
    a_erf = 1.06069630;
    b_erf = 0.55643831;
    c_erf = 1.06198960;
    d_erf = 1.72456090;
    f_erf = 0.56498823;
    g_erf = 0.06651874;

    erfcy = Ierf*(a_erf + b_erf*y)/(c_erf + d_erf*y + y*y) + (1-Ierf)*f_erf/(g_erf + y);
    
    IcosChi = 0.5*(1 + tanh(100*cosChi));
    IsinChi = 0.5*(1 + tanh(100*(r*absSinChi - R0)));
    
    alpha1 = Rp*erfcy*sqrt(0.5*pi*Xp);
    auxXp = (1-IcosChi)*IsinChi*Xp*(1-absSinChi);
    Rg = rho*H*exp(auxXp);
    alpha2 = (2*Rg - Rp*erfcy)*sqrt(0.5*pi*Xp);
    
    alpha = IcosChi*alpha1 + (1-IcosChi)*(IsinChi*alpha2 + (1-IsinChi)*1e2);
    
    %Compute weighted density compositions (n_i/rho = Chi/mi)

    
    Chi = sym(zeros(nspecies,1));
    Chi(1) = 1.0;
    for iSpecies = 2:nspecies
        coeffsDensity = eta((3+nspecies)*nWaves+4*(iSpecies-2)+1:(3+nspecies)*nWaves+4*(iSpecies-1));
        Chi(iSpecies) = coeffsDensity(1)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*exp(coeffsDensity(4)*(r-R0)*H0);
        Chi(1) = Chi(1) - Chi(iSpecies);
    end
    
    mass = eta((3+nspecies)*nWaves+4*(nspecies-1)+1:(3+nspecies)*nWaves+4*(nspecies-1)+nspecies);
    mw = sym(0.0);
    for iSpecies = 1:nspecies
        mw = mw + mass(iSpecies)*Chi(iSpecies);
    end
    
    %Compute EUV
    s_EUV = sym(0.0);
    lambda = eta(1:nWaves);
    AFAC = eta(nWaves+1:2*nWaves);
    F74113 = eta(2*nWaves+1:3*nWaves);
    
    for iSpecies = 1:nspecies
        crossSection = eta((3+iSpecies-1)*nWaves+1:(3+iSpecies)*nWaves);
        
        tau = M0*Chi(iSpecies)*crossSection*alpha/mass(iSpecies);
        slope0 = 1 + AFAC*(F10p7_mean-80);
        Islope = 0.5*(1+tanh(1000*(slope0-0.8)));
        slopeIntensity =  slope0.*Islope + 0.8*(1-Islope);
        Intensity0 = F74113.*slopeIntensity;
        Intensity = Intensity0.*exp(-tau);
        
        Q = sum(crossSection.*(Intensity./lambda));
%         for iWave = 1:nWaves
%             lambda = eta(iWave);
%             AFAC = eta(nWaves+iWave);
%             F74113 = eta(2*nWaves+iWave);
%             crossSection = eta((3+ispecies-1)*nWaves+iWave);
% 
%             tau = M0*crossSection*alpha;
% 
%             slope0 = 1 + AFAC*(F10p7_mean-80);
%             Islope = 0.5*(1+tanh(1000*(slope0-0.8)));
%             slopeIntensity =  slope0*Islope + 0.8*(1-Islope);
%             Intensity0 = F74113*slopeIntensity;
%             Intensity = Intensity0*exp(-tau);
% 
%             Q = Q + crossSection*Intensity/lambda;
%         end
        s_EUV = s_EUV + Chi(iSpecies)*Q*mw/mass(iSpecies);
    end
    
    eff = mu(7);
%     eff = eff0*(1+2*exp(-(((r-R0)*H0-80000)/100000).^2));
    
    s_EUV = gam*gam1*eff*s_EUV/K0;