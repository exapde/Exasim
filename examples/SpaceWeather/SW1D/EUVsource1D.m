function s_EUV = EUVsource1D(u, x, t, mu, eta)
    r = x(1);
    
    gam = mu(1);
    gam1 = gam - 1.0;
    Minf = mu(4);
    M2 = Minf^2;
    
    Fr2 = mu(5);
    omega = mu(6);
    Q0 = mu(12);
    M0 = mu(13);
    
    R0 = mu(10);
    R1 = mu(11);
    Ldim = mu(14);
    
    longitude = mu(17)*pi/180;
    latitude = mu(18)*pi/180;
    declinationSun = mu(19)*pi/180;
    
    %% computation of angles
    %define local time
    long_offset = omega*t;
    localTime = longitude + long_offset;
    cosChi = sin(declinationSun)*sin(latitude) + cos(declinationSun)*cos(latitude)*cos(localTime);
%     cosChi = cos(localTime);
    
    absSinChi = sqrt(1-cosChi^2);
    
    %Computation F10.7 (let's assume it constant at first, the variation is at another scale)
    F10p7 = 100;
    F10p7_81 = 100;
    F10p7_mean = 0.5*(F10p7 + F10p7_81);
    
    rtilde = u(1);
    rho = exp(rtilde);
    T = u(3)/sqrt(rho);

    % Quantities
    gravity = Fr2*(R0^2/(r^2));
    H = T/(gam*M2*gravity);
    
    %Chapman integral 
    Rp = rho*H;
    Xp = r/H;
    y = sqrt(Xp/2)*abs(cosChi);
    
    Ierf = 0.5*(1+tanh(1000*(8-y)));
    a_erf = 1.06069630;
    b_erf = 0.55643831;
    c_erf = 1.06198960;
    d_erf = 1.72456090;
    f_erf = 0.56498823;
    g_erf = 0.06651874;

    erfcy = Ierf*(a_erf + b_erf*y)/(c_erf + d_erf*y + y*y) + (1-Ierf)*f_erf/(g_erf + y);
    
%     erfy = erf(y);
%     erfcy = exp(y^2)*(1-erfy);
    
    IcosChi = 0.5*(1 + tanh(100000*cosChi));
    IsinChi = 0.5*(1 + tanh(100000*(r*absSinChi - R0)));
    
    alpha1 = Rp*erfcy*sqrt(0.5*pi*Xp);
    auxXp = (1-IcosChi)*IsinChi*Xp*(1-absSinChi);
    Rg = rho*H*exp(auxXp);
    alpha2 = (2*Rg - Rp*erfcy)*sqrt(0.5*pi*Xp);
    
    alpha = IcosChi*alpha1 + (1-IcosChi)*(IsinChi*alpha2 + (1-IsinChi)*1e32);
    
    Q = 0;
    for iWave = 1:37
        lambda = eta(iWave);
        crossSection = eta(37+iWave);
        AFAC = eta(2*37+iWave);
        F74113 = eta(3*37+iWave);
        
        tau = M0*crossSection*alpha;
        
        slope0 = 1 + AFAC*(F10p7_mean-80);
        Islope = 0.5*(1+tanh(1000*(slope0-0.8)));
        slopeIntensity =  slope0*Islope + 0.8*(1-Islope);
        Intensity0 = F74113*slopeIntensity;
        Intensity = Intensity0*exp(-tau);
        
        Q = Q + crossSection*Intensity/lambda;
    end
    
    eff = mu(15);
    
    s_EUV = gam*gam1*M2*eff*Q0*Q;