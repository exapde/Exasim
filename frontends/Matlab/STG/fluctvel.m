function velprime = fluctvel(xdg, t,  Uinf, lemax, kcute, waveno, randno, Amean)

M = size(xdg,1);
N = size(waveno,1);

sn = zeros(M,1);
un = zeros(M,1);
vn = zeros(M,1);
wn = zeros(M,1);
% loop over wavenumbers
for n = 1:N
    % wavenumber at mode n
    kn = waveno(n,1);
    dkn = waveno(n,2);
    
    % random numbers
    phi = randno(n,1); 
    dx = randno(n,2);
    dy = randno(n,3); 
    dz = randno(n,4);
    sigmax  = randno(n,5);
    sigmay  = randno(n,6); 
    sigmaz = randno(n,7);
    
    % Kinetic energy 
    En = calculateEk(kn, kcute(:,1), kcute(:,2), kcute(:,3));    
    % wave amplitude
    qn = En*dkn;
    % sum of N wave amplitudes
    sn = sn + qn;
        
    % position 
    rx = (2.0 * pi / (kn * lemax)) * (xdg(:,1) - Uinf*t);
    % angle
    an = kn*(dx*rx + dy*xdg(:,2) + dz*xdg(:,3)) + phi;
    % Fourier mode
    bn = sqrt(qn).*cos(an);
    
    % fluctuating velocity field
    un = un + sigmax*bn;
    vn = vn + sigmay*bn;
    wn = wn + sigmaz*bn;        
end

% apply scaling
scale = 2.0 * sqrt(3.0 / 2.0);
un = scale*un./sqrt(sn);
vn = scale*vn./sqrt(sn);
wn = scale*wn./sqrt(sn);

% fluctuating velocity field
velprime = zeros(M,3);
velprime(:,1) = Amean(:,1).*un;
velprime(:,2) = Amean(:,2).*un + Amean(:,3).*vn;
velprime(:,3) = Amean(:,4).*un + Amean(:,5).*vn + Amean(:,6).*wn;





