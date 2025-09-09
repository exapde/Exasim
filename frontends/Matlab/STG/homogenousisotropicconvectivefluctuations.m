function [un, vn, wn] = homogenousisotropicconvectivefluctuations(xdg, stgdata, U, t)

fluctamp = stgdata(:,2);
awaveno = stgdata(:,1);
randno = stgdata(:,3:end);

M = size(xdg,1);
N = length(fluctamp);

fluctamp = 2*fluctamp;

% % % random numbers 
% randno = randomgen(N);

un = zeros(M,1);
vn = zeros(M,1);
wn = zeros(M,1);
% loop over wavenumbers
for n = 1:N
    % wave number
    kn = awaveno(n);
    
    % random numbers
    phi = randno(n,1); 
    dx = randno(n,2);
    dy = randno(n,3); 
    dz = randno(n,4);
    sigmax  = randno(n,5);
    sigmay  = randno(n,6); 
    sigmaz = randno(n,7);    
    omega = randno(n,8);
                
    % angle
    an = kn*(dx*(xdg(:,1)-U*t) + dy*xdg(:,2) + dz*xdg(:,3)) + phi + omega*t;
        
    % Fourier mode
    bn = fluctamp(n).*cos(an);
    
    % fluctuating velocity field
    un = un + sigmax*bn;
    vn = vn + sigmay*bn;
    wn = wn + sigmaz*bn;        
end


