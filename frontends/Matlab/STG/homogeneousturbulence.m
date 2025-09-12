function [fluctamp,pars,Ek,Em,awaveno,dwaveno,waveno] = homogeneousturbulence(gridLength, turbLength, visc, urms, N)

% gridLength: Smallest grid size
% turbLength: Turbulence length scale
% visc: Molecular viscosity
% urms: the root-mean-square of the turbulent velocity fluctuations for a unit freestream velocity
% N: number of wave numbers 


% Initial guess for alpha
alpha = 1.453;

while (1)
    
% turbulent kinetic energy
turbKE = (3/2)*urms*urms;

% turbulence dissipation rate
turbDR = ((0.09)^(3/4))*turbKE^1.5/turbLength;

% Kolmogorov wave number
keta = (turbDR/visc^3)^0.25;

% maximum wavenumber
kmax = keta;

%  the wave number ke corresponding to the most energy containing eddies at the peak in the spectrum.
ke = 9*pi*alpha/(55*turbLength);

while (1)
    % minimum wavenumber
    kmin = min(ke(:))/(5); 

    % generate wave number distribution
    [waveno, ~, dwaveno] = wavenumberdistr(N, kmin, kmax);
    
    % compute the integral length scale
    kke = waveno/ke;
    kke2 = kke.*kke;
    Ek = ((kke2 .* kke)./((1.0 + kke2).^(17/6))).*exp(-2.0*(waveno/keta).^2);
    Em = 0.5*(Ek(2:end) + Ek(1:end-1));
    
    % The integral length scale
    Lambda = (alpha*pi/(2*ke*ke))*sum(Em.*dwaveno);
    if abs(Lambda-turbLength)<(1e-3*turbLength)                    
        break;
    else                        
        ke = (Lambda/turbLength)*ke;        
    end        

%     % compute the integral dissipation rate
%     kke = waveno/ke;
%     kke2 = kke.*kke;
%     Ek = ((waveno.*waveno.*kke2 .* kke2)./((1.0 + kke2).^(17/6))).*exp(-2.0*(waveno/keta).^2);
%     Em = 0.5*(Ek(2:end) + Ek(1:end-1));
%     
%     % The integral dissipation rate
%     Lambda = (2*visc*alpha*urms*urms/ke)*sum(Em.*dwaveno);
%     if abs(Lambda-turbDR)<(1e-3*turbDR)                    
%         break;
%     else                        
%         ke = (turbDR/Lambda)*ke;        
%     end        

%    ke*turbLength
end

% minimum wavenumber
kmin = min(ke(:))/(5); 

% generate wave number distribution
[waveno, awaveno, dwaveno] = wavenumberdistr(N, kmin, kmax);

% Energy spectrum amplitude
Eamp = alpha*(2*turbKE/3)/ke;

% compute the energy spectrum at the mid wave numbers
Ek = vonkarmanspectrum(waveno, ke, keta, Eamp);
%Em = 0.5*(Ek(2:end) + Ek(1:end-1));
Em = vonkarmanspectrum(awaveno, ke, keta, Eamp);

% The amplitude of each mode
fluctamp = sqrt(Em.*dwaveno);

% the integral kinetic energy
turbk2 = sum(fluctamp.^2);

if abs(turbk2-turbKE)<(1e-3*turbKE)
    %[turbk2 turbKE]
    % parameters
    pars = [alpha turbKE turbDR ke kmin keta Eamp];    
    break;
else
    alpha = (turbKE/turbk2)*alpha;        
end

end


% k = linspace(kmin, kmax, 1e4);
% E = vonkarmanspectrum(k, ke, keta, Eamp);
% figure(1); %clf;
% hold on;
% plot(k,E);
% hold on;
% plot(ke, vonkarmanspectrum(ke, ke, keta, Eamp),'ko');
% plot(sqrt(12/5)*ke, vonkarmanspectrum(sqrt(12/5)*ke, ke, keta, Eamp),'ro');
% ka = 2*pi/turbLength;
% plot(ka, vonkarmanspectrum(ka, ke, keta, Eamp),'rs');
% ka = 2*pi/gridLength;
% plot(ka, vonkarmanspectrum(ka, ke, keta, Eamp),'rd');
% plot(keta, vonkarmanspectrum(keta, ke, keta, Eamp),'r*');

% k = linspace(kmin, kmax, 1e4);
% E = vonkarmanspectrum(k, ke, keta, Eamp);
% figure(2); clf;
% semilogy(k,E);
% hold on;
% semilogy(ke, vonkarmanspectrum(ke, ke, keta, Eamp),'ko');
% semilogy(sqrt(12/5)*ke, vonkarmanspectrum(sqrt(12/5)*ke, ke, keta, Eamp),'ro');
% ka = 2*pi/turbLength;
% semilogy(ka, vonkarmanspectrum(ka, ke, keta, Eamp),'rs');
% ka = 2*pi/gridLength;
% semilogy(ka, vonkarmanspectrum(ka, ke, keta, Eamp),'rd');
% semilogy(keta, vonkarmanspectrum(keta, ke, keta, Eamp),'r*');

k = linspace(kmin, kmax, 1e4);
E = vonkarmanspectrum(k, ke, keta, Eamp);
figure(1); clf; 
%hold on;
loglog(k,E);
hold on;
loglog(ke, vonkarmanspectrum(ke, ke, keta, Eamp),'ko');
loglog(sqrt(12/5)*ke, vonkarmanspectrum(sqrt(12/5)*ke, ke, keta, Eamp),'ro');
ka = 2*pi/turbLength;
[ke ka ka/ke]
loglog(ka, vonkarmanspectrum(ka, ke, keta, Eamp),'rs');
ka = 2*pi/gridLength;
loglog(ka, vonkarmanspectrum(ka, ke, keta, Eamp),'rd');
loglog(keta, vonkarmanspectrum(keta, ke, keta, Eamp),'r*');

% DNS directly solves the Navier-Stokes equations capturing all eddies from the length scale of the grid geometry right down to the Kolmogorov length scales (relating to the smallest eddies). The dx,dy,dz (=dL) of the mesh needs to be small enough to capture eddies down to the Kolmogorov length scale.
% 
% The argument for the cell resolution, (and thus dL) goes something like this:
% 
% computational box of length L
% 
% number of grid points in ONE direction, N
% 
% grid spacing dL
% 
% Kolmogorov length scale, eta
% 
% Molecular viscosity, mu
% 
% Energy dissipation rate, epsilon
% 
% rms turbulent velocity scale, u'
% 
% -------------------
% 
% For a box of length L, the number of points depends on dL:
% 
% N = L / dL --------------(1)
% 
% dL must be small enough to resolve the smallest eddies, which have the length scale eta. Thus dL=eta is the maximum value for dL to capture the smallest eddies without them `dropping through' the grid (idealy dL= 0.5 * eta for better resolution).
% 
% Thus (1) becomes:
% 
% N(min) = L / dL(max) = L / eta --------(2)
% 
% Eta is defined as: eta = ( mu^3 / epsilon )^(1/4) --------(3)
% 
% Epsilon defined as: epsilon = u'^3 / L --------(4)
% 
% Substituting (3) & (4) into (2) gives:
% 
% N = ( u' * L / mu )^(3/4)
% 
% Noting that u'L/mu is a form of Reynolds Number this gives: N^3 = Re^(9/4)
% 
% Knowing the Reynolds number, Re and the geometry size, L enables you to make a rough estimate of dL.
% 
% ---------------------------
% 
% The papers I've seen are concerned with making the cell size near the wall small enough to get the first few mesh points (this is related to pi*eta (or pi*dL) thus the first 3 points) in the viscous sub-layer. This is why they are always mentioning the y+ values of the near wall cells.
% 
% (3)&(4) See standard turbulence textbooks for definitions.
% 
% For papers discussing this matter:
% 
% Eggels et al JOURNAL OF FLUID MECH. VOL 268 pp175-209 (page 179 specifically). A paper on DNS of turbulent pipe flow
% 
% Kim, Moin, Moser JOURNAL OF FLUID MECH. VOL 177 pp133-166 (p135, there is also a reference to Moser&Moin 1984 (internal Stanford report)). A paper on DNS of channel flow
% 
% An LES computational grid only needs a dL small enough to resolve the large scale flow structures (e.g a recirculation bubble). Any structures smaller than this are passed on to the subgrid scale (SGS) model.
% 
% Hope this helps, Denver
