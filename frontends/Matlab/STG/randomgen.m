function randno = randomgen(N)

% generate random angles
% phi = unifrnd(0,2*pi,N,1);
% alpha = unifrnd(0,2*pi,N,1);
% psi = unifrnd(0,2*pi,N,1);
% angle = unifrnd(0,1,N,1);

phi   = 2*pi*rand(N,1);
alpha = 2*pi*rand(N,1);
psi   = 2*pi*rand(N,1);
angle = rand(N,1);
theta = acos(1.0 - 2.0 * angle);

% wave directions
dx = sin(theta) .* cos(psi);
dy = sin(theta) .* sin(psi);
dz = cos(theta); 

% amplitudes
sigmax = cos(psi) .* cos(theta) .* cos(alpha) - sin(psi) .* sin(alpha);
sigmay = sin(psi) .* cos(theta) .* cos(alpha) + cos(psi) .* sin(alpha);
sigmaz = -sin(theta) .* cos(alpha);

% 
randno = [phi dx dy dz sigmax sigmay sigmaz];

