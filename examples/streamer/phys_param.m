function param = phys_param()
    % Returns all of the simulation physical parameters that don't depend on the electric field

    l_ref = 1e-4;   % m
    mu_ref = 0.05;  % m^2/(V-s)
    E_ref = 3e6;    % V/m
    e_eps0 = 1.80955e-8;    % Quantity e/epsilon0
    phi0 = 18.75e3;     % V
    N0 = 5e18; % 1/m^3
    z0 = 1e-2; % m
    sigma0 = 4e-4;  % m
    n_background = 1e13;    % 1/m^3

    %          1       2      3       4      5   6    7     8       9
    param = [l_ref, mu_ref, E_ref, e_eps0, phi0, N0, z0, sigma0, n_background];
end

% Reference list of physics parameters
% l_ref = param{1};
% mu_ref = param{2};
% E_ref = param{3};
% e_eps0 = param{4};
% phi0 = param{5};
% N0 = param{6};
% z0 = param{7};
% sigma0 = param{8};
% n_background = param{9};
