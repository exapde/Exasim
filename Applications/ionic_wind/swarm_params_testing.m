E = logspace(3,7,100);
N = 2.686780111798444e+25;

swarm = get_swarm_params(3e6, N);
disp(swarm)

alpha_vec = zeros(100);
eta_vec = zeros(100);
D_vec = zeros(100);
mue_vec = zeros(100);
mun_vec = zeros(100);


for i=1:size(E,2)
    swarm = get_swarm_params(E(i), N);
    swarm = num2cell(swarm);
    [alpha, eta, beta, D, mue, mup, mun] = swarm{:};
    disp(mup);
    alpha_vec(i) = alpha;
    eta_vec(i) = eta;
    D_vec(i) = D;
    mue_vec(i) = mue;
    mun_vec(i) = mun;
end

% figure;
% plot(E, alpha_vec);
% title('alpha');
% figure;
% plot(E, eta_vec);
% figure;
% plot(E, D_vec);
% figure;
% plot(E, mue_vec);
% figure;
% plot(E, mun_vec);
% disp(beta)


function params = get_swarm_params(normE, N)
    % Using the values in the appendix of this paper: https://iopscience.iop.org/article/10.1088/0022-3727/30/4/017/pdf
    % Note that the reduced E field is converted to V*cm2 for performing the calculations, then the resultant quantity is converted back from cm->m.
    % Expects the reduced E field in V*m2. Should be computed by dividing the L2 norm of the E field by the number density of neutrals.

    cm2m = 0.01;
    m2cm = 100;
    E_N = normE/N*m2cm^2;   % Reduced electric field in V*cm2. Expects E to be the L2 norm of the E field.

    % First ionization coefficient
    if E_N <= 1.5e-15
        alpha = 6.619e-17*exp(-5.593e-15/E_N)*cm2m^2;
    elseif E_N > 1.5e-15
        alpha = 2e-16*exp(-7.248e-15/E_N)*cm2m^2;
    end

    % Two-body attachment coefficient
    if E_N <= 1.05e-15
        eta2 = 6.089e-4*E_N - 2.893e-19*cm2m^2;
    elseif E_N > 1.05e-15
        eta2 = 8.889e-5*E_N + 2.567e-19*cm2m^2;
    end

    % Recombination coefficient
    beta = 2e-7*cm2m^3;

    % Electron mobility
    % Removing the sign term because these formulas were designed for computing the velocity, not the mobility. The mobility always has the same sign.
    if E_N <= 2.6e-17
        mue = -(3.38e4 + 6.87e22*E_N) * cm2m / normE;   % Formula is for the velocity, so need to divide by the E field strength to get the mobility
    elseif ((E_N > 2.6e-17) && (E_N <= 1e-16))
        mue = -(1.63e6 + 7.2973e21*E_N) * cm2m / normE;
    elseif ((E_N > 1e-16) && (E_N <= 2e-15))
        mue = -(1.3e6 + 1.03e22*E_N) * cm2m / normE;
    elseif E_N > 2e-15
        mue = -(7.1e6 + 7.4e21*E_N) * cm2m / normE;
    end
    
    % Negative ion mobility
    if E_N <= 5e-16
        mun = -1.86*cm2m^2;
    elseif E_N > 5e-16
        mun = -2.7*cm2m^2;
    end

    % Positive ion mobility
    % Assuming that P0/P is approximately 1 and that the units for the E field are V/cm
    mup = 2.34*cm2m^2;

    % Electron diffusion coefficient
    % Unit conversion (cm->m) is "baked in" through the mue calculated previously
    D = mue*0.3341e9*E_N^0.54069;

    params = [alpha, eta2, beta, D, mue, mup, mun];
end