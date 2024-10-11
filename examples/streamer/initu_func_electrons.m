function u0 = initu_func_electrons(dgnodes, param)

    % Physics parameters
    l_ref = param(1);
    n_background = param(9);

    n_background_tilde = n_background*(l_ref^3);

    u0 = n_background_tilde;     % Uniform background density
    % u0 = log(n_background_tilde);     % Uniform background density
end
