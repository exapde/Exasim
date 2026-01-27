function out = f_T(T, rho_tilde, rhoe, alpha, species_thermo_structs)
    RU = 8.314471468617452;
    out = alpha;
    % HoverRT = h_n(T)
    for i = 1:5
        out = out + rho_tilde(i) * enthalpyEval(species_thermo_structs{i}, T);
        % out = out + rho_tilde(i) * HoverRT(i)
    end
    % out = out + sum(rho_tilde .* HoverRT)
    % sum += rho_tilde[5] * enthalpyEval2int(tmp[5], T)
    out = out * T;
    out = out - rhoe / RU;
end