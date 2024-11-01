function out = wilkeMixture(mu_i, X_i, phi_i)
    out =  sum( mu_i .* X_i ./ phi_i);
end
