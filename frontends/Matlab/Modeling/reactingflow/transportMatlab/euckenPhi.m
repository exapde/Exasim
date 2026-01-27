function phi = euckenPhi(mu_i, Mw, X_i)
    ns = length(mu_i);
    phi = zeros(ns,1,class(mu_i(1)));
    for i = 1:ns
        for j = 1:ns
            if i == j
                phi(i) = phi(i) + X_i(j);
            else
                mu_ratio = mu_i(i) / mu_i(j);
                M_ratio = Mw(i) / Mw(j);
                tmp = 1.0 + sqrt(mu_ratio / sqrt(M_ratio));
                phi(i) = phi(i) + X_i(j) * tmp^2 / sqrt(8.0 *(1.0 + M_ratio));
                % phi(i) += X_i(j) * (1.0 + sqrt(mu_ratio) * (1.0 / M_ratio)^(0.25))^2 / sqrt(8.0 * (1.0 + M_ratio))
            end
        end
    end
end
