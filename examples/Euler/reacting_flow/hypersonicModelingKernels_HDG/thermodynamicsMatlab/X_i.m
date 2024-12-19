function out = X_i(rho_i, Mw)
    conc = conc_i(reshape(rho_i, size(Mw)), Mw);
    out =  conc ./ sum(conc);
end