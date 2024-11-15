function T = temperature(p, r_i, Mw)
    RU = 8.314471468617452;
    T = p / (RU * sum(r_i ./ Mw));
end