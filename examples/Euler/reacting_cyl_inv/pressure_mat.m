function p = pressure(T, r_i, Mw)
    RU = 8.314471468617452;
    p = T .* sum(r_i ./ Mw , 2) * RU;
end
