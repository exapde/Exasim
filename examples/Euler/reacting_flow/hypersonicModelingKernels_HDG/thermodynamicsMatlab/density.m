function out = density(T, P, X, Mw)
    RU = 8.314471468617452;
    density = sum(X .* Mw);
    density = density .* P ./ (RU .* T);
    out =  density;
end
