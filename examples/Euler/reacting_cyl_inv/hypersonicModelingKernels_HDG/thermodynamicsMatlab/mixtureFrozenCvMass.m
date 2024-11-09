function out = mixtureFrozenCvMass(T, Mw, Y, thermo_structs)
    cv = getCvsMass(T, Mw, thermo_structs);
    % out =  sum(cv .* Y);
    out = dot(cv, Y);
end
