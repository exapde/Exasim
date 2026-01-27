function out = mixtureFrozenCvMass(T, Mw, Y, thermo_structs)
    cv = getCvsMass(T, Mw, thermo_structs);
    out =  sum(cv .* reshape(Y, size(cv)));
    % out = dot(cv, Y);
end
