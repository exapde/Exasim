function out = G_n(T, species_thermo_structs)
    out = zeros(5, 1, class(T));
%     tmp = {Nstruct, Ostruct, NOstruct, N2struct, O2struct}
    for i = 1:5
        out(i) = gibbsEval(species_thermo_structs{i}, T);
    end
end