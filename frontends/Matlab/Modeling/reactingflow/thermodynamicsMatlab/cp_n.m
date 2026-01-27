function out = cp_n(T, species_thermo_structs)
    out = zeros(1, 5, class(T));
    for i = 1:5
        out(i) = cpEval(species_thermo_structs{i}, T);
    end
end