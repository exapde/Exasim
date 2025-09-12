function h = h_n(T, species_thermo_structs)
    h = zeros(1,5, class(T));
    for i = 1:5
        h(i) = enthalpyEval(species_thermo_structs{i}, T);
    end
end