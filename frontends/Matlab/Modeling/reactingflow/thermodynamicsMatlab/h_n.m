function h = h_n(T, species_thermo_structs)
    h = zeros(5,1, class(T));
    for i = 1:5
        h(i) = enthalpyEval(species_thermo_structs{i}, T);
    end
end