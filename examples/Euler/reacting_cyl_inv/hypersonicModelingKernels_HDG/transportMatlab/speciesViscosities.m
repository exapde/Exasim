function mu_out = speciesViscosities(T, species_structs)
    ns = length(species_structs);
    mu_out = zeros(ns,1,class(T));
    for i = 1:ns
        mu_out(i) = evalCurveFit_guptaMu(T, species_structs{i});
    end
end
