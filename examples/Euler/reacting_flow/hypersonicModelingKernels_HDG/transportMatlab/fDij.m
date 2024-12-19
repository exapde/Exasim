function Dijout = fDij(T, gupta_structs)
    ns = 5;
    Dijout = zeros(ns, ns, class(T));
    for i = 1:ns
        for j = 1:ns
            Dijout(i,j) = evalCurveFit_guptaDij(T, gupta_structs{i,j});
        end
    end
end
