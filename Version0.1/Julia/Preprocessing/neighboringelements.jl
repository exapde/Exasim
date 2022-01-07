function neighboringelements(t2t, elem)

t2te = t2t[:,elem];
nbelem = sort(unique(t2te[:]));
if nbelem[1] == 0
    nbelem = nbelem[2:end];
end

return nbelem

end
