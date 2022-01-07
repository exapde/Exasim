function sortrows(a)

n = size(a,1);
b = sortslices([a collect(1:n)], dims=1);
ind = Int.(b[:,end]);
b = b[:,1:end-1];

return b,ind

end
