function sortcolumns(a)

n = size(a,2);
b = sortslices([a; reshape(collect(1:n),1,n)], dims=2);
ind = Int.(b[end,:]);
b = b[1:end-1,:];

return b,ind

end
