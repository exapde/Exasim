function mkelconcg(dgnodes)

# remove duplicate nodes in mesh.p1
ns,dim,nt = size(dgnodes);
A = reshape(permutedims(dgnodes,[1,3,2]),ns*nt,dim);

# B = flip(A,1);
# [~,I]=unique(B,'rows');
# B = B(sort(I),:);
# B = flip(B,1);
# B = unique(A,dims=1);
# b = xiny(A,B,1);

B,~,b = uniquerows(A);
# display(maximum(abs.(b-b1)))
# ii = 100:200;
# display(b[ii]')
# display(b1[ii]')
# #display(b1[31:60]')
# error("here")

# CG mesh
cgnodes = B';
cgelcon = reshape(b,ns, nt);

return cgnodes,cgelcon

end
