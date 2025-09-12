function mkf2e(t,elemtype,nd)

ne = size(t,2);

# get local faces of an element
face = getelemface(nd,elemtype);
nvf,nfe = size(face);

# sort faces for all elements to make two identical faces have the same numbering order
N = nfe*ne;
tf = reshape(t[face[:],:],nvf, N);
tf = sort(tf, dims=1);

tf, jx = sortcolumns(tf); # do this so that two identical faces are next each other
jx = Int.(jx);
dx = tf[:,2:end]-tf[:,1:end-1]; # do this to find two identical faces
dx = sum(dx.*dx,dims=1);           # identical faces corresponds to zero entries
in1 = findall(dx[:].==0);           # interior global faces connected to 1st elements (1st identical faces)
in2 = in1.+1;                 # interior global faces connected to 2nd elements (2nd identical faces)

in0 = setdiff(collect(1:length(jx)), unique([in1[:]; in2[:]])); # boundary global faces connected to 1st elements

nf = length(in0)+length(in1);
f2e = zeros(Int,4,nf);   # allocate memory
e2e = zeros(Int,nfe,ne); # allocate memory

# interior faces
e1 = Int.(ceil.(Float64.(jx[in1])./Float64(nfe)));      # 1st elements
l1 = jx[in1] .- (e1.-1).*nfe;   # 1st local faces
e2 = Int.(ceil.(Float64.(jx[in2])./Float64(nfe)));      # 2nd elements
l2 = jx[in2] .- (e2.-1).*nfe;   # 2nd local faces
g = 1:length(in1);           # indices for interior faces
f2e[1,g] = e1;
f2e[2,g] = l1;
f2e[3,g] = e2;
f2e[4,g] = l2;
for i = 1:length(e1)
    e2e[l1[i],e1[i]] = e2[i];
    e2e[l2[i],e2[i]] = e1[i];
end

# boundary faces
e1 = Int.(ceil.(Float64.(jx[in0])./Float64(nfe)));      # 1st elements
l1 = jx[in0] .- (e1.-1).*nfe;   # 1st local faces
g = (length(in1)+1):nf;      # indices for boundary faces
f2e[1,g] = e1;
f2e[2,g] = l1;

return f2e, e2e

end
