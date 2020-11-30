function mkfaceblocks(mf,bcm,ns::Int=4096)

i = 1;
nf = mf[i+1]-mf[i];
nmf,nb = mkelemblocks(nf,ns);
tm = vcat(mf[i] .+ nmf, bcm[i]*ones(Int,1,size(nmf,2)));
nm = tm;

for i = 2:length(mf)-1
    nf = mf[i+1]-mf[i];
    nmf,nb = mkelemblocks(nf,ns);
    tm = vcat(mf[i] .+ nmf, bcm[i]*ones(Int,1,size(nmf,2)));
    nm = [nm tm];    
end
nb = size(nm,2);

return nm,nb

end
