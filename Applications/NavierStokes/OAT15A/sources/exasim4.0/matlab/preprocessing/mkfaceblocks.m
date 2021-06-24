function [nm,nb] = mkfaceblocks(mf,bcm,ns)

if nargin<3
    ns = 4096; % default number of faces per block
end

nm = [];
for i = 1:length(mf)-1
    nf = mf(i+1)-mf(i);    
    nmf = mkelemblocks(nf,ns);
    tm = mf(i)+nmf;
    tm(3,:) = bcm(i);      
    nm = [nm tm];    
end
nb = size(nm,2);

