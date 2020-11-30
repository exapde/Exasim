function [nm,nb] = mkelemblocks(ne,ns)

if nargin<2
    ns = 512; % default number of elements per block
end

if ne==1
    nm = [1; 1];
    nb = 1;
    return;
end

ns = min(ns,ne);
nb = floor(ne/ns);  % number of blocks
na = round(ne/nb); % number of elements per block
nk = 1:na:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];
if (nm(2,end)-nm(1,end))<na/2
    nm(2,end-1) = nm(2,end);
    nm(:,end) = [];
end

nb = size(nm,2);
while min(nm(2,:)-nm(1,:))<ns/2-1 || max(nm(2,:)-nm(1,:))>ns
    nb = nb+1;
    na = round(ne/nb); % number of elements per block
    nk = 1:na:ne;
    nm = [nk(1:end); [nk(2:end)-1,ne]];  
    if (nm(2,end)-nm(1,end))<na/2
        nm(2,end-1) = nm(2,end);
        nm(:,end) = [];
    end        
end
nb = size(nm,2);

if nm(2,end)~=ne
    error('something wrong');
end
if max(nm(2,:)-nm(1,:))>ns
    error('something wrong');
end
if min(nm(2,:)-nm(1,:))<ns/2-1
    error('something wrong');
end
