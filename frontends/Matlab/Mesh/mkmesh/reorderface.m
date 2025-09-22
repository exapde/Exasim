function [f,t2f,mf] = reorderface(f,t2f,bcm)


ind = f(:,end)<0; 
a = -unique(f(ind,end));
a = sort(a);
b = (1:length(bcm));
if length(a)~=length(b)
    error('Something wrong. Check bcm or f');
end
if max(abs(a(:)-b(:)))~=0
    error('Something wrong. Check bcm or f');
end

bcn = unique(bcm);
nbc = length(bcn);
inb = cell(nbc,1);
for i = 1:nbc
    inb{i} = find(bcm==bcn(i));
end

nf = size(f,1);
mp = zeros(nf,1);
mf = zeros(nbc+2,1);

% interior faces
ind = find(f(:,end)>0); 
mp(1:length(ind))=ind;

% boundary faces
n = length(ind);
mf(2) = n;
for i = 1:nbc
    for j = 1:length(inb{i})       
        ind = find(f(:,end)==-inb{i}(j));
        mp((n+1):(n+length(ind)))=ind;
        n = n + length(ind);
        mf(i+2) = n;
    end
end

[nt,nfv] = size(t2f);

f = f(mp,:);
a = zeros(nf,1);
a(mp) = (1:nf)';
t2f = reshape(t2f,nfv*nt,1);
ii = find(t2f > 0);
t2f(ii) = a(t2f(ii));
ii = find(t2f < 0);
t2f(ii) = -a(-t2f(ii));
t2f = reshape(t2f,nt,nfv);

