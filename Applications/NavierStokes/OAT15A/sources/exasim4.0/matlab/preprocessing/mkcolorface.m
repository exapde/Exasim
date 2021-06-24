function [fcolor,ncolor] = mkcolorface(facecon,npe,f1,f2)

if ndims(facecon)==3
    [npf,~,nf] = size(facecon);
    facecon = reshape(permute(facecon,[2 1 3]),[2 npf*nf]);
end

nf = f2-f1+1;
ndf = npf*nf;
uind = zeros(ndf,2);
for i=1:ndf
    m = npf*(f1-1)+i;
    k1 = facecon(1,m); m1 = rem(k1-1,npe)+1; n1 = (k1-m1)/npe+1;
    k2 = facecon(2,m); m2 = rem(k2-1,npe)+1; n2 = (k2-m2)/npe+1;      
    uind(i,1) = m1+(n1-1)*npe; 
    uind(i,2) = m2+(n2-1)*npe;         
end

uind = reshape(uind,[npf nf 2]);
[fcolor,ncolor] = colorface(uind);
f = f1:1:f2;
fcolor = f(fcolor);


function [fcolor,ncolor] = colorface(uind)

nf = size(uind,2);
c = 1;
ncolor = [0 0];
while (1)    
    if c == 1
        color = setcolor(uind);
        fcolor = color;
        ind = fcolor;
        ncolor(c+1) = length(fcolor);
    else
        ind = setdiff((1:nf), fcolor);        
        if ~isempty(ind)
            color = setcolor(uind(:,ind,:));
            fcolor = [fcolor ind(color)];
            ncolor(c+1) = length(fcolor);
        end
    end    
    if isempty(ind)
        break;
    else
        c = c + 1;
    end
end

for i = 1:c-1
    for j = 1:2
        tm = uind(:,fcolor(ncolor(i)+1:ncolor(i+1)),j);
        if length(unique(tm(:))) ~= numel(tm)
            error('something wrong');
        end
    end
end


function [color,all] = setcolor(uind)
% find all faces that have the same color

nf = size(uind,2);
all = [uind(:,1,1) uind(:,1,2)];
j = 1;
color = 1; % list of faces
for i = 2:nf
    ui = [uind(:,i,1) uind(:,i,2)];
    if any(ismember(ui(:),all(:)))==0 % if all elements of ui is not in uall
        j = j + 1;
        all = [all ui];
        color(j) = i;         
    end
end

