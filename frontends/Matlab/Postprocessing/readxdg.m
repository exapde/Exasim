function xdg = readxdg(filesol)
    
tm = readbin(filesol);
sz = tm(1);
k1 = 2;
k2 = k1+(sz)-1;
nsize = tm(k1:k2);

k1 = k2+1;
k2 = k1+nsize(1)-1;
ndims = tm(k1:k2);

k1 = k2+1;
k2 = k1+nsize(2)-1;
xdg = tm(k1:k2);

ne = ndims(1);
npe = ndims(4);
ncx = length(xdg(:))/(npe*ne);
xdg = reshape(xdg, [npe ncx ne]);

% if (nsize(3) > 0) && (nargout > 1)
%     k1 = k2+1;
%     k2 = k1+nsize(3)-1;
%     udg = tm(k1:k2);
%     ncu = length(udg(:))/(npe*ne);
%     udg = reshape(udg, [npe ncu ne]);
% else
%     udg = []; 
% end
% 

