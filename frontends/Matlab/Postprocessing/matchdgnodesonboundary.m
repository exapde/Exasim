function [in, im] = matchdgnodesonboundary(xb1, xb2)

x1 = squeeze(sum(xb1,1));
x2 = squeeze(sum(xb2,1));
in = xiny(x1',x2');

npf = size(xb1,1);
nfb = length(in);
im = zeros(npf, nfb);
for i1 = 1:nfb
  i2 = in(i1);
  if i2 > 0
    x1 = xb1(:,:,i1);
    x2 = xb2(:,:,i2);
    im(:,i1) = xiny(x1,x2);    
  end
end
