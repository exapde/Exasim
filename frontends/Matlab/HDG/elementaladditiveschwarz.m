function BE = elementaladditiveschwarz(AE, f2e, elcon)

[npf, nfe, ne] = size(elcon);
ndf = npf*nfe;
nf = size(f2e,2);
n = size(AE,1);
nch = n/ndf;
ncf = nch*npf;

BE  = reshape(AE,[nch npf nfe nch npf nfe ne]); 
for i = 1:nf  
  e2 = f2e(3,i); 
  if e2>0 
      e1 = f2e(1,i);
      i1 = f2e(2,i);  
      i2 = f2e(4,i);        
      j2 = elcon(:,i2,e2) - npf*(i-1);   
      tm = BE(:,:,i1,:,:,i1,e1) + BE(:,j2,i2,:,j2,i2,e2);            
      BE(:,:,i1,:,:,i1,e1) = tm;      
      BE(:,j2,i2,:,j2,i2,e2) = tm;            
  end
end

BE  = reshape(BE,[ncf*nfe ncf*nfe ne]); 
for i = 1:ne
  BE(:,:,i) = inv(BE(:,:,i));
end


