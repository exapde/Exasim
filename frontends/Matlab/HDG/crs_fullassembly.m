function val = crs_fullassembly(Ae, elcon, f2e, face, row_ptr, col_ind, ncu, npf, nfe)
% Assemble block CRS matrix and return reusable local index map
% Inputs:
%   Ae: [ncu x npf x nfe x ncu x npf x nfe x ne] elemental matrices
% Outputs:
%   row_ptr, col_ind: CRS format
%   val: block values, [K x K x nnz]
%   local_idx_map: [M x M x N] CRS index per (i,j,e)

M = ncu*npf*nfe;
ne = numel(Ae)/(M*M);
Ae = reshape(Ae, [ncu npf nfe ncu npf nfe ne]);
elcon = reshape(elcon, [npf nfe ne]);

[nb, nf] = size(face);
nblocks = row_ptr(end);

val = zeros(ncu, npf, ncu, npf, nb, nblocks);

for i = 1:nf % loop over each local face i
  for n = 1:nb
    fi = face(n,i);    
    ie1 = f2e(1,fi);
    il1 = f2e(2,fi);        
    ie2 = f2e(3,fi);   
    row_start = row_ptr(i) + 1;
    row_end = row_ptr(i+1);    
    in1 = elcon(:,il1,ie1) - npf*(fi-1);   
    if (ie2>0)
      il2 = f2e(4,fi);     
      in2 = elcon(:,il2,ie2) - npf*(fi-1);
      val(:,:,:,:,n,row_start) = Ae(:,in1,il1,:,in1,il1,ie1) + Ae(:,in2,il2,:,in2,il2,ie2);    
    else
      val(:,:,:,:,n,row_start) = Ae(:,in1,il1,:,in1,il1,ie1);    
    end    
    for k = (row_start+1):row_end
        j = col_ind(k);  % neighboring face j                
        fj = face(n,j); 
        je1 = f2e(1,fj);       
        je2 = f2e(3,fj);               
        if (je1 == ie1)
          me = ie1;                    
          m1 = il1;
          m2 = f2e(2,fj);          
        elseif (je1 == ie2)
          me = ie2;                    
          m1 = il2;
          m2 = f2e(2,fj);    
        elseif (je2 == ie1)
          me = ie1;                    
          m1 = il1;
          m2 = f2e(4,fj);          
        elseif (je2 == ie2)
          me = ie2;                    
          m1 = il2;
          m2 = f2e(4,fj);             
        end        
        in = elcon(:,m1,me) - npf*(fi-1);     
        jn = elcon(:,m2,me) - npf*(fj-1);     
        val(:,:,:,:,n,k) = Ae(:,in,m1,:,jn,m2,me);    
    end    
  end
end

val = reshape(val, [ncu*npf ncu*npf nb nblocks]);

