function [val, A] = crs_assembly(Ae, elcon, f2e, face, row_ptr, col_ind, ncu, npf, nfe)
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

nblocks = row_ptr(end);
val = zeros(ncu, npf, ncu, npf, nblocks);

nf = length(face);
A = zeros(ncu*npf*nf, ncu*npf*nf);
for i = 1:nf % loop over each local face i
    fi = face(i);    
    ie1 = f2e(1,fi);
    il1 = f2e(2,fi);        
    ie2 = f2e(3,fi);   
    row_start = row_ptr(i) + 1;
    row_end = row_ptr(i+1);    
    in1 = elcon(:,il1,ie1) - npf*(fi-1);   
    if (ie2>0)
      il2 = f2e(4,fi);     
      in2 = elcon(:,il2,ie2) - npf*(fi-1);
      val(:,:,:,:,row_start) = Ae(:,in1,il1,:,in1,il1,ie1) + Ae(:,in2,il2,:,in2,il2,ie2);    
    else
      val(:,:,:,:,row_start) = Ae(:,in1,il1,:,in1,il1,ie1);    
    end    
    I = (ncu*npf*(i-1)+1):(ncu*npf*i);
    A(I,I) = reshape(val(:,:,:,:,row_start), [ncu*npf ncu*npf]);
    for k = (row_start+1):row_end
        j = col_ind(k);  % neighboring face j                
        fj = face(j); 
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
        val(:,:,:,:,k) = Ae(:,in,m1,:,jn,m2,me);    
        J = (ncu*npf*(j-1)+1):(ncu*npf*j);
        A(I,J) = reshape(val(:,:,:,:,k), [ncu*npf ncu*npf]);
    end    
end

val = reshape(val, [ncu*npf ncu*npf nblocks]);

% nf = length(row_ptr) - 1;
% A = zeros(ncu*npf*nf, ncu*npf*nf);
% for i = 1:nf % loop over each local face i
%     ie1 = f2e(1,i);
%     il1 = f2e(2,i);        
%     ie2 = f2e(3,i);   
%     row_start = row_ptr(i) + 1;
%     row_end = row_ptr(i+1);    
%     e1 = elem(ie1);
%     fi = face(i);
%     in1 = elcon(:,il1,e1) - npf*(fi-1);   
%     if (ie2>0)
%       il2 = f2e(4,i);     
%       e2 = elem(ie2);
%       in2 = elcon(:,il2,e2) - npf*(fi-1);
%       val(:,:,:,:,row_start) = Ae(:,in1,il1,:,in1,il1,e1) + Ae(:,in2,il2,:,in2,il2,e2);    
%     else
%       val(:,:,:,:,row_start) = Ae(:,in1,il1,:,in1,il1,e1);    
%     end    
%     I = (ncu*npf*(i-1)+1):(ncu*npf*i);
%     A(I,I) = reshape(val(:,:,:,:,row_start), [ncu*npf ncu*npf]);
%     for k = (row_start+1):row_end
%         j = col_ind(k);  % neighboring face j                
%         je1 = f2e(1,j);       
%         je2 = f2e(3,j);               
%         if (je1 == ie1)
%           me = ie1;                    
%           m1 = il1;
%           m2 = f2e(2,j);          
%         elseif (je1 == ie2)
%           me = ie2;                    
%           m1 = il2;
%           m2 = f2e(2,j);    
%         elseif (je2 == ie1)
%           me = ie1;                    
%           m1 = il1;
%           m2 = f2e(4,j);          
%         elseif (je2 == ie2)
%           me = ie2;                    
%           m1 = il2;
%           m2 = f2e(4,j);             
%         end        
%         e = elem(me);        
%         fj = face(j); 
%         in = elcon(:,m1,e) - npf*(fi-1);     
%         jn = elcon(:,m2,e) - npf*(fj-1);     
%         val(:,:,:,:,k) = Ae(:,in,m1,:,jn,m2,e);    
%         J = (ncu*npf*(j-1)+1):(ncu*npf*j);
%         A(I,J) = reshape(val(:,:,:,:,k), [ncu*npf ncu*npf]);
%     end    
% end
% 
% 
% 
