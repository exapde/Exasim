function [A1, A2, A3, B1, B2, B3, C1, C2, C3, D] = matrix_assembly(Ae, elcon, f2e, face, row_ptr, col_ind, color, count1, count2, count3, count4, ncu, npf, nfe)
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

A1 = zeros(ncu, npf, 1, ncu, npf, 1, nb, count1);
A2 = zeros(ncu, npf, 2, ncu, npf, 2, nb, count2);
A3 = zeros(ncu, npf, 3, ncu, npf, 3, nb, count3);
B1 = zeros(ncu, npf, 1, ncu, npf, nfe-1, nb, count1);
B2 = zeros(ncu, npf, 2, ncu, npf, nfe-2, nb, count2);
B3 = zeros(ncu, npf, 3, ncu, npf, nfe-3, nb, count3);
C1 = zeros(ncu, npf, nfe-1, ncu, npf, 1, nb, count1);
C2 = zeros(ncu, npf, nfe-2, ncu, npf, 2, nb, count2);
C3 = zeros(ncu, npf, nfe-3, ncu, npf, 3, nb, count3);
D = zeros(ncu, npf, 1, ncu, npf, 1, nb, count4);

for n = 1:nb
  for i = 1:nf % loop over each local face i
      fi = face(n, i);    
      ie1 = f2e(1,fi);
      il1 = f2e(2,fi);        
      ie2 = f2e(3,fi);   
      row_start = row_ptr(i) + 1;
      row_end = row_ptr(i+1);    
      in1 = elcon(:,il1,ie1) - npf*(fi-1);   
      if (ie2>0)
        il2 = f2e(4,fi);     
        in2 = elcon(:,il2,ie2) - npf*(fi-1);
        val = Ae(:,in1,il1,:,in1,il1,ie1) + Ae(:,in2,il2,:,in2,il2,ie2);    
      else
        val = Ae(:,in1,il1,:,in1,il1,ie1);    
      end 
      c = color(1:4,row_start);
      if c(1)==1
        A1(:,:,1,:,:,1,n,c(2)) = val;      
      elseif c(1)==2
        A2(:,:,c(3),:,:,c(4),n,c(2)) = val;      
      elseif c(1)==3
        A3(:,:,c(3),:,:,c(4),n,c(2)) = val;        
      elseif c(1)==10
        D(:,:,1,:,:,1,n,c(3)) = val;          
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
          val = Ae(:,in,m1,:,jn,m2,me);    
          c = color(1:4,k);
          if c(1)==2
            A2(:,:,c(3),:,:,c(4),n,c(2)) = val;      
          elseif c(1)==3
            A3(:,:,c(3),:,:,c(4),n,c(2)) = val;        
          elseif c(1)==4
            B1(:,:,c(3),:,:,c(4),n,c(2)) = val;          
          elseif c(1)==5
            B2(:,:,c(3),:,:,c(4),n,c(2)) = val;            
          elseif c(1)==6
            B3(:,:,c(3),:,:,c(4),n,c(2)) = val;            
          elseif c(1)==7
            C1(:,:,c(3),:,:,c(4),n,c(2)) = val;            
          elseif c(1)==8
            C2(:,:,c(3),:,:,c(4),n,c(2)) = val;              
          elseif c(1)==9
            C3(:,:,c(3),:,:,c(4),n,c(2)) = val;             
          elseif c(1)==10
            D(:,:,1,:,:,1,n,c(3)) = val;          
          end
      end    
  end
end

A1 = reshape(A1, [ncu*npf ncu*npf nb count1]);
A2 = reshape(A2, [ncu*npf*2 ncu*npf*2 nb count2]);
A3 = reshape(A3, [ncu*npf*3 ncu*npf*3 nb count3]);
B1 = reshape(B1, [ncu*npf ncu*npf*(nfe-1) nb count1]);
B2 = reshape(B2, [ncu*npf*2 ncu*npf*(nfe-2) nb count2]);
B3 = reshape(B3, [ncu*npf*3 ncu*npf*(nfe-3) nb count3]);
C1 = reshape(C1, [ncu*npf*(nfe-1) ncu*npf nb count1]);
C2 = reshape(C2, [ncu*npf*(nfe-2) ncu*npf*2 nb count2]);
C3 = reshape(C3, [ncu*npf*(nfe-3) ncu*npf*3 nb count3]);
D = reshape(D, [ncu*npf ncu*npf nb count4]);

