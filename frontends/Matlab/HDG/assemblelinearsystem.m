function [AE,FE] = assemblelinearsystem(AE, FE, elcon, f, t2f)

[ndf, ne] = size(elcon);
ncu = size(FE,1) / ndf;

if nargin <=3    
  il = zeros(ndf*ncu,ndf*ncu,ne);
  jl = zeros(ndf*ncu,ndf*ncu,ne);
  for i=1:ne
      con = repmat((elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,ndf);
      con = reshape(con,ndf*ncu,1);    
      il(:,:,i) = repmat(con ,1,ndf*ncu);
      jl(:,:,i) = repmat(con',ndf*ncu,1);        
  end

  AE = sparse(reshape(il,ndf*ncu*ndf*ncu*ne,1),reshape(jl,ndf*ncu*ndf*ncu*ne,1),reshape(AE,ndf*ncu*ndf*ncu*ne,1));        
  FE = sparse(reshape(il(:,1,:),ndf*ncu*ne,1),ones(ndf*ncu*ne,1),reshape(FE,ndf*ncu*ne,1));                                                                      
else    
   AE = reshape(AE, [ncu, ndf, ncu, ndf, ne]);
  [AE,FE] = densesystem(AE, FE, f, t2f, elcon);
end

function [K,F] = densesystem(AE, FE, f, t2f, elcon)
% mkBlockMatrix returns 4-dimensional array for the Jacobian matrix 
%  
%   AE    :  5-dimensional array for the elemental matrices
%   FE    :  3-dimensional array for the elemental vectors
%   f     :  Face to element connectivity
%   t2f   :  Element to face connectivity
%
%   K     :  4-dimensional array for the Jacobian matrix
%   F     :  2-dimensional array for the RHS vector
%   f2f   :  Face to face connectivity

nch = size(AE,1);   % number of components of UH
ndf = size(AE,2);   % number of points per face times number of faces per element
ne  = size(AE,5);   % number of elements
nf  = size(f,1);    % number of faces
nfe = size(t2f,2);  % number of faces per element
npf = ndf/nfe;      % number of points per face
ncf = nch*npf;      % number of components of UH times number of points per face
nbf = 2*nfe-1;      % number of neighboring faces

FE  = reshape(FE,[nch npf nfe ne]);
AE  = reshape(AE,[nch npf nfe nch npf nfe ne]);
elcon = reshape(elcon,[npf nfe ne]);


K   = zeros(ncf,ncf,nbf-1,nf);
F   = zeros(ncf,nf);
for i = 1:nf
    fi = f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face
        kf = t2f(fi,:);         % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
        i2 = find(kf(2,:)==i);  % obtain the index of face i in the 2nd element            
                        
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;        
        
        % the first block        
        k = 1;
        M = reshape(AE(:,j1,i1,:,j1,i1,fi(1)) + AE(:,j2,i2,:,j2,i2,fi(2)), [ncf ncf]);        
        F(:,i)   = reshape(FE(:,j1,i1,fi(1)) + FE(:,j2,i2,fi(2)), [ncf 1]);                
        
        for is=1:nfe % loop over each faces of the 1st element
            if is ~= i1  
                k = k + 1;                 
                j = kf(1,is);
                j3 = elcon(:,is,fi(1)) - (j-1)*npf;                
                K(:,:,k-1,i) = reshape(AE(:,j1,i1,:,j3,is,fi(1)), [ncf ncf]);                                                                    
            end
        end
        
        for is=1:nfe % loop over faces of the 2nd element
            if is ~= i2                                                
                k = k + 1;                 
                j = kf(2,is);
                j4 = elcon(:,is,fi(2)) - (j-1)*npf;                
                K(:,:,k-1,i) = reshape(AE(:,j2,i2,:,j4,is,fi(2)), [ncf ncf]);                                                    
            end
        end        
    else % face i is a boundary face
        kf = t2f(fi(1),:);      % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element         
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;                
        
        % the first block        
        k = 1;
        M = reshape(AE(:,j1,i1,:,j1,i1,fi(1)), [ncf ncf]);        
        F(:,i)   = reshape(FE(:,j1,i1,fi(1)), [ncf 1]);                
        
        for is=1:nfe % loop over each faces of the 1st element
            if is ~= i1  
                k = k + 1;          
                j = kf(1,is);
                j3 = elcon(:,is,fi(1)) - (j-1)*npf;
                K(:,:,k-1,i) = reshape(AE(:,j1,i1,:,j3,is,fi(1)), [ncf ncf]);                                                    
            end
        end        
    end    
    
    Mi = inv(M);
    F(:,i) = Mi*F(:,i);      
    for j=1:nbf-1
       K(:,:,j,i) = Mi*K(:,:,j,i);
    end    
end
K = reshape(K,[ncf,ncf*(nbf-1),nf]);


% function [K,F] = densesystem2(AE, FE, f2t, perm)
% % mkBlockMatrix returns 4-dimensional array for the Jacobian matrix 
% %  
% %   AE    :  5-dimensional array for the elemental matrices
% %   FE    :  3-dimensional array for the elemental vectors
% %   f     :  Face to element connectivity
% %   t2f   :  Element to face connectivity
% %
% %   K     :  4-dimensional array for the Jacobian matrix
% %   F     :  2-dimensional array for the RHS vector
% %   f2f   :  Face to face connectivity
% 
% nfe = size(perm,2);  % number of faces per element
% npf = size(perm,1);  % number of points per face
% ndf = npf*nfe;
% nch = size(AE,1)/ndf; % number of components of UH
% ne  = size(AE,3);     % number of elements
% nf  = size(f2t,2);    % number of faces
% ncf = nch*npf;      % number of components of UH times number of points per face
% nbf = 2*nfe-1;      % number of neighboring faces
% 
% FE  = reshape(FE,[nch npf nfe ne]);
% AE  = reshape(AE,[nch npf nfe nch npf nfe ne]);
% 
% K   = zeros(ncf,ncf,nbf-1,nf);
% F   = zeros(ncf,nf);
% for i = 1:nf
%     fi = f2t(:,i); % obtain two elements sharing the same face i      
%     if fi(3)>0           % face i is an interior face
%         %kf = t2f(fi,:);         % obtain neighboring faces 
%         e1 = fi(1);
%         e2 = fi(3);
%         l1 = fi(2);  % obtain the index of face i in the 1st element
%         l2 = fi(4);  % obtain the index of face i in the 2nd element                                    
%         
%         % the first block        
%         k = 1;
%         M = reshape(AE(:,:,l1,:,:,l1,e1) + AE(:,j2,l2,:,j2,l2,e2), [ncf ncf]);        
%         F(:,i)   = reshape(FE(:,j1,l1,e1) + FE(:,j2,l2,e2), [ncf 1]);                
%         
%         for is=1:nfe % loop over each faces of the 1st element
%             if is ~= l1  
%                 k = k + 1;                 
%                 j = kf(1,is);
%                 j3 = elcon(:,is,e1) - (j-1)*npf;                
%                 K(:,:,k-1,i) = reshape(AE(:,j1,l1,:,j3,is,e1), [ncf ncf]);                                                                    
%             end
%         end
%         
%         for is=1:nfe % loop over faces of the 2nd element
%             if is ~= l2                                                
%                 k = k + 1;                 
%                 j = kf(2,is);
%                 j4 = elcon(:,is,e2) - (j-1)*npf;                
%                 K(:,:,k-1,i) = reshape(AE(:,j2,l2,:,j4,is,e2), [ncf ncf]);                                                    
%             end
%         end        
%     else % face i is a boundary face
%         kf = t2f(e1,:);      % obtain neighboring faces 
%         l1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element         
%         j1 = elcon(:,l1,e1) - (i-1)*npf;                
%         
%         % the first block        
%         k = 1;
%         M = reshape(AE(:,j1,l1,:,j1,l1,e1), [ncf ncf]);        
%         F(:,i)   = reshape(FE(:,j1,l1,e1), [ncf 1]);                
%         
%         for is=1:nfe % loop over each faces of the 1st element
%             if is ~= l1  
%                 k = k + 1;          
%                 j = kf(1,is);
%                 j3 = elcon(:,is,e1) - (j-1)*npf;
%                 K(:,:,k-1,i) = reshape(AE(:,j1,l1,:,j3,is,e1), [ncf ncf]);                                                    
%             end
%         end        
%     end    
%     
%     Mi = inv(M);
%     F(:,i) = Mi*F(:,i);      
%     for j=1:nbf-1
%        K(:,:,j,i) = Mi*K(:,:,j,i);
%     end    
% end
% K = reshape(K,[ncf,ncf*(nbf-1),nf]);
