function [x, flags, iter, rev] = hdggmres(AE, b, BE, f2e, elcon, x, restart, tol, maxit)

% check the number of input arguments
if nargin < 2
  error('Not enough input arguments.');
end

% get the dimension from b
N = length(b(:));

% Default parameters
if nargin < 6 || isempty(x),       x = zeros(N,1);  end
if nargin < 7 || isempty(restart), restart = 100;    end
if nargin < 8 || isempty(tol),     tol = 1e-6;      end
if nargin < 9 || isempty(maxit),   maxit = N;       end

fprintf("||b|| = %g\n", norm(b(:)));

    % [fpath, lpath, fintf, lintf, epath] = pathreordering2(reshape(epath,[nep length(epath)/nep])', mesh1.t2f);
% [A, B1, B2, C1, C2, D1, D2, DL, DU] = pathcompute2(AE, epath, fpath, lpath, fintf, lintf, mesh1.f2t, mesh1.t2f, mesh1.elcon);
% uh2 = pathapply2(A, B1, B2, C1, C2, D1, D2, DL, DU, F, fpath, fintf, nep);

if isstruct(BE)  
  if isfield(BE,'D') && isfield(BE,'A')
    ncf = size(BE.D,1)/2;  
    b = pathapply(BE.A, BE.B, BE.C, BE.D, reshape(b,ncf,[]), BE.fpath, BE.fintf, BE.nelems);  
    b = b(:);
  elseif isfield(BE,'D2')
    ncf = size(BE.D2, 1);  
    b = pathapply2(BE.A, BE.B1, BE.B2, BE.C1, BE.C2, BE.D1, BE.D2, BE.DL, BE.DU, reshape(b,ncf,[]), BE.fpath, BE.fintf, BE.nep);  
    b = b(:);
  elseif isfield(BE,'row_ptr2')  
    [neb, nfeb] = size(BE.face);
    ncf = size(BE.D,1);
    b1 = faceextract(reshape(b,ncf,[]), BE.face);    
    [r1, r2, r3, r4] = vector_compute(BE.C1, BE.C2, BE.C3, b1, BE.idr1, BE.idr2, BE.idr3, neb);    
    for i = 1:neb
      u4 = block_ilu0_solve(BE.row_ptr2, BE.col_ind2, squeeze(BE.D(:,:,i,:)), reshape(r4(:,i,:), [], 1));
      r4(:,i,:) = reshape(u4, size(r4,1), 1, size(r4,3));
    end   
    b1 = vector_apply(BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, r1, r2, r3, r4, BE.face, BE.idr1, BE.idr2, BE.idr3);
    b = faceinsert(reshape(0*b,ncf,[]), b1, BE.face);
    b = b(:);
  elseif isfield(BE,'row_ptr')
    [neb, nfeb] = size(BE.face);
    ncf = size(BE.A,1);   
    b1 = faceextract(reshape(b,ncf,[]), BE.face); % ncf * nb * nfb
%     for i = 1:neb
%       tm = block_ilu0_solve(BE.row_ptr, BE.col_ind, BE.A(:,:,:,i), reshape(b1(:,i,:), [ncf*nfeb 1]));
%       b1(:,i,:) = reshape(tm, [ncf 1 nfeb]);
%     end
    %b1 = crs_parblockilu0_solve(BE.row_ptr, BE.col_ind, BE.A, b1);    
    b1 = crs_parblockilu0_solve2(BE.Lind_ji, BE.Uind_ji, BE.Lnum_ji, BE.Unum_ji, BE.A, b1);    
    b = faceinsert(reshape(0*b,ncf,[]), b1, BE.face);
    b = b(:);    
  elseif isfield(BE,'face')
    [neb, nfeb] = size(BE.face);
    ncf = size(BE.A,1)/nfeb;
    b1 = faceextract(reshape(b,ncf,[]), BE.face);
    for i = 1:neb
      tm = BE.A(:,:,i)\reshape(b1(:,i,:), [ncf*nfeb 1]);
      b1(:,i,:) = reshape(tm, [ncf 1 nfeb]);
    end
    b = faceinsert(reshape(0*b,ncf,[]), b1, BE.face);
    b = b(:);
  end  
elseif size(BE,1) == size(AE,1) && size(BE,3) == size(AE,3)
  % elemental additive Schwarz preconditioner
  b = hdgmatvec(BE, b, f2e, elcon);
else
  % block jacobi preconditioner
  b = applyblockjacobi(BE, full(b));
end

fprintf("||P*r|| = %g\n", norm(b(:)));

% initialization
nrmb   = norm(b); 
flags  = inf; 
iter   = 0; 
cycle  = 0;

% allocate memory
hh = zeros(restart+1,restart);
v  = zeros(N,restart+1);
e1 = zeros(restart+1,1); e1(1)=1;
rev = zeros(restart,1);

while (1) 
    % perform matrix-vector multiplication      
    d = hdgmatvec(AE, x, f2e, elcon);
    
    if isstruct(BE)        
      if isfield(BE,'D') && isfield(BE,'A')
        d = pathapply(BE.A, BE.B, BE.C, BE.D, reshape(d,ncf,[]), BE.fpath, BE.fintf, BE.nelems);        
      elseif isfield(BE,'D2')
        d = pathapply2(BE.A, BE.B1, BE.B2, BE.C1, BE.C2, BE.D1, BE.D2, BE.DL, BE.DU, reshape(d,ncf,[]), BE.fpath, BE.fintf, BE.nep);  
      elseif isfield(BE,'row_ptr2')  
        b1 = faceextract(reshape(d,ncf,[]), BE.face);    
        [r1, r2, r3, r4] = vector_compute(BE.C1, BE.C2, BE.C3, b1, BE.idr1, BE.idr2, BE.idr3, neb);    
        for i = 1:neb
          u4 = block_ilu0_solve(BE.row_ptr2, BE.col_ind2, squeeze(BE.D(:,:,i,:)), reshape(r4(:,i,:), [], 1));
          r4(:,i,:) = reshape(u4, size(r4,1), 1, size(r4,3));
        end   
        b1 = vector_apply(BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, r1, r2, r3, r4, BE.face, BE.idr1, BE.idr2, BE.idr3);
        d = faceinsert(reshape(0*d,ncf,[]), b1, BE.face);                    
      elseif isfield(BE,'row_ptr')
        b1 = faceextract(reshape(d,ncf,[]), BE.face);
%         for i = 1:neb
%           tm = block_ilu0_solve(BE.row_ptr, BE.col_ind, BE.A(:,:,:,i), reshape(b1(:,i,:), [ncf*nfeb 1]));
%           b1(:,i,:) = reshape(tm, [ncf 1 nfeb]);
%         end
        %b1 = crs_parblockilu0_solve(BE.row_ptr, BE.col_ind, BE.A, b1); 
        b1 = crs_parblockilu0_solve2(BE.Lind_ji, BE.Uind_ji, BE.Lnum_ji, BE.Unum_ji, BE.A, b1);    
        d = faceinsert(reshape(0*d,ncf,[]), b1, BE.face);        
      elseif isfield(BE,'face')
        b1 = faceextract(reshape(d,ncf,[]), BE.face);
        for i = 1:neb
          tm = BE.A(:,:,i)\reshape(b1(:,i,:), [ncf*nfeb 1]);
          b1(:,i,:) = reshape(tm, [ncf 1 nfeb]);
        end
        d = faceinsert(reshape(0*d,ncf,[]), b1, BE.face);        
      end
    elseif size(BE,1) == size(AE,1) && size(BE,3) == size(AE,3)
      d = hdgmatvec(BE, d, f2e, elcon);
    else
      d = applyblockjacobi(BE, d);     
    end
    
    % compute the residual vector      
    r = b(:) - d(:);        
    beta = norm(r);
    v(:,1) = r/beta;

    %fprintf("||v|| = %g\n", norm(v(:,1)));
    
    res  = beta;
    iter = iter+1;
    rev(iter) = res;    
    for j = 1:restart

        % set flag=0 if convergence
        if res/nrmb <= tol
            flags = 0;
            break;
        end        
        
        % set flag=1 if reaching maximum iteration 
        if iter>=maxit
            flags = 1;
            break;
        end                
                 
        % perform matrix-vector multiplication        
        d = hdgmatvec(AE, v(:,j), f2e, elcon);
%         fprintf("||v|| = %g\n", norm(d));
        
        if isstruct(BE)  
          if isfield(BE,'D') && isfield(BE,'A')
            d = pathapply(BE.A, BE.B, BE.C, BE.D, reshape(d,ncf,[]), BE.fpath, BE.fintf, BE.nelems);            
          elseif isfield(BE,'D2')
            d = pathapply2(BE.A, BE.B1, BE.B2, BE.C1, BE.C2, BE.D1, BE.D2, BE.DL, BE.DU, reshape(d,ncf,[]), BE.fpath, BE.fintf, BE.nep);                
          elseif isfield(BE,'row_ptr2')  
            b1 = faceextract(reshape(d,ncf,[]), BE.face);    
            [r1, r2, r3, r4] = vector_compute(BE.C1, BE.C2, BE.C3, b1, BE.idr1, BE.idr2, BE.idr3, neb);    
            for i = 1:neb
              u4 = block_ilu0_solve(BE.row_ptr2, BE.col_ind2, squeeze(BE.D(:,:,i,:)), reshape(r4(:,i,:), [], 1));
              r4(:,i,:) = reshape(u4, size(r4,1), 1, size(r4,3));
            end   
            b1 = vector_apply(BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, r1, r2, r3, r4, BE.face, BE.idr1, BE.idr2, BE.idr3);
            d = faceinsert(reshape(0*d,ncf,[]), b1, BE.face);            
          elseif isfield(BE,'row_ptr')
            b1 = faceextract(reshape(d,ncf,[]), BE.face);
%             for i = 1:neb
%               tm = block_ilu0_solve(BE.row_ptr, BE.col_ind, BE.A(:,:,:,i), reshape(b1(:,i,:), [ncf*nfeb 1]));
%               b1(:,i,:) = reshape(tm, [ncf 1 nfeb]);
%             end
            %b1 = crs_parblockilu0_solve(BE.row_ptr, BE.col_ind, BE.A, b1);    
            b1 = crs_parblockilu0_solve2(BE.Lind_ji, BE.Uind_ji, BE.Lnum_ji, BE.Unum_ji, BE.A, b1);    
            d = faceinsert(reshape(0*d,ncf,[]), b1, BE.face);                    
          elseif isfield(BE,'face')
            b1 = faceextract(reshape(d,ncf,[]), BE.face);
            for i = 1:neb
              tm = BE.A(:,:,i)\reshape(b1(:,i,:), [ncf*nfeb 1]);
              b1(:,i,:) = reshape(tm, [ncf 1 nfeb]);
            end
            d = faceinsert(reshape(0*d,ncf,[]), b1, BE.face);                    
          end
        elseif size(BE,1) == size(AE,1) && size(BE,3) == size(AE,3)
          d = hdgmatvec(BE, d, f2e, elcon);
        else
          d = applyblockjacobi(BE, d);     
        end
%         fprintf("||v|| = %g\n", norm(d));
%         error("Matlab");
    
        % Arnoldi process (i.e., GS orthogonalization on the Krylov supspace)        
        v(:,j+1) = d(:);
        for i = 1:j        
            hh(i,j) = v(:,i)'*v(:,j+1);
            v(:,j+1) = v(:,j+1) - hh(i,j)*v(:,i);        
        end
%         hh(1:j,j) = v(:,1:j)'*v(:,j+1);
%         v(:,j+1) = v(:,j+1) - v(:,1:j)*hh(1:j,j);        
        hh(j+1,j) = norm(v(:,j+1));     
        if (hh(j+1,j) ~= 0.0)          
            v(:,j+1) = v(:,j+1)/hh(j+1,j);   
        else            
            break;
        end

        % solve the reduced system 
        y = hh(1:j+1,1:j)\(beta*e1(1:j+1));
        
        % compute the residual norm
        res = norm(beta*e1(1:j+1)-hh(1:j+1,1:j)*y);
        
        iter = iter + 1; 
        rev(iter) = res;
        
        if rem(iter,500)==0
          res
        end
    end 
          
    % compute the solution    
    x(:) = x(:) + v(:,1:length(y))*y;    

    cycle = cycle + 1;

    % stop if converging or reaching the maximum iteration
    if flags < inf 
        fprintf('gmres(%d) converges at %d iterations to a solution with relative residual %g\n', [restart iter res/nrmb]);       
        break; 
    end     
end
rev = rev/nrmb;
