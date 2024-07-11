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

% block jacobi preconditioner
b = applyblockjacobi(BE, full(b));

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
    d = applyblockjacobi(BE, d);     
    
    % compute the residual vector      
    r = b(:) - d(:);        
    beta = norm(r);
    v(:,1) = r/beta;

    fprintf("||v|| = %g\n", norm(v(:,1)));
    
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
        
        d = applyblockjacobi(BE, d); 
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
