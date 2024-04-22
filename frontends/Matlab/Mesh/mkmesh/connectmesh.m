function [p,t] = connectmesh(p1,t1,p2,t2,tol)
% connect (p1,t1) and (p2,t2) to form (p,t)

if nargin<5
    tol = 1e-6;
end

nd1 = size(p1,2);
np1 = size(p1,1);
np2 = size(p2,1);

% find duplicated nodes in p2 
in1 = [];
in2 = [];
for i=1:np2    
    e = (p1(:,1)-p2(i,1)).^2;
    for k=2:nd1 
        e = e + (p1(:,k)-p2(i,k)).^2;
    end
    [tm,k] = min(e);
    if tm<tol^2   
        in2 = [in2 i];        
        in1 = [in1 k];
    end
end

%[p1(in1,:) p2(in2,:)]

% remove duplicated nodes
p2(in2,:) = [];

% fix t2
for i = 1:size(t2,1)
    for j = 1:size(t2,2)
        m = t2(i,j);
        e = m-in2;
        [tm,k] = min(abs(e));
        if tm==0            
            t2(i,j) = in1(k);
        else
            n = nnz(e > 0);
            t2(i,j) = t2(i,j)+np1-n;
        end        
    end
end
    
% form (p,t)
p   = [p1;p2];
t   = [t1;t2];
