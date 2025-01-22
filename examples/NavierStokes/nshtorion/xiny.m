function in = xiny(x,y,tol)
% Determine if each row of x is a member of y
% If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% Else in(j) = 0

if nargin < 3
    tol = 1e-12;
end

[m,dim] = size(x);
in = zeros(m,1);
if dim==1
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2;
        [md,id] = min(d2);
        if md<tol, in(j)=id; end
    end
elseif dim==2
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
        [md,id] = min(d2);
        if md<tol, in(j)=id; end
    end
elseif dim==3
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
        [md,id] = min(d2);
        if md<tol, in(j)=id; end
    end
else
    n = size(y,1);    
    for j=1:m
        d2 = sum((y - repmat(x(j,:),[n 1])).^2,2);
        [md,id] = min(d2);
        if md<tol, in(j)=id; end
    end    
end


