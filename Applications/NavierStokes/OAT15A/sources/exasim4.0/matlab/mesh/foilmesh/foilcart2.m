function [X, Y] = foilcart2(X,Y,nr2,porder,b2,R1,R2)

p = masternodes(porder,1,1,0);
r = zeros(porder+1,nr2);
z = loginc(linspace(R1,R2,nr2+1)',b2); 
for i = 1:length(z)-1
    r(:,i) = z(i) + (z(i+1)-z(i))*p;
end
r = r(1:end-1,:);
r = [r(:); R2];

%[nx,ny] = foilnormal(X(:,end),Y(:,end));
n = size(X,2);
for i=n+1:n+length(r)-1
%     x_tmp = X(:,i-1) + (r(i+1-n)-r(i-n))*nx;
%     y_tmp = Y(:,i-1) + (r(i+1-n)-r(i-n))*ny;
%     X(:,i) = x_tmp;
%     Y(:,i) = y_tmp;
    
    alpha = linspace(0,pi,size(X,1));
    x_tmp2 = X(:,i-1) - (r(i+1-n)-r(i-n))*sin(alpha(:));
    y_tmp2 = Y(:,i-1) + (r(i+1-n)-r(i-n))*cos(alpha(:));
    
    X(:,i) = 1.0*x_tmp2;% + 0.5*x_tmp2;
    Y(:,i) = 1.0*y_tmp2;% + 0.5*y_tmp2;
    
end

