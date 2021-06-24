function [X3, Y3] = foilcart3(X,Y,nw,nt,porder,c,L)

x1 = X(1,:);
x2 = X(end,:);
y1 = Y(1,:);
y2 = Y(end,:);

if (abs(x1(1)-x2(1))<1e-10) && (abs(y1(1)-y2(1))<1e-10) 
    % closed trailing edge
    x = [x1 x2(2:end)];
    y = [y1 y2(2:end)];
else
    % open trailing edge    
    xt = linspace(x1(1),x2(1),porder*nt+1); 
    yt = linspace(y1(1),y2(1),porder*nt+1);   
    
    roundTE = 1;
    flattingRatio = 0.75; % 0 for flat TE, 1 for rounded TE
    if roundTE == 1
        TExcenter = (x1(1)+x2(1))/2;
        TEycenter = (y1(1)+y2(1))/2;
        TEradious = sqrt((x1(1)-x2(1))^2+(y1(1)-y2(1))^2) / 2;

        angle = acos( sqrt((xt-TExcenter).^2+(yt-TEycenter).^2) / TEradious );
        angle(yt<TEycenter) = pi-angle(yt<TEycenter);

        xt = TExcenter + flattingRatio*TEradious*sin(angle);
        yt = TEycenter + TEradious*cos(angle);
    end
    
    x = [x1 xt(2:end) x2(2:end)];
    y = [y1 yt(2:end) y2(2:end)];
end
[y,ind] = sort(y);
x = x(ind);
m = length(x);


p = masternodes(porder,1,1,0);
r = zeros(porder+1,nw);
z = loginc(linspace(0,L,nw+1)',c); 
for i = 1:length(z)-1
    r(:,i) = z(i) + (z(i+1)-z(i))*p;
end
r = r(1:end-1,:);
r = [r(:); L];
n = length(r);

X3 = zeros(m,n);
Y3 = zeros(m,n);
X3(:,1) = x(:);
Y3(:,1) = y(:);
yUniform = linspace(y(1),y(end),m);
nx = 1; ny = 0;
uniformWeight = 1.0;
yEnd = uniformWeight*yUniform(:) + (1-uniformWeight)*y(:);
iThreshold = ceil(1.0*n);
for i = 2:n
    x_tmp1 = X3(:,1) + r(i)*nx;
    x_tmp2 = X3(:,1) + r(n)*(i / n)*nx;
    y_tmp = Y3(:,1) + r(i)*ny; 
    blendingAlongY = (abs(y_tmp) / max(abs(y_tmp(:))));
    X3(:,i) = (1-blendingAlongY).*x_tmp1+blendingAlongY.*x_tmp2;% + 0.5*x_tmp2;
%     Y3(:,i) = Y3(:,i-1) + (r(i)-r(i-1))*ny;  
    if i>iThreshold
        Y3(:,i) = Y3(:,i-1);
    else
        Y3(:,i) = ((L+x(:)-X3(:,i)) .* Y3(:,1) + (X3(:,i) - x(:)) .* yEnd(:)) / L;
    end
end

% X3 = X3';
% Y3 = Y3';
