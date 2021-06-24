function [X,Y] = foilcart1(xfl,yfl,xfu,yfu,nl,nu,nr,porder,a,b,R)

[x,y] = foilx(xfl,yfl,xfu,yfu,nl,nu,porder,a);

[nx,ny] = foilnormal(x,y);

p = masternodes(porder,1,1,0);
r = zeros(porder+1,nr);
z = R*loginc(linspace(0,1,nr+1)',b); 
for i = 1:length(z)-1
    r(:,i) = z(i) + (z(i+1)-z(i))*p;
end
r = r(1:end-1,:);
r = [r(:); R];

%r = R*loginc(linspace(0,1,nr*porder+1)',b);


Na = length(x);
Nr = length(r);
X = zeros(Na,Nr);
Y = zeros(Na,Nr);

X(:,1) = x;
Y(:,1) = y;
for i=2:Nr
    %[nx,ny] = foilnormal(X(:,i-1),Y(:,i-1));
    X(:,i) = X(:,i-1) + (r(i)-r(i-1))*nx;
    Y(:,i) = Y(:,i-1) + (r(i)-r(i-1))*ny;
end

a = loginc(linspace(1, 0, round(nr*0.25)*porder+1),2);
for i=0:length(a)-1
    [X(:,Nr-i),Y(:,Nr-i)] = fixl(X(:,Nr-i),Y(:,Nr-i),a(i+1));
    [X(:,Nr-i),Y(:,Nr-i)] = fixu(X(:,Nr-i),Y(:,Nr-i),a(i+1));
end

% figure(1); clf; 
% plot(X(:,1),Y(:,1),'-b','LineWidth',1);
% hold on;
% for i=2:Nr
%     plot(X(:,i),Y(:,i),'-b','LineWidth',1);
% end
% axis equal; axis tight;

function [x,y] = fixl(x,y,a)

ind = (y<=0) & (x<0.75);
x1 = x(ind);
y1 = y(ind);
dy = abs(y1(1:end-1)-y1(2:end));
[~,ii] = min(dy);
x0 = x1(ii);

ind = find(x>x0 & y<0);
%x1 = x(ind(1));
y1 = y(ind(1));
%x2 = x(ind(end));
y2 = y(ind(end));
ymin = min([y1 y2]);

%c = [x1 1; x2 1]\[y1; y2];
y(ind) = a*ymin + (1-a)*y(ind);


function [x,y] = fixu(x,y,a)

ind = y>=0;
x1 = x(ind);
y1 = y(ind);
dy = abs(y1(1:end-1)-y1(2:end));
[~,ii] = min(dy);
x0 = x1(ii);

ind = find(x>x0 & y>0);
%x1 = x(ind(1));
y1 = y(ind(1));
%x2 = x(ind(end));
y2 = y(ind(end));

%c = [x1 1; x2 1]\[y1; y2];
ymax = max([y1 y2]); % c(1)*x(ind) + c(2);

y(ind) = a*ymax + (1-a)*y(ind);




