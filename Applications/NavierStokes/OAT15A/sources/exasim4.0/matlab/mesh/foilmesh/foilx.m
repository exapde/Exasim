function [x,y] = foilx(xfl,yfl,xfu,yfu,nl,nu,porder,a)

[xfl,ind] = sort(xfl);
yfl = yfl(ind); 

[xfu,ind] = sort(xfu);
yfu = yfu(ind); 

xlmin = xfl(1);
xlmax = xfl(end);
xumin = xfu(1);
xumax = xfu(end);

if abs(xfl(1)-xfu(1))>1e-12
    error('Airfoil geometry is not valid');
end
if abs(yfl(1)-yfu(1))>1e-12
    error('Airfoil geometry is not valid');
end
    
% xl = loginc(linspace(xlmin,xlmax,nl*porder+1)',a);
% xu = loginc(linspace(xumin,xumax,nu*porder+1)',a);

xl = zeros(porder+1,nl);
xu = zeros(porder+1,nu);
p = masternodes(porder,1,1,0);

x = loginc(linspace(xlmin,xlmax,nl+1)',a); 
% z = loginc(linspace(x(1),x(2),6)',2);   
% z = z(2:end-1)
z = x(2)*[1/240 1/40 1/16 1/8 1/5 1/3.25 0.45 0.625 0.8]';
x = sort([x(1); z; x(2); ...
    x(2)+[0.2; 0.435; 0.70]*(x(3)-x(2));...
    x(3); x(3)+0.45*(x(4)-x(3)); x(4);...
    x(4)+0.50*(x(5)-x(4)); x(5:end)]);

for i = 1:length(x)-1
    xl(:,i) = x(i) + (x(i+1)-x(i))*p;
end
xl = xl(1:end-1,:);
xl = [xl(:); xlmax];

x = loginc(linspace(xumin,xumax,nu+1)',a);
% z = [x(2)/16; x(2)/4];
% x = [x(1); x(2:end)];
% z = x(2)*[1/240 1/40 1/16 1/8 1/5 1/3 0.5 0.73]';
% x = sort([x(1); z; x(2); ...
%     x(2)+[0.25; 0.5; 0.75]*(x(3)-x(2));...
%     x(3); x(3)+0.45*(x(4)-x(3)); x(4);...
%     x(4)+0.50*(x(5)-x(4)); x(5:end)]);
z = x(2)*[1/240 1/40 1/16 1/8 1/5 1/3.25 0.45 0.625 0.8]';
x = sort([x(1); z; x(2); ...
    x(2)+[0.2; 0.435; 0.70]*(x(3)-x(2));...
    x(3); x(3)+0.45*(x(4)-x(3)); x(4);...
    x(4)+0.50*(x(5)-x(4)); x(5:end)]);
% x = sort([x(1); z; x(2); ...
%     x(2)+[0.25; 0.5; 0.75]*(x(3)-x(2));...
%     x(3); x(3)+0.45*(x(4)-x(3)); x(4:end)]);

for i = 1:length(x)-1
    xu(:,i) = x(i) + (x(i+1)-x(i))*p;
end
xu = xu(1:end-1,:);
xu = [xu(:); xumax];

%xu = loginc(linspace(xumin,xumax,nu+1)',a);


yl = interp1(xfl,yfl,xl,'linear','extrap');
yu = interp1(xfu,yfu,xu,'linear','extrap');
    
xu = xu(end:-1:2);
yu = yu(end:-1:2);

x = [xu; xl];
y = [yu; yl];

% xl = xl(end:-1:2);
% yl = yl(end:-1:2);
% 
% x = [xl; xu];
% y = [yl; yu];

%figure(1); clf; plot(x,y,'-o'); axis equal; axis tight;
