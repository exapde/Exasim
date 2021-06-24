function [nx,ny] = foilnormal(x,y)

n = length(x);
nx = 0*x;
ny = nx;

for i = 1:n
    if (abs(x(i))<=1e-10)
        nx(i) = -1; ny(i) = 0;
    elseif (i==1)
        [nx(i),ny(i)] = normal(x(i),y(i),x(i+1),y(i+1));
    elseif (i==n)
        [nx(i),ny(i)] = normal(x(i-1),y(i-1),x(i),y(i));
    else
        [nx1,ny1] = normal(x(i),y(i),x(i+1),y(i+1));
        [nx2,ny2] = normal(x(i-1),y(i-1),x(i),y(i));
        nx3 = (nx1+nx2)/2;
        ny3 = (ny1+ny2)/2;
        s  = sqrt(nx3^2+ny3^2);
        nx(i) = nx3/s;
        ny(i) = ny3/s;
    end
end

% figure(1); clf; 
% plot(x,y,'-b');
% hold on;
% plot(x+0.01*nx,y+0.01*ny,'-r');
% axis equal; axis tight;

function [nx,ny] = normal(x1,y1,x2,y2)

nx = (y2-y1);
ny = -(x2-x1);
s  = sqrt(nx^2+ny^2);
nx = nx/s;
ny = ny/s;

% figure(1); clf;
% plot([x1 x2],[y1 y2]);
% hold on;
% plot([x1 x1+0.01*nx],[y1 y1+0.01*ny],'-r');
% axis equal; axis tight;
% pause



