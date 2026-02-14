function [pn, tn] = sphericalfrustummesh(slant, nc)
% slant: coordinates of an inclined line of the frustum
% nc: number of elements along theta-direction

% the radii of the cross sections of the frustum
R = abs(slant(2,:));

% the radius of the left base
R1 = R(1);    

% the radius of the right base
R2 = R(end);  

% the heigh of the frustum
H = abs(slant(1,1) - slant(1,end)); 

% the height of the inner cone formed by the left base
L = H/(R2/R1-1); 

% a quarter of the unit circle mesh
[p,t] = quartercirclemesh(nc);
% [p,t] = circlemesh(nc);

ns = size(slant,2);
np = size(p,2);
points = zeros(3, np, ns);
for i = 1:ns
  r = R(i);        %  radius of the conical frustum cross section
  x = slant(1,i);  %  x-coordinate of the conical frustum cross section  
  y = r*p(2,:);    %  y-coordinate of points on the conical frustum cross section   
  z = r*p(1,:);    %  z-coordinate of points on the conical frustum cross section   
  a = sqrt(y.^2 + z.^2);   % radius of the original circle 
  d = abs(x - slant(1,1)); % distance from the top base of radius r
  h = L+d;                 % the height of the cone whose base is radius r
  l = sqrt(r^2 + h^2);     % slant height of the cone whose base is radius r
  m = sqrt(a.^2 + h^2);    % slant height of the inner cone whose base is radius a
  b = (a*l)./m;            % radius of the projected circle is computed by triangle similarity 
  s = (h*l)./m;            % slant height of the projected cone whose base is radius b
  points(2,:,i) = b.*p(2,:); % y-coordinate of points on the projected circle of radius b
  points(3,:,i) = b.*p(1,:); % z-coordinate of points on the projected circle of radius b
  points(1,:,i) = x - (s - h); % x-coordinate of points on the projected circle of radius b
end

pn = squeeze(points(:,:,1));
tn = [];
np = size(pn,2);
t0 = t;
for i = 2:ns
    pn = cat(2, pn, squeeze(points(:,:,i)));
    t1 = t0 + np;
    tloc = cat(1, t0, t1);
    tn = cat(2, tn, tloc);
    t0 = t1;
end

figure(1); clf; simpplot(pn',tn'); axis on; view(3)

% figure(3); clf; simpplot(pn',tn(1:8,100)'); 
% x=pn(:,tn(:,100))';
% hold on;
% simpplot(pn',tn(5:8,100)'); 
% for i = 1:size(x,1)  
%   plot3(x(i,1),x(i,2),x(i,3),'o','LineWidth',2);
%   text(x(i,1), x(i,2), x(i,3), ...
%         ['  ' num2str(i)], ...   % node number
%         'FontSize', 12, ...
%         'FontWeight', 'bold');
% end
% xlabel('x'); ylabel('y'); zlabel('z');
% axis on;

