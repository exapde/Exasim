function [p1,t1,c1,p2,t2,c2,p3,c3,t3] = marchingtet(xdg, ydg, zdg, udg, isoval, cdg)

nd = 3;
np = size(xdg,1);
if np~=4
    error('The first dimension of the input must be 4.');
end

umin = min(udg,[],1);
umax = max(udg,[],1);

% find the tetrahedral that contain the isosurface
inde = find((umin<isoval) & (isoval < umax));

% get the number of elements containing the isosurface
ne   = length(inde);

% get the scalar field on these tetrahedral
udge = udg(:,inde);
cdge = cdg(:,inde);

% coordinates of the tetrahedral containing the isosurface
pdgex = reshape(xdg(:,inde),[np ne]);
pdgey = reshape(ydg(:,inde),[np ne]);
pdgez = reshape(zdg(:,inde),[np ne]);

% sort the scalar field for each tet 
[udgs,inds] = sort(udge,1);

% get the indices such that udgs = udge(indx)
indx = inds+repmat((0:(ne-1))*np,[np 1]);

% coordinates of the vertex-sorted tetrahedral containing the isosurface
px = pdgex(indx);
py = pdgey(indx);
pz = pdgez(indx);

% the scalar field on the isosurface 
c = cdge(indx);

% clear the redundant variables
clear udge cedge pdgex pdgey pdgez;

%%%%%%%%%%%%%%%%%%%%%%%% Case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
ind1 = find(udgs(2,:) > isoval);
ne1 = length(ind1);

% the triangular mesh and scalar field on the isosurface for the first case  
p1 = zeros(np-1,ne1,nd);
t1 = reshape(1:(np-1)*ne1,[np-1,ne1])';
c1 = zeros(np-1,ne1);

p1x = px(1,ind1); p1y = py(1,ind1); p1z = pz(1,ind1); valp1 = udgs(1,ind1); col1 = c(1,ind1);
for j = 1:3 % for each vertex
    k = j + 1;
    p2x = px(k,ind1); p2y = py(k,ind1); p2z = pz(k,ind1); valp2 = udgs(k,ind1); col2 = c(k,ind1);
    p = InterpolateVertices(isoval,p1x',p1y',p1z',p2x',p2y',p2z',valp1',valp2',col1',col2');
    p1(j,:,1:nd) = p(:,1:nd);
    c1(j,:) = p(:,end);
end

p1 = reshape(p1,[(np-1)*ne1,nd]);

%%%%%%%%%%%%%%%%%%%%%%%% Case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
ind2 = find(udgs(3,:) < isoval);
ne2 = length(ind2);

if ne2 > 0
    % the triangular mesh and scalar field on the isosurface for the second case  
    p2 = zeros(np-1,ne2,nd);
    t2 = reshape(1:(np-1)*ne2,[np-1,ne2])';
    c2 = zeros(np-1,ne2);

    p1x = px(4,ind2); p1y = py(4,ind2); p1z = pz(4,ind2); valp1 = udgs(4,ind2); col1 = c(4,ind2);
    for j = 1:3 % for each vertex
        k = j;
        p2x = px(k,ind2); p2y = py(k,ind2); p2z = pz(k,ind2); valp2 = udgs(k,ind2); col2 = c(k,ind2);
        p = InterpolateVertices(isoval,p1x',p1y',p1z',p2x',p2y',p2z',valp1',valp2',col1',col2');
        p2(j,:,1:nd) = p(:,1:nd);
        c2(j,:) = p(:,end);
    end

    p2 = reshape(p2,[(np-1)*ne2,nd]);

else
    p2 = [];
    t2 = [];
    c2 = [];
end

%%%%%%%%%%%%%%%%%%%%%%%% Case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
ind3 = setdiff((1:ne),[ind1 ind2]);
ne3 = length(ind3);

if ne3 > 0
    % the quadrilateral mesh and scalar field on the isosurface for the third case  
    p3 = zeros(np,ne3,nd);
    t3 = reshape(1:np*ne3,[np,ne3])';
    c3 = zeros(np,ne3);

    for j = 1:4 % for each vertex
        if j==1
            i1 = 1; i2 = 3;
        elseif j==2
            i1 = 1; i2 = 4;
        elseif j==3
            i1 = 2; i2 = 4;
        elseif j==4
            i1 = 2; i2 = 3;    
        end
        p1x = px(i1,ind3); p1y = py(i1,ind3); p1z = pz(i1,ind3); valp1 = udgs(i1,ind3); col1 = c(i1,ind3);
        p2x = px(i2,ind3); p2y = py(i2,ind3); p2z = pz(i2,ind3); valp2 = udgs(i2,ind3); col2 = c(i2,ind3);
        p = InterpolateVertices(isoval,p1x',p1y',p1z',p2x',p2y',p2z',valp1',valp2',col1',col2');
        p3(j,:,1:nd) = p(:,1:nd);
        c3(j,:) = p(:,end);
    end
    p3 = reshape(p3,[np*ne3,nd]);
else
    p3 = [];
    t3 = [];
    c3 = [];
end

if nargout<1
    figure(1); clf;   
    hold on; 
    patch('vertices',p1,'faces',t1,'cdata',c1(:), ...
               'facecol','interp','edgec','none');    
    patch('vertices',p2,'faces',t2,'cdata',c2(:), ...
               'facecol','interp','edgec','none');    
    patch('vertices',p3,'faces',t3,'cdata',c3(:), ...
               'facecol','interp','edgec','none');                  
    hold off;       
    axis equal; axis tight;
    camlight('headlight'), lighting gouraud
    set(gcf,'rend','z');
    fancycamera;
end

%    set(gcf,'rend','z');
    %colorbar,axis equal,drawnow;
    %cameramenu;   

function p = InterpolateVertices(isolevel,p1x,p1y,p1z,p2x,p2y,p2z,valp1,valp2,col1,col2)

if nargin == 9
    p = zeros(length(p1x), 3);
elseif nargin == 11
    p = zeros(length(p1x), 4);
else
    error('Wrong number of arguments');
end
mu = zeros(length(p1x), 1);
id = abs(valp1-valp2) < (10*eps) .* (abs(valp1) + abs(valp2));
if ( any(id) )
    p(id, 1:3) = [ p1x(id), p1y(id), p1z(id) ];
    if ( nargin == 11 )
        p(id, 4) = col1(id);
    end
end
nid = ~id;
if any(nid)
    mu(nid) = (isolevel - valp1(nid)) ./ (valp2(nid) - valp1(nid));    
    p(nid, 1:3) = [p1x(nid) + mu(nid) .* (p2x(nid) - p1x(nid)), ...
        p1y(nid) + mu(nid) .* (p2y(nid) - p1y(nid)), ...
        p1z(nid) + mu(nid) .* (p2z(nid) - p1z(nid))];
    if nargin == 11
        p(nid, 4) = col1(nid) + mu(nid) .* (col2(nid) - col1(nid));
    end
end









