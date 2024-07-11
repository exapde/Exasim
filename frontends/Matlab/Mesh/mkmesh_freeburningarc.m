function mesh = mkmesh_freeburningarc(porder,nd)

if nargin==0
    porder = 1; nd  = 2;
elseif nargin == 1
    nd = 2;
end

G = 1;
D = 26-G; % diameter
L = 10; % gap size
T = 3;  % tip 
H = 13; % height

% h = 0.5; % mesh size
% n1 = 70; n2 = 60; n3 = 54; n4 = 50; n5 = 30; n6 = 20;

% h = 1; % mesh size
% n1 = 44; n2 = 60; n3 = 31; n4 = 29; n5 = 15; n6 = 10;
% 
% a1 = logdec(linspace(0,-L,n1)',0.7);
% a2 = loginc(linspace(0,D,n2)',0.6); a2 = a2(2:end);
% a3 = linspace(-L,T+H,n3)'; a3 = a3(2:end);
% a4 = linspace(D,T*tan(pi/6),n4)'; a4 = a4(2:end);
% a5 = linspace(T+H,T,n5)'; a5 = a5(2:end);
% a6 = linspace(T*tan(pi/6),0,n6)'; a6 = a6(2:end-1);
% a7 = linspace(3,0,n6)'; a7 = a7(2:end-1);
% pv{1} = [0*a1 a1; a2 -L*ones(length(a2),1); D*ones(length(a3),1) a3;...
%          a4 (T+H)*ones(length(a4),1); T*tan(pi/6)*ones(length(a5),1) a5; a6 a7];
% [p,t]=polymesh(pv,[1,0],[1,0;0,1],[h,1.1]);
% 
% elemtype = 0;
% nodetype = 1;
% 
% bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-4)',...
%            'all(p(:,1)>max(p0(:,1))-1e-4)',...
%            'all(p(:,2)>max(p0(:,2))-1e-4)',...
%            'all(p(:,1)<min(p0(:,1))+1e-4)',...
%            'true'};     
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

elemtype = 1;
nodetype = 1;

n1 = 11; n2 = 11; n3 = 4; n4=7; n7 = 16;
a1 = 2.5; a2 = 0.5;

% n1 = 16; n2 = 16; n3 = 6; n4=12; n7 = 24;
% a1 = 2.5; a2 = 0.5;

% n1 = 30; n2 = 31; n3 = 12; n4=20; n7 = 40;
% a1 = 2.5; a2 = 0.5;

[p1,t1] = squaremesh(n1,n2,0,elemtype);
p1(:,1) = loginc(p1(:,1),a1);
p1(:,2) = logdec(p1(:,2),a2);
p1(:,1) = D*p1(:,1);
p1(:,2) = -L+L*p1(:,2);

[p2,t2] = squaremesh(n1,n3,0,elemtype);
p2(:,1) = loginc(p2(:,1),a1);
p2(:,2) = loginc(p2(:,2),a2);
p2(:,2) = T*p2(:,2);
p2(:,1) = p2(:,2)*tan(pi/6) + (D-p2(:,2)*tan(pi/6)).*p2(:,1);

[p3,t3] = squaremesh(n1,n4,1,elemtype);
p3(:,1) = loginc(p3(:,1),a1);
p3(:,2) = loginc(p3(:,2),a2+0.5);
p3(:,1) = T*tan(pi/6)+(D-T*tan(pi/6))*p3(:,1);
p3(:,2) = T+H*p3(:,2);

[p,t] = connectmesh(p1,t1,p2,t2);
[p,t] = connectmesh(p,t,p3,t3);
p(:,1) = p(:,1)+G;

if nd==2
    n5 = 5;
    [p4,t4] = squaremesh(n5,n2,0,elemtype);      
    p4(:,2) = logdec(p4(:,2),a2);
    p4(:,1) = G*p4(:,1);
    p4(:,2) = -L+L*p4(:,2);    
    [p,t] = connectmesh(p,t,p4,t4);
    
    %[p,t] = fixmesh(p,t);
%     pl = p; 
%     pl(:,1)=-pl(:,1);    
%     tl = t(:,[2 1 4 3]);    
%     [p,t] = connectmesh(p,t,pl,tl);
    %p = pl; t = tl;
end

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-4)',...
           'all(p(:,1)>max(p0(:,1))-1e-4)',...
           'all(p(:,2)>max(p0(:,2))-1e-4)',...
           'all(p(:,1)<min(p0(:,1))+1e-4)',...
           'true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

if nd==3    
    
    mesh2d = mesh;
    %mesh2d = mkmesh_square(11,11,porder,0,1,1,1,1);    
    zz = linspace(0,1,n7+1);
    z1 = min(zz);
    z2 = max(zz);    
    bndexpr = {['p(:,3)<' num2str(z1) '+1e-3'],...
               ['p(:,3)>' num2str(z2) '-1e-3'],...
               'true'};                     
    mesh1 = mkmesh_extrude2dmesh(mesh2d,zz,bndexpr);
    mesh1.p = mesh1.p(:,[3 1 2]);
    mesh1.dgnodes = mesh1.dgnodes(:,[3 1 2],:);
    
    alfa = 2*pi;
    R = mesh1.p(:,2);
    t = alfa*mesh1.p(:,1);
    mesh1.p(:,1) = R.*sin(t);
    mesh1.p(:,2) = R.*cos(t);
    mesh1.p(:,3) = mesh1.p(:,3);
        
    R = mesh1.dgnodes(:,2,:);
    t = alfa*mesh1.dgnodes(:,1,:);
    mesh1.dgnodes(:,1,:) = R.*sin(t);
    mesh1.dgnodes(:,2,:) = R.*cos(t);
    mesh1.dgnodes(:,3,:) = mesh1.dgnodes(:,3,:);    
    
    % remove duplicate nodes
    %dgnodes = mesh1.dgnodes;
    snap=1e-6;
    p = round(mesh1.p/snap)*snap;      
    [~,ix,jx]=unique(p,'rows');
    p=p(ix,:);
    t=jx(mesh1.t); 
%     for i = 1:mesh1.ne
%         e = p(t(i,:),:)-mesh1.p(mesh1.t(i,:),:);
%         d = sqrt(e(:,1).^2+e(:,2).^2+e(:,3).^2);
%         max(d)
%     end
%     bndexpr = {'all(p(:,3)<min(p0(:,3))+1e-4)',...
%                'all(p(:,3)>max(p0(:,3))-1e-4)',...               
%                'true'};         
%      mesh1 = mkmesh(p,t,porder,bndexpr,elemtype,nodetype); 
    mesh1.p = p;
    mesh1.t = t;
    %mesh1.dgnodes = dgnodes;
     %mesh1.p = round(mesh1.p/snap)*snap;   
    
    mesh2 = mkmesh_circle(porder,elemtype,nodetype,G,0.6*G,n7/4+1,2);    
    zz = unique(p1(:,2));
    z1 = min(zz);
    z2 = max(zz);    
    bndexpr = {['p(:,3)<' num2str(z1) '+1e-3'],...
               ['p(:,3)>' num2str(z2) '-1e-3'],...
               'true'};                     
    mesh2 = mkmesh_extrude2dmesh(mesh2,zz,bndexpr);
               
    mesh2.p = round(mesh2.p/snap)*snap;
    [p,t] = connectmesh(mesh2.p,mesh2.t,mesh1.p,mesh1.t);
    %[p,t] = connectmesh(mesh1.p,mesh1.t,mesh2.p,mesh2.t);
    bndexpr = {'all(p(:,3)<min(p0(:,3))+1e-4)',...
               'all(p(:,3)>max(p0(:,3))-1e-4)',...
               'all(sqrt(p(:,1).^2+p(:,2).^2)>10)',...
               'true'};         
    mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype); 
    mesh.dgnodes = cat(3,mesh2.dgnodes,mesh1.dgnodes);    
    %mesh.dgnodes = cat(3,mesh1.dgnodes,mesh2.dgnodes);    
    
%     mesh2 = mkmesh_circle(porder,elemtype,nodetype,G,0.6*G,n7/4+1,2);    
%     zz
%     zz = linspace(-10,0,11)    
%     bndexpr = {['p(:,3)<' num2str(z1) '+1e-3'],...
%                ['p(:,3)>' num2str(z2) '-1e-3'],...
%                'true'};                     
%     mesh2 = mkmesh_extrude2dmesh(mesh2,zz,bndexpr);
    
%     figure(1);clf; meshplot(mesh1,0); %simpplot(mesh1.p,mesh1.t);
%     figure(2);clf; meshplot(mesh2,0); %simpplot(mesh2.p,mesh2.t);
end

%close all; simpplot(p,t); hold on; plot(pv{1}(:,1),pv{1}(:,2),'or'); axis on;


% h = 0.5;
% a1 = logdec(linspace(0,-10,70)',1);
% a2 = loginc(linspace(0,15,40)',1); a2 = a2(2:end);
% a3 = linspace(-10,3,35)'; a3 = a3(2:end);
% a4 = loginc(linspace(15,3*tan(pi/6),35)',1); a4 = a4(2:end);
% a5 = linspace(3*tan(pi/6),0,25)'; a5 = a5(2:end-1);
% a6 = linspace(3,0,25)'; a6 = a6(2:end-1);
% pv{1} = [0*a1 a1; a2 -10*ones(length(a2),1); 15*ones(length(a3),1) a3;...
%          a4 3*ones(length(a4),1); a5 a6];
% [p,t]=polymesh(pv,[1,0],[1,0;0,1],[h,1.1]);
% 
% elemtype = 0;
% nodetype = 1;
% 
% bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-4)',...
%            'all(p(:,1)>max(p0(:,1))-1e-4)',...
%            'all(p(:,2)>max(p0(:,2))-1e-4)',...
%            'all(p(:,1)<min(p0(:,1))+1e-4)',...
%            'true'};     
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
% 
% close all; meshplot(mesh);

%pv{1} = [0 0; 0 -10; 15 -10; 15 3; 3*tan(pi/6) 3];

