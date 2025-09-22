function [p,t] = lesmeshwithwake2d(xf, yf, dlay, dwall, nx, ny, nw, xref, yref)

% (xf, yf): airfoil coordinates
% dlay: BL thickness
% dwall: thickness of the first element at the wall
% nx: number of elements in x-direction
% ny: number of elements in y-direction
% xref: refinement along the x-direction
% yref: refinement along the y-direction
  
% Example: 
% dlay=0.1; dwall=2e-5; nx=96; ny = 25;
% [p,t] = lesmesh2d(xf, yf, dlay, dwall, nx, ny, [1.3 2; 1.55 1.78], [0.03 0.01 0.003]);

if size(xref,2)~=2
    error('xref must have dimension Nx times 2');
end

% calculate the mesh ratio
c = 1 - dlay/dwall;
rat = fsolve(@(x) scalingfun(x,ny,c),[1;3]);
rat = rat(1); 

% scaling distribution over the normal direction
yv = zeros(ny+1,1);
yv(2) = dwall;
for i = 1:(ny-1)
    yv(i+2) = yv(i+1) + dwall*(rat^i);
end

if abs(yv(end)-dlay)>1e-8
    error("Something wrong with the input parameters (dlay, dwall, ny)");
end

% Uniform distribution over foil
ns = length(xf);
t = 0:ns-1;
spx = spline(t,xf);
spy = spline(t,yf);
ttp = distribute(nx,spx,spy,ns);
xv = zeros(nx+1,1);
xv(1:nx+1,1) = 2*ttp/(ns-1);

% refine according to xref
for i = 1:size(xref,1)
    ind1 = (xv < xref(i,1));
    ind2 = (xv >= xref(i,1)) & (xv <= xref(i,2));
    ind3 = (xv > xref(i,2));
    x1 = xv(ind1); % no refinement
    x2 = xv(ind2); % yes refinement
    x3 = xv(ind3); % no refinement
    % refine x2
    xw = 0.5*(x2(1:end-1)+x2(2:end));
    x2 = sort([x2; xw]);
    % merge grid points 
    xv = unique([x1; x2; x3]);
end

% trailing edge
xv = [xv(1); 0.5*(xv(1)+xv(2)); xv(2);  0.5*(xv(2)+xv(3)); xv(3:end)];
%xv = [xv(1:end-1);  0.5*(xv(end-1)+xv(end)); xv(end)];
xv = [xv(1:end-2); 0.5*(xv(end-2)+xv(end-1)); xv(end-1);  0.5*(xv(end-1)+xv(end)); xv(end)];

if rem(length(xv),2)==0
    xv = [xv(1:end-1); 0.5*(xv(end-1)+xv(end)); xv(end)];
end
    
% make the grid for the foil
[p,t] = quadgrid(xv,yv);

% refine according to yref
n = length(yref);
if n>0
    yref = sort(yref,'descend');
    for i = 1:n
        [p,t] = refineaty(p,t,yref(i));
    end
    [p,t] = fixmesh(p,t);
end
 
[p,t] = removeelemement(p, t, ['y>' num2str(dlay/1.75)]);
yv = yv(yv<=dlay/1.75);

% wake region
dlay = max(yv);
ny = length(yv)-1;
ttwp = 0:nw;
d0 = 0.025;
rat = 1.05;
ttwp = d0*(rat.^ttwp -1)/(rat-1);
ttwp = ttwp/ttwp(end);
xwl = -fliplr(ttwp);
xwr = ttwp + 2;
X = xwr(:)*ones(1,ny+1);

wg = (1:-1/nw:0)';
wg1 = 1-wg; 
al = 0.01;
yvu = linspace(0,dlay,ny+1);
Y = (wg+al*wg1)*(yv(:)') + (1-al)*wg1*(yvu(:)');
pwr = [X(:) Y(:)];

X = xwl(:)*ones(1,ny+1);
wg = (0:1/nw:1)';
wg1 = 1-wg;
Y = (wg+al*wg1)*(yv(:)') + (1-al)*wg1*(yvu(:)');
pwl = [X(:) Y(:)];

m = length(xwr);
n = length(yv);
tw = [1 2 m+2 m+1];
tw = kron(tw,ones(n-1,1))+kron(ones(size(tw)),(0:n-2)'*m);
tw = kron(tw,ones(m-1,1))+kron(ones(size(tw)),(0:m-2)');

% make grid for foil and wake
[p,t] = connectmesh(pwl,tw,p,t,1e-8);
[p,t] = connectmesh(p,t,pwr,tw,1e-8);

ilb = (p(:,2) == 0 & p(:,1) <= 0);
ilt = (p(:,2) == 0 & p(:,1) >= 2);

% plot grid 
figure(1);clf;simpplot(p,t);axis on;

lf = (p(:,1) >=0 & p(:,1) <=2);
lb = p(:,1) < 0;
lt = p(:,1) > 2;

p(lf,:) = map2foil(p(lf,:),xf,yf);
p(lb,:) = map2bottomwake(p(lb,:),xf,yf);
p(lt,:) = map2topwake(p(lt,:),xf,yf);

figure(2);clf;simpplot(p,t);axis on;

pt = p(ilt,:);
[~,ind] = sort(pt,1);
pt = pt(ind(:,1),:);
pb = p(ilb,:);
[~,ind] = sort(pb,1);
pb = pb(ind(:,1),:);
nx = size(pt,1)-1;

nlw = 8;
pw = pb;
for i = 1:nlw
    pw = [pw; i*pt/nlw + (nlw-i)*pb/nlw];
end

tw = [1, 2, nx+3, nx+2];
tw = kron(tw,ones(nx,1)) + kron(ones(size(tw)),(0:nx-1)');
tw = kron(tw,ones(nlw,1)) + kron(ones(size(tw)),(0:nlw-1)'*(nx+1));

np = size(p,1);
p = [p; pw];
t = [t; tw+np];

[p,t] = fixmesh(p,t);

% n = 5;
% a = [loginc(linspace(0,0.5,n)',2); logdec(linspace(0.5,1.0,n)',3)];
% pw = points(pb, pt, unique(a));


figure(3);clf;simpplot(p,t);axis on;

return;

% % % wake region
% ttwp = 0:nw;
% d0 = 0.025;
% rat = 1.15;
% ttwp = d0*(rat.^ttwp -1)/(rat-1);
% ttwp = ttwp/ttwp(end);
% xwl = -fliplr(ttwp);
% xwr = ttwp + 2;
% [pwl,twl] = quadgrid(xwl,yv);
% [pwr,twr] = quadgrid(xwr,yv);
% 
% % make grid for foil and wake
% [p,t] = connectmesh(pwl,twl,p,t,1e-8);
% [p,t] = connectmesh(p,t,pwr,twr,1e-8);

% 
% % refine according to yref
% n = length(yref);
% if n>0
%     yref = sort(yref,'descend');
%     for i = 1:n
%         [p,t] = refineaty(p,t,yref(i));
%     end
%     [p,t] = fixmesh(p,t);
% end
% 
% lf = (p(:,1) >=0 & p(:,1) <=2);
% lb = p(:,1) < 0;
% lt = p(:,1) > 2;
% 
% x = p(lt,1);
% y = p(lt,2);
% xa = unique(x);
% ya = unique(y);
% yu = linspace(0,dlay,length(ya))';
% al = 0.3;
% for i = 1:length(xa)
%     xi = xa(i);
%     in = abs(x-xi)<1e-8;
%     wg1 = xi-2; %[0,1]
%     wg  = 1-wg1;%[1,0]
%     a = (wg+al*wg1);
%     %yu = linspace(0,dlay,sum(in))';
%     [~,im] = ismember(y(in),ya);
%     y(in) = a*y(in) + (1-a)*yu(im);
% end
% p(lt,2) = y;

% a=min(y(:));
% b=max(y(:));
% alpha  = 3*(x-2);
% p(lt,2) = a + (b-a)*(1-exp(-alpha.*(y-a)/(b-a)))./(1-exp(-alpha));

% p(lt,2) = logdec(x,alpha);
% wg1 = p(lt,1)-2; %[0,1]
% wg  = 1-wg1;
% al = 0.3;
% a = (wg+al*wg1);
% p(lt,2) = a.*p(lt,2) + (1-a).*pu;

% y = a(x)*y + (1-a(x))*yu
% x = 0 -> a(x) = 1   -> y = y
% x = 1 -> a(x) = 0.3 -> y = 0.3*y + 0.7*yu

% map rect to foil
% p = map2foil(p,xf,yf);
% 
% % plot grid 
% figure(2);clf;simpplot(p,t);axis on;


