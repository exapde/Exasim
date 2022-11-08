close all

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

polarCoordinates = 1;
plotEUV = 0;

% generate input files and store them in datain folder
load('inputFiles_n1');
pde.soltime = [60:60:1440,1440:25:(1440+2600)]; % steps at which solution are collected
pde.timestepOffset = 0;

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd,'dataout');
logrho = sol(:,1,:,:);
v0 = pde.physicsparam(18)/pde.physicsparam(21);
vx = v0*sol(:,2,:,:)./sqrt(exp(logrho));
vy = v0*sol(:,3,:,:)./sqrt(exp(logrho));
vz = v0*sol(:,4,:,:)./sqrt(exp(logrho));
T = pde.physicsparam(19)*sol(:,5,:,:)./sqrt(exp(logrho));

nt = length(pde.soltime);
nPoints = size(mesh.dgnodes,1);
nDim = size(mesh.dgnodes,2);
nElem = size(mesh.dgnodes,3);

if plotEUV
    sol2 = zeros(size(sol,1),6,size(sol,3),size(sol,4));
    t = 0;
    for iT = 1:nt
        indexT = pde.soltime(iT);
        t = t + pde.saveSolFreq*pde.dt(indexT);
        fprintf(('Computing EUV for time %d of %d\n'),iT,nt);
        for iElem = 1:nElem
            for iPoint = 1:nPoints
                xg = mesh.dgnodes(iPoint,:,iElem);
                u0 = sol(iPoint,:,iElem,iT);
                sol2(iPoint,6,iElem,iT) = EUVsource3D(u0,xg,t,pde.physicsparam,pde.externalparam);
            end
        end
    end
else
    sol2 = zeros(size(sol,1),5,size(sol,3),size(sol,4));
end

if polarCoordinates
    dgnodesPolar = zeros(size(mesh.dgnodes));
   [a,e,r] = cart2sph(mesh.p(1,:),mesh.p(2,:),mesh.p(3,:));
   p2 = [a' e' r']';
   
   ne = pde.ne;
   for ie = 1:ne
       dge = mesh.dgnodes(:,:,ie);
       [adg,edg,rdg] = cart2sph(dge(:,1),dge(:,2),dge(:,3));
       adg = mod(adg+pi,2*pi)-pi;
       dgnodesPolar(:,:,ie) = [adg edg rdg];
       
       te = mesh.t(:,ie);
       pe = p2(:,te);
       pe(1,:) = mod(pe(1,:)+pi,2*pi)-pi;

       sa = unique((pe(1,:)>0));
       if (length(sa)>1 && max(abs(pe(1,:))>2))
           p2(1,te) = pe(1,:) + (pe(1,:)<0)*2*pi;
           dgnodesPolar(:,1,ie) = adg + (adg<0)*2*pi;
       end
   end
   
   sina = sin(dgnodesPolar(:,1,:));
   cosa = cos(dgnodesPolar(:,1,:));
   sine = sin(dgnodesPolar(:,2,:));
   cose = cos(dgnodesPolar(:,2,:));
   
   er = [cosa.*sine,sina.*sine,cose];
   ea = [-sina,cosa,zeros(size(dgnodesPolar(:,1,:)))];
   ee = [cosa.*cose,sina.*cose,-sine];
   
   vr = vx.*er(:,1,:) + vy.*er(:,2,:) + vz.*er(:,3,:);
   va = vx.*ea(:,1,:) + vy.*ea(:,2,:) + vz.*ea(:,3,:);
   ve = vx.*ee(:,1,:) + vy.*ee(:,2,:) + vz.*ee(:,3,:);
   
    sol2(:,1,:,:) = log10(pde.physicsparam(20)*exp(logrho));
    sol2(:,2,:,:) = va;
    sol2(:,3,:,:) = ve;
    sol2(:,4,:,:) = vr;
    sol2(:,5,:,:) = T;
   
   dgnodesPolar(:,1,:) = dgnodesPolar(:,1,:)*180/pi;
   dgnodesPolar(:,2,:) = dgnodesPolar(:,2,:)*180/pi;
   dgnodesPolar(:,3,:) = (dgnodesPolar(:,3,:)-pde.physicsparam(16))*pde.physicsparam(18)/1000;

   p2(1,:) = p2(1,:)*180/pi;
   p2(2,:) = p2(2,:)*180/pi;
   p2(3,:) = (p2(3,:)-pde.physicsparam(16))*pde.physicsparam(18)/1000;
   
   mesh.dgnodes = dgnodesPolar;
   mesh.p = p2;
   
else
    
    sol2(:,1,:,:) = logrho;
    sol2(:,2,:,:) = vx;
    sol2(:,3,:,:) = vy;
    sol2(:,4,:,:) = vz;
    sol2(:,5,:,:) = T;
end

% nComps = size(sol2,2);
% xDG = reshape(mesh.dgnodes(:,1,:),nPoints*nElem,1);
% yDG = reshape(mesh.dgnodes(:,2,:),nPoints*nElem,1);
% zDG = reshape(mesh.dgnodes(:,3,:),nPoints*nElem,1);
% uDG = reshape(sol2,nPoints*nElem,nComps);
% 
% f1 = scatteredInterpolant(xDG,yDG,zDG,reshape(sol2(:,1,:,:),nPoints*nElem,1));
% f2 = scatteredInterpolant(xDG,yDG,zDG,reshape(sol2(:,2,:,:),nPoints*nElem,1));
% f3 = scatteredInterpolant(xDG,yDG,zDG,reshape(sol2(:,3,:,:),nPoints*nElem,1));
% f4 = scatteredInterpolant(xDG,yDG,zDG,reshape(sol2(:,4,:,:),nPoints*nElem,1));
% f5 = scatteredInterpolant(xDG,yDG,zDG,reshape(sol2(:,5,:,:),nPoints*nElem,1));
% 
% sol3 = sol2;
% sol3(:,1,:,:) = reshape(f1([xDG,yDG,zDG]),nPoints,1,nElem);
% sol3(:,2,:,:) = reshape(f2([xDG,yDG,zDG]),nPoints,1,nElem);
% sol3(:,3,:,:) = reshape(f3([xDG,yDG,zDG]),nPoints,1,nElem);
% sol3(:,4,:,:) = reshape(f4([xDG,yDG,zDG]),nPoints,1,nElem);
% sol3(:,5,:,:) = reshape(f5([xDG,yDG,zDG]),nPoints,1,nElem);




% % visualize the numerical solution of the PDE model using Paraview
if plotEUV
    pde.visscalars = {"density", 1, "temperature", 5, "EUV", 6};  % list of scalar fields for visualization
else
    pde.visscalars = {"density", 1, "temperature", 5};  % list of scalar fields for visualization
end
pde.visvectors = {"velocity", [2, 3, 4]}; % list of vector fields for visualization
xdg = vis(sol2,pde,mesh); % visualize the numerical solution
% xdg = vis(sol3,pde,mesh); % visualize the numerical solution
disp("Done!");
