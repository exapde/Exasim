function sol = polarTransform(udg,pde,dgnodesPolar)

sol = zeros(size(udg,1),5,size(udg,3),size(udg,4));


dgnodesPolar(:,1,:) = dgnodesPolar(:,1,:)*pi/180;
dgnodesPolar(:,2,:) = dgnodesPolar(:,2,:)*pi/180;
dgnodesPolar(:,3,:) = dgnodesPolar(:,3,:)*1000/pde.physicsparam(18) + pde.physicsparam(16);

logrho = udg(:,1,:,:);
v0 = pde.physicsparam(18)/pde.physicsparam(21);
vx = v0*udg(:,2,:,:)./sqrt(exp(logrho));
vy = v0*udg(:,3,:,:)./sqrt(exp(logrho));
vz = v0*udg(:,4,:,:)./sqrt(exp(logrho));
T = pde.physicsparam(19)*udg(:,5,:,:)./sqrt(exp(logrho));
   
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

sol(:,1,:,:) = pde.physicsparam(20)*exp(logrho);
sol(:,2,:,:) = va;
sol(:,3,:,:) = ve;
sol(:,4,:,:) = vr;
sol(:,5,:,:) = T;