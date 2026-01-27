function [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata(master,mesh,app,UDG,UH,wid,elemAvg,deltaT,WDG)
w_flag = 1;
if nargin < 7; elemAvg = 0; end
if nargin < 8; deltaT = 1; end
if nargin < 9; WDG = 0*UDG; w_flag = 0; end

    % gam = app.arg{1};
    % Re = app.arg{3};
    % Pr = app.arg{4};
    % Minf = app.arg{2};
    % pinf = 1/(gam*Minf^2);
    % pinf = 1.69061708; % For T106C with Min = 0.65
    % pinf = 1.269841269841270; % For T106C with Min = 0.75
    % pinf = 476
    nd  = master.nd;
    ngf = master.ngf;
    npf = master.npf;
    ne = mesh.ne;
nfe = size(master.perm,2);

elcon  = mesh.elcon;
perm   = master.perm(:,:,1);
shapft = squeeze(master.shapft(:,:,1));
dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);

xc = [];
yc = [];
zc = [];
cpData = [];
cfData = [];
chData = [];

for i=1:ne
    dg = mesh.dgnodes(:,:,i);       
    uh = UH(:,elcon(:,i))';          % uh(3*nps,nc)
    udg = UDG(:,:,i);
    wdg = WDG(:,:,i);


    fc = mesh.f(:,i);
    
    for is = 1:nfe
        if fc(is) == wid                
                pn    = dg(perm(:,is),:);
                pg    = shapft*pn;
                udgg  = shapft*udg(perm(:,is),:);
                wdgg  = shapft*wdg(perm(:,is),:);
                uhg   = shapft*uh((is-1)*npf+1:is*npf,:);   
%                 uhg = udgg(:,1:size(UH,1));
                dpg   = dshapft*pn(:,1:nd);
                dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  

                if nd==2
                    jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
                    nlg   = [dpg(:,2),-dpg(:,1)];
                    nlg   = bsxfun(@rdivide, nlg, jac);
                elseif nd==3
                    nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
                    nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
                    nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
                    jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
                    nlg   = bsxfun(@rdivide, nlg, jac);
                end             
                
                [fh] = fhat_dim(nlg,pg,udgg,uhg,app.arg,0,wdgg,w_flag);                                 
                                
                pr = fh(:,6).*nlg(:,1) + fh(:,7).*nlg(:,2);
                if any(abs(nlg(:,2)) < 1.0e-2)
                    1;
                end
%                 pres = eulereval(udgg,'p',1.4);
%                 prew = eulereval(uhg,'p',1.4);
%                 [pr pres prew] 
%                 [udgg(:,2:3) uhg(:,2:3)]
%                 pause
                
                q = fh(:,8);
                if nlg(:,2) > 0     % This is only an approximation for the boundary layer side
                    tau = fh(:,2).*nlg(:,2) - fh(:,3).*nlg(:,1);
                else
                    tau = - fh(:,2).*nlg(:,2) + fh(:,3).*nlg(:,1);
                end
                
                if elemAvg == 0
                    xNew = pg(:,1);
                    yNew = pg(:,2);
                    if nd == 3
                        zNew = pg(:,3);
                    end
                    cpNew = pr;
                    cfNew = tau;
                    chNew = q;
                elseif elemAvg == 1
                    xNew = sum(abs(jac).*master.gwfc.*pg(:,1))/sum(abs(jac).*master.gwfc);
                    yNew = sum(abs(jac).*master.gwfc.*pg(:,2))/sum(abs(jac).*master.gwfc);
                    if nd == 3
                        zNew = sum(abs(jac).*master.gwfc.*pg(:,3))/sum(abs(jac).*master.gwfc);
                    end
                    cpNew = sum(abs(jac).*master.gwfc.*pr)/sum(abs(jac).*master.gwfc);
                    cfNew = sum(abs(jac).*master.gwfc.*tau)/sum(abs(jac).*master.gwfc);
                    chNew = sum(abs(jac).*master.gwfc.*q)/sum(abs(jac).*master.gwfc);
                else
                    error('elemAvg has invalid value');
                end
                
                xc = [xc; xNew];
                yc = [yc; yNew];
                if nd == 3
                    zc = [zc; zNew];
                end
                cpData = [cpData; cpNew];
                cfData = [cfData; cfNew];
                chData = [chData; chNew];
            end
        end
    end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%  AIRFOIL DEPENDENT %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % xcmax = max(xc);
% % % % % % ind = yc < -(1/xcmax^2)*0.1*xc.*xc+0.1;    % Trefftz foil dividing curve
% % 
% % NACA 65-(18)10:
% % alpha = (pi/180)*16.73;
% % theta = (pi/180)*45 - alpha;
% theta = 0;      % This for Eppler 387
% [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,0,0,1,18/10,1);
% ind = yWithNoRotation < yMeanLine;       % % NACA 65-(18)10
% % 
% ind = yc < 0;                             % $ Symemtric foil
%ind = (yc < 0.0015*xc.^2 & xc<0.99) | (xc>=0.99 & yc<=-0.1*xc+0.1);    
xWithNoRotation = xc;
yWithNoRotation = yc;
ind = yWithNoRotation < 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % T106C:
% theta = -0.527999463;
% [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,0,0,1);
% ind = yWithNoRotation < yMeanLine;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Cylinder
% xWithNoRotation = xc;
% yWithNoRotation = yc;
% ind = yWithNoRotation < 0;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xclower = xWithNoRotation(ind);
yclower = yWithNoRotation(ind);
if nd == 3; zclower = zc(ind); end
datalower = cpData(ind);
data1lower = cfData(ind);
data2lower = chData(ind);

xcupper = xWithNoRotation(~ind);
ycupper = yWithNoRotation(~ind);
if nd == 3; zcupper = zc(~ind); end
dataupper = cpData(~ind);
data1upper = cfData(~ind);
data2upper = chData(~ind);

[xcl,ij] = sort(xclower);
ycl = yclower(ij);
if nd == 3; zcl = zclower(ij); end
datal = datalower(ij);
data1l = data1lower(ij);
data2l = data2lower(ij);

[xcu,ij] = sort(xcupper);
ycu = ycupper(ij);
if nd == 3; zcu = zcupper(ij); end
datau = dataupper(ij);
data1u = data1upper(ij);
data2u = data2upper(ij);

% x = [flipud(xcl); xcu];
% xmin = min(x);
% xmax = max(x);
% x = (x-xmin)/(xmax-xmin);
% y = [flipud(ycl); ycu]/(xmax-xmin);
% z = [flipud(zcl); zcu]/(xmax-xmin);
x = [flipud(xcl); xcu];
y = [flipud(ycl); ycu];
if nd == 3; z = [flipud(zcl); zcu]; end 
if nd == 2
    x = [x(:) y(:)];
elseif nd == 3
    x = [x(:) y(:) z(:)];
end

Cp = -([flipud(datal); datau]);
Cf = -2*[flipud(data1l); data1u];
Ch = [flipud(data2l); data2u];

% if nd == 3
%     [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch);
% else
    Cp2d = [];
    Cf2d = [];
    x2d = [];
    Ch2d = [];
% end

% figure(1)
% plot(x,Cp,'-'); 
% grid on;
% axis equal;
% axis tight;
% axis on;
% 
% figure(2) 
% plot(x,Cf); grid on;
% %axis square;
% axis tight;
% axis on;

