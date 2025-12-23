function [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata2_hdgcode(master,mesh,app,UDG,UH,wid,elemAvg,deltaT)

    if nargin < 7; elemAvg = 0; end
    if nargin < 8; deltaT = 1; end
    fhat   = str2func(app.fhat);
    
    
    gam = app.arg{1};
    % Re = app.arg{3};
    % Pr = app.arg{4};
    Minf = app.arg{4};
    pinf = 1/(gam*Minf^2);
    % pinf = 1.69061708; % For T106C with Min = 0.65
    % pinf = 1.269841269841270; % For T106C with Min = 0.75
    
    
    nd  = master.nd;
    ngf = master.ngf;
    npf = master.npf;
    % ne  = size(mesh.t,1);
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
        
        fc = mesh.f(:,i);
    
        for is = 1:nfe
            if fc(is) == wid
                pn    = dg(perm(:,is),:);
                pg    = shapft*pn;
                udgg  = shapft*udg(perm(:,is),:);
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
                
                [fh] = fhat2(nlg,pg,udgg,uhg,app.arg,0);                                 
                                
                pr = fh(:,2).*nlg(:,1) + fh(:,3).*nlg(:,2);
                if any(abs(nlg(:,2)) < 1.0e-2)
                    1;
                end
%                 pres = eulereval(udgg,'p',1.4);
%                 prew = eulereval(uhg,'p',1.4);
%                 [pr pres prew] 
%                 [udgg(:,2:3) uhg(:,2:3)]
%                 pause
                
                q = fh(:,nd+2);
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
    
    
    Cp = -2*([flipud(datal); datau]-pinf);
    Cf = -2*[flipud(data1l); data1u];
    Ch = [flipud(data2l); data2u] / (deltaT*gam);
    
    
    if nd == 3
        [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch);
    else
        Cp2d = [];
        Cf2d = [];
        x2d = [];
        Ch2d = [];
    end
    
    
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
    
    
    end
    
    
    function [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch)
    % z = x(:,3);
    % z = z(:);snap = 1.0e-8;
    snap = 1.0e-9;
    xyPoints = unique(round(squeeze(x(:,1:2))/snap)*snap,'rows');
    num_xyPoints = size(xyPoints,1);
    Cp2d = zeros(num_xyPoints,1);
    Cf2d = zeros(num_xyPoints,1);
    Ch2d = zeros(num_xyPoints,1);
    
    
    tol = 1.0e-5;
    for i=1:num_xyPoints
        indices = find(sqrt((x(:,1)-xyPoints(i,1)).^2 + (x(:,2)-xyPoints(i,2)).^2) < tol);
        Cp2d(i) = sum(Cp(indices))/length(indices);
        Cf2d(i) = sum(Cf(indices))/length(indices);
        Ch2d(i) = sum(Ch(indices))/length(indices);
    end
    % % NACA 65-(18)10
    % [yMeanLine,~,yWithNoRotation] = getNACA65_1810meanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1,18/10,1);
    % ind = yWithNoRotation < yMeanLine;
    % % T106C
    % [yMeanLine,~,yWithNoRotation] = getT106CmeanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1);
    % ind = yWithNoRotation < yMeanLine;
    % % Cylinder
    ind = xyPoints(:,2) < 0;
    
    
    xclower2d = xyPoints(ind,1);
    yclower2d = xyPoints(ind,2);
    Cp2dlower = Cp2d(ind);
    Cf2dlower = Cf2d(ind);
    Ch2dlower = Ch2d(ind);
    
    
    xcupper2d = xyPoints(~ind,1);
    ycupper2d = xyPoints(~ind,2);
    Cp2dupper = Cp2d(~ind);
    Cf2dupper = Cf2d(~ind);
    Ch2dupper = Ch2d(~ind);
    
    
    [xcl2d,ij] = sort(xclower2d);
    ycl2d = yclower2d(ij);
    Cp2dlower = Cp2dlower(ij);
    Cf2dlower = Cf2dlower(ij);
    Ch2dlower = Ch2dlower(ij);
    
    
    [xcu2d,ij] = sort(xcupper2d);
    ycu2d = ycupper2d(ij);
    Cp2dupper = Cp2dupper(ij);
    Cf2dupper = Cf2dupper(ij);
    Ch2dupper = Ch2dupper(ij);
    
    
    x2d = [flipud(xcl2d); xcu2d];
    y2d = [flipud(ycl2d); ycu2d];
    x2d = [x2d(:) y2d(:)];
    Cp2d = [flipud(Cp2dlower); Cp2dupper];
    Cf2d = [flipud(Cf2dlower); Cf2dupper];
    Ch2d = [flipud(Ch2dlower); Ch2dupper];
    
    
    end
    
    
    function [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,xTranslation,yTranslation,scaling,cl,a)
    
    
    if nargin < 5; a = 1; end
    if a == 1; a = 0.999999; end
    
    
    A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
    xyzWithNoRotation = [xc(:) yc(:)];
    xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
    xWithNoRotation = xyzWithNoRotation(:,1);
    yWithNoRotation = xyzWithNoRotation(:,2);
    xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
    yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
    g = -(1/(1-a)) * (a^2*(0.5*log(a)-0.25) + 0.25);
    h = (1/(1-a)) * (0.5*(1-a)^2*log(1-a) - 0.25*(1-a)^2) + g;
    yMeanLine = (cl/(2*pi*(a+1))) * ( (1/(1-a))* (0.5*(a-xWithNoRotation).^2.*log(abs(a-xWithNoRotation)) - 0.5*(1-xWithNoRotation).^2.*log(1-xWithNoRotation) + ...
            0.25*(1-xWithNoRotation).^2 - 0.25*(a-xWithNoRotation).^2) - xWithNoRotation.*log(xWithNoRotation) + g - h*xWithNoRotation);
    
    
    end
    
    
    function [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,xTranslation,yTranslation,scaling)
    
    
    A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
    xyzWithNoRotation = [xc(:) yc(:)];
    xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
    xWithNoRotation = xyzWithNoRotation(:,1);
    yWithNoRotation = xyzWithNoRotation(:,2);
    xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
    yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
    
    
    % xyMidChord = [-0.004850771437097   0.008104085834677
    %    0.026267527614243   0.040006331527312
    %    0.058813891977280   0.154717136626771
    %    0.119764563959287   0.220726452638423
    %    0.187900992975752   0.274415205427912
    %    0.272808035650415   0.299349380953031
    %    0.368800317722688   0.305277014484574
    %    0.473686474141603   0.295955379128744
    %    0.600733470684711   0.248637186068960
    %    0.747001884020969   0.168362315002401
    %    0.931784011091559   0.022052558159896
    %    0.989643954790913   0.008104085275412];
    
    
    % xyMidChord = [0.003   0.003
    %    0.011900899226019   0.022010534473781
    %    0.022099897892833   0.089780617294087
    %    0.043846090146471   0.137752100197530
    %    0.073184363539016   0.172706344247624
    %    0.103314843692981   0.206302285554901
    %    0.136214284225593   0.235150619721944
    %    0.171832517935802   0.259337362483664
    %    0.203472806334644   0.290344615105094
    %    0.250885642315712   0.294308555695582
    %    0.294730428985119   0.304390206210480
    %    0.346638258603599   0.300647115652073
    %    0.391564089454086   0.308875228242563
    %    0.445012961862948   0.302489894836824
    %    0.501235349883211   0.291349143974849
    %    0.565068672054260   0.267158829241619
    %    0.634686967796641   0.233049708141524
    %    0.707440700461415   0.193564626144587
    %    0.788313850923384   0.140158143825069
    %    0.867798377388409   0.089132569969907
    %    0.965128210899839   0.007509769875418
    %    0.982   -4.8e-3];
    
    
    xyMidChord = [-0.004850771437097   0.008104085834677
       0.011900899226019   0.022010534473781
       0.022099897892833   0.089780617294087
       0.043846090146471   0.137752100197530
       0.073184363539016   0.172706344247624
       0.103314843692981   0.206302285554901
       0.136214284225593   0.235150619721944
       0.171832517935802   0.259337362483664
       0.203472806334644   0.290344615105094
       0.250885642315712   0.294308555695582
       0.294730428985119   0.304390206210480
       0.346638258603599   0.300647115652073
       0.391564089454086   0.308875228242563
       0.445012961862948   0.302489894836824
       0.501235349883211   0.291349143974849
       0.565068672054260   0.267158829241619
       0.634686967796641   0.233049708141524
       0.707440700461415   0.193564626144587
       0.788313850923384   0.140158143825069
       0.867798377388409   0.089132569969907
       0.965128210899839   0.007509769875418
       0.989643954790913   0.008104085275412];
    
    
    yMeanLine = interp1(xyMidChord(:,1),xyMidChord(:,2),xWithNoRotation);
    
    
    % % Mid-chord (x,y) coordinates with rotation:
    %   -0.000107288164668   0.009444293515160
    %    0.042845768527915   0.021324391866838
    %    0.128751881913082   0.104016436967154
    %    0.214657995298249   0.130329055583308
    %    0.300564108683416   0.142378736834755
    %    0.386470222068583   0.121140594544775
    %    0.472376335453749   0.077899471296392
    %    0.558282448838916   0.017004985691412
    %    0.644188562224083  -0.087876330128531
    %    0.730094675609250  -0.230910044153280
    %    0.816000788994417  -0.450389243885509
    %    0.858953845687000  -0.491588368894000
    %     
    % % Mid-chord (x,y) coordinates with no rotation:
    %     -0.004850771437097   0.008104085834677
    %    0.026267527614243   0.040006331527312
    %    0.058813891977280   0.154717136626771
    %    0.119764563959287   0.220726452638423
    %    0.187900992975752   0.274415205427912
    %    0.272808035650415   0.299349380953031
    %    0.368800317722688   0.305277014484574
    %    0.473686474141603   0.295955379128744
    %    0.600733470684711   0.248637186068960
    %    0.747001884020969   0.168362315002401
    %    0.931784011091559   0.022052558159896
    %    0.989643954790913   0.008104085275412
        
    end
    
    
    % 
    % % setapplicationpath('FM/ns');
    % % wid = 2;
    % % tStart = 140;
    % % tEnd = 250;
    % % tFrequency = 10;
    % % nproc = 128;
    % % nd = app.nd;
    % % if rem(tEnd - tStart,tFrequency) ~= 0; error('Start and end time-steps not consistent with time-step frequency'); end
    % % nt = (tEnd - tStart)/tFrequency + 1;
    % % if strcmp(mesh.hybrid,'hdg'); UH = reshape(UH, [app.nch size(mesh.perm,1) mesh.nf]); end
    % % timeStepNo = 0;
    % % Cp_avg = 0;
    % % Cf_avg = 0;
    % % Cp2d_avg = 0;
    % % Cf2d_avg = 0;
    % % Cp = cell(0);
    % % Cf = cell(0);
    % % x = cell(0);
    % % if nd == 3
    % %   Cp2d = cell(0);
    % %   Cf2d = cell(0);
    % %   x2d = cell(0);
    % % end
    % % for n = tStart:tFrequency:tEnd
    % %     timeStepNo = timeStepNo + 1;
    % %     disp(['Time-step No. ',num2str(timeStepNo),' / ',num2str(nt)]);
    % %     for i = 1:nproc
    % %         filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA45/3DLESIEDG2_t' num2str(n) '_np' num2str(i-1)];
    % %         fileID = fopen([filename,'.bin'],'r');
    % %         data = fread(fileID,'double');    
    % %         fclose(fileID);
    % %         
    % %         nei = sum(dmd{i}.elempartpts(1:2));
    % %         nfi = sum(dmd{i}.entpartpts(1:2)); 
    % %         inde = dmd{i}.elempart(1:nei);               
    % %         indf = dmd{i}.entpart(1:nfi);
    % %         
    % %         npf = size(mesh.perm,1);
    % %         nc = app.nc;    
    % %         nch = app.nch;    
    % %         npe = size(mesh.plocal,1);        
    % %         n1 = npe*nc*nei;        
    % %         
    % %         UDG(:,:,inde) = reshape(data(1:n1),[npe nc nei]);       
    % %         if strcmp(app.hybrid,'hdg')
    % %             UH(:,:,indf) = reshape(data(n1+1:end),[nch npf nfi]);    
    % %         elseif strcmp(app.hybrid,'edg') || strcmp(app.hybrid,'iedg')
    % %             UH(:,indf) = reshape(data(n1+1:end),[nch nfi]);
    % %         else
    % %             error('Invalid hybrid flag.');
    % %         end
    % %     end
    % %     if nd ==2
    % %         [Cp{timeStepNo},Cf{timeStepNo},x{timeStepNo}] = getsurfacedata(master,mesh,app,UDG,UH,wid);
    % %     elseif nd == 3
    % %         [Cp{timeStepNo},Cf{timeStepNo},x{timeStepNo},Cp2d{timeStepNo},Cf2d{timeStepNo},x2d{timeStepNo}] = getsurfacedata(master,mesh,app,UDG,UH,wid);
    % %     end
    % %     if timeStepNo == 1
    % %         Cp_avg = zeros(size(Cp{1}));
    % %         Cf_avg = zeros(size(Cf{1}));
    % %         if nd == 3
    % %             Cp2d_avg = zeros(size(Cp2d{1}));
    % %             Cf2d_avg = zeros(size(Cf2d{1}));
    % %         end
    % %     end
    % %     Cp_avg = Cp_avg + Cp{timeStepNo};
    % %     Cf_avg = Cf_avg + Cf{timeStepNo};
    % %     if nd == 3
    % %         Cp2d_avg = Cp2d_avg + Cp2d{timeStepNo};
    % %         Cf2d_avg = Cf2d_avg + Cf2d{timeStepNo};
    % %     end
    % % end
    % % Cp_avg = Cp_avg / nt;
    % % Cf_avg = Cf_avg / nt;
    % % if nd == 3
    % %     Cp2d_avg = Cp2d_avg / nt;
    % %     Cf2d_avg = Cf2d_avg / nt;
    % % end
    % 
    % function [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata(master,mesh,app,UDG,UH,wid,elemAvg,deltaT)
    % 
    % if nargin < 7; elemAvg = 0; end
    % if nargin < 8; deltaT = 1; end
    % fhat   = str2func(app.fhat);
    % 
    % gam = app.arg{1};
    % % Re = app.arg{3};
    % % Pr = app.arg{4};
    % Minf = app.arg{5};
    % pinf = 1/(gam*Minf^2);
    % % pinf = 1.69061708; % For T106C with Min = 0.65
    % % pinf = 1.269841269841270; % For T106C with Min = 0.75
    % 
    % nd  = master.nd;
    % ngf = master.ngf;
    % npf = master.npf;
    % ne  = size(mesh.t,1);
    % nfe = size(master.perm,2);
    % 
    % elcon  = mesh.elcon;
    % perm   = master.perm(:,:,1);
    % shapft = squeeze(master.shapft(:,:,1));
    % dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
    % 
    % xc = [];
    % yc = [];
    % zc = [];
    % cpData = [];
    % cfData = [];
    % chData = [];
    % 
    % for i=1:ne
    %     dg = mesh.dgnodes(:,:,i);       
    %     uh = UH(:,elcon(:,i))';          % uh(3*nps,nc)
    %     udg = UDG(:,:,i);
    % 
    % 
    %     fc=mesh.f(abs(mesh.t2f(i,:)),end);  % fc=mesh.f(abs(mesh.t2f(i,:)),4);
    %     bf  = (fc<0);
    % 
    %     for is = 1:nfe
    %         if bf(is)==1 % boundary face
    % %             ib = app.bcm(-fc(is));
    %             if fc(is) == -wid %ib == wid
    % 
    %                 pn    = dg(perm(:,is),:);
    %                 pg    = shapft*pn;
    %                 udgg  = shapft*udg(perm(:,is),:);
    %                 % uhg   = shapft*uh((is-1)*npf+1:is*npf,:);   
    %                 uhg = udgg(:,1:size(UH,1));
    %                 dpg   = dshapft*pn(:,1:nd);
    %                 dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  
    % 
    %                 if nd==2
    %                     jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
    %                     nlg   = [dpg(:,2),-dpg(:,1)];
    %                     nlg   = bsxfun(@rdivide, nlg, jac);
    %                 elseif nd==3
    %                     nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
    %                     nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
    %                     nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
    %                     jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
    %                     nlg   = bsxfun(@rdivide, nlg, jac);
    %                 end             
    % 
    %                 [fh] = fhat(nlg,pg,udgg,uhg,app.arg,0);                                 
    % 
    %                 pr = fh(:,2).*nlg(:,1) + fh(:,3).*nlg(:,2);
    %                 if any(abs(nlg(:,2)) < 1.0e-2)
    %                     1;
    %                 end
    % %                 pres = eulereval(udgg,'p',1.4);
    % %                 prew = eulereval(uhg,'p',1.4);
    % %                 [pr pres prew] 
    % %                 [udgg(:,2:3) uhg(:,2:3)]
    % %                 pause
    % 
    %                 q = fh(:,nd+2);
    %                 if nlg(:,2) > 0     % This is only an approximation for the boundary layer side
    %                     tau = fh(:,2).*nlg(:,2) - fh(:,3).*nlg(:,1);
    %                 else
    %                     tau = - fh(:,2).*nlg(:,2) + fh(:,3).*nlg(:,1);
    %                 end
    % 
    %                 if elemAvg == 0
    %                     xNew = pg(:,1);
    %                     yNew = pg(:,2);
    %                     if nd == 3
    %                         zNew = pg(:,3);
    %                     end
    %                     cpNew = pr;
    %                     cfNew = tau;
    %                     chNew = q;
    %                 elseif elemAvg == 1
    %                     xNew = sum(abs(jac).*master.gwfc.*pg(:,1))/sum(abs(jac).*master.gwfc);
    %                     yNew = sum(abs(jac).*master.gwfc.*pg(:,2))/sum(abs(jac).*master.gwfc);
    %                     if nd == 3
    %                         zNew = sum(abs(jac).*master.gwfc.*pg(:,3))/sum(abs(jac).*master.gwfc);
    %                     end
    %                     cpNew = sum(abs(jac).*master.gwfc.*pr)/sum(abs(jac).*master.gwfc);
    %                     cfNew = sum(abs(jac).*master.gwfc.*tau)/sum(abs(jac).*master.gwfc);
    %                     chNew = sum(abs(jac).*master.gwfc.*q)/sum(abs(jac).*master.gwfc);
    %                 else
    %                     error('elemAvg has invalid value');
    %                 end
    % 
    %                 xc = [xc; xNew];
    %                 yc = [yc; yNew];
    %                 if nd == 3
    %                     zc = [zc; zNew];
    %                 end
    %                 cpData = [cpData; cpNew];
    %                 cfData = [cfData; cfNew];
    %                 chData = [chData; chNew];
    %             end
    %         end
    %     end
    % end
    % 
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%%%%%%%%%%%%%%%%%%%%%%%  AIRFOIL DEPENDENT %%%%%%%%%%%%%%%%%%%%%%%%
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % % xcmax = max(xc);
    % % % % % % % ind = yc < -(1/xcmax^2)*0.1*xc.*xc+0.1;    % Trefftz foil dividing curve
    % % % 
    % % % NACA 65-(18)10:
    % % % alpha = (pi/180)*16.73;
    % % % theta = (pi/180)*45 - alpha;
    % % theta = 0;      % This for Eppler 387
    % % [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,0,0,1,18/10,1);
    % % ind = yWithNoRotation < yMeanLine;       % % NACA 65-(18)10
    % % % 
    % % ind = yc < 0;                             % $ Symemtric foil
    % %ind = (yc < 0.0015*xc.^2 & xc<0.99) | (xc>=0.99 & yc<=-0.1*xc+0.1);    
    % xWithNoRotation = xc;
    % yWithNoRotation = yc;
    % ind = yWithNoRotation < 0;
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % T106C:
    % % theta = -0.527999463;
    % % [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,0,0,1);
    % % ind = yWithNoRotation < yMeanLine;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % Cylinder
    % % xWithNoRotation = xc;
    % % yWithNoRotation = yc;
    % % ind = yWithNoRotation < 0;
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % xclower = xWithNoRotation(ind);
    % yclower = yWithNoRotation(ind);
    % if nd == 3; zclower = zc(ind); end
    % datalower = cpData(ind);
    % data1lower = cfData(ind);
    % data2lower = chData(ind);
    % 
    % xcupper = xWithNoRotation(~ind);
    % ycupper = yWithNoRotation(~ind);
    % if nd == 3; zcupper = zc(~ind); end
    % dataupper = cpData(~ind);
    % data1upper = cfData(~ind);
    % data2upper = chData(~ind);
    % 
    % [xcl,ij] = sort(xclower);
    % ycl = yclower(ij);
    % if nd == 3; zcl = zclower(ij); end
    % datal = datalower(ij);
    % data1l = data1lower(ij);
    % data2l = data2lower(ij);
    % 
    % [xcu,ij] = sort(xcupper);
    % ycu = ycupper(ij);
    % if nd == 3; zcu = zcupper(ij); end
    % datau = dataupper(ij);
    % data1u = data1upper(ij);
    % data2u = data2upper(ij);
    % 
    % % x = [flipud(xcl); xcu];
    % % xmin = min(x);
    % % xmax = max(x);
    % % x = (x-xmin)/(xmax-xmin);
    % % y = [flipud(ycl); ycu]/(xmax-xmin);
    % % z = [flipud(zcl); zcu]/(xmax-xmin);
    % x = [flipud(xcl); xcu];
    % y = [flipud(ycl); ycu];
    % if nd == 3; z = [flipud(zcl); zcu]; end 
    % if nd == 2
    %     x = [x(:) y(:)];
    % elseif nd == 3
    %     x = [x(:) y(:) z(:)];
    % end
    % 
    % Cp = -2*([flipud(datal); datau]-0);
    % Cf = -2*[flipud(data1l); -data1u];
    % Ch = [flipud(data2l); data2u] / (deltaT*gam);
    % 
    % if nd == 3
    %     [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch);
    % else
    %     Cp2d = [];
    %     Cf2d = [];
    %     x2d = [];
    %     Ch2d = [];
    % end
    % 
    % % figure(1)
    % % plot(x,Cp,'-'); 
    % % grid on;
    % % axis equal;
    % % axis tight;
    % % axis on;
    % % 
    % % figure(2) 
    % % plot(x,Cf); grid on;
    % % %axis square;
    % % axis tight;
    % % axis on;
    % 
    % end
    % 
    % function [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch)
    % % z = x(:,3);
    % % z = z(:);snap = 1.0e-8;
    % snap = 1.0e-9;
    % xyPoints = unique(round(squeeze(x(:,1:2))/snap)*snap,'rows');
    % num_xyPoints = size(xyPoints,1);
    % Cp2d = zeros(num_xyPoints,1);
    % Cf2d = zeros(num_xyPoints,1);
    % Ch2d = zeros(num_xyPoints,1);
    % 
    % tol = 1.0e-5;
    % for i=1:num_xyPoints
    %     indices = find(sqrt((x(:,1)-xyPoints(i,1)).^2 + (x(:,2)-xyPoints(i,2)).^2) < tol);
    %     Cp2d(i) = sum(Cp(indices))/length(indices);
    %     Cf2d(i) = sum(Cf(indices))/length(indices);
    %     Ch2d(i) = sum(Ch(indices))/length(indices);
    % end
    % % % NACA 65-(18)10
    % % [yMeanLine,~,yWithNoRotation] = getNACA65_1810meanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1,18/10,1);
    % % ind = yWithNoRotation < yMeanLine;
    % % % T106C
    % % [yMeanLine,~,yWithNoRotation] = getT106CmeanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1);
    % % ind = yWithNoRotation < yMeanLine;
    % % % Cylinder
    % ind = xyPoints(:,2) < 0;
    % 
    % xclower2d = xyPoints(ind,1);
    % yclower2d = xyPoints(ind,2);
    % Cp2dlower = Cp2d(ind);
    % Cf2dlower = Cf2d(ind);
    % Ch2dlower = Ch2d(ind);
    % 
    % xcupper2d = xyPoints(~ind,1);
    % ycupper2d = xyPoints(~ind,2);
    % Cp2dupper = Cp2d(~ind);
    % Cf2dupper = Cf2d(~ind);
    % Ch2dupper = Ch2d(~ind);
    % 
    % [xcl2d,ij] = sort(xclower2d);
    % ycl2d = yclower2d(ij);
    % Cp2dlower = Cp2dlower(ij);
    % Cf2dlower = Cf2dlower(ij);
    % Ch2dlower = Ch2dlower(ij);
    % 
    % [xcu2d,ij] = sort(xcupper2d);
    % ycu2d = ycupper2d(ij);
    % Cp2dupper = Cp2dupper(ij);
    % Cf2dupper = Cf2dupper(ij);
    % Ch2dupper = Ch2dupper(ij);
    % 
    % x2d = [flipud(xcl2d); xcu2d];
    % y2d = [flipud(ycl2d); ycu2d];
    % x2d = [x2d(:) y2d(:)];
    % Cp2d = [flipud(Cp2dlower); Cp2dupper];
    % Cf2d = [flipud(Cf2dlower); Cf2dupper];
    % Ch2d = [flipud(Ch2dlower); Ch2dupper];
    % 
    % end
    % 
    % function [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,xTranslation,yTranslation,scaling,cl,a)
    % 
    % if nargin < 5; a = 1; end
    % if a == 1; a = 0.999999; end
    % 
    % A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
    % xyzWithNoRotation = [xc(:) yc(:)];
    % xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
    % xWithNoRotation = xyzWithNoRotation(:,1);
    % yWithNoRotation = xyzWithNoRotation(:,2);
    % xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
    % yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
    % g = -(1/(1-a)) * (a^2*(0.5*log(a)-0.25) + 0.25);
    % h = (1/(1-a)) * (0.5*(1-a)^2*log(1-a) - 0.25*(1-a)^2) + g;
    % yMeanLine = (cl/(2*pi*(a+1))) * ( (1/(1-a))* (0.5*(a-xWithNoRotation).^2.*log(abs(a-xWithNoRotation)) - 0.5*(1-xWithNoRotation).^2.*log(1-xWithNoRotation) + ...
    %         0.25*(1-xWithNoRotation).^2 - 0.25*(a-xWithNoRotation).^2) - xWithNoRotation.*log(xWithNoRotation) + g - h*xWithNoRotation);
    % 
    % end
    % 
    % function [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,xTranslation,yTranslation,scaling)
    % 
    % A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
    % xyzWithNoRotation = [xc(:) yc(:)];
    % xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
    % xWithNoRotation = xyzWithNoRotation(:,1);
    % yWithNoRotation = xyzWithNoRotation(:,2);
    % xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
    % yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
    % 
    % % xyMidChord = [-0.004850771437097   0.008104085834677
    % %    0.026267527614243   0.040006331527312
    % %    0.058813891977280   0.154717136626771
    % %    0.119764563959287   0.220726452638423
    % %    0.187900992975752   0.274415205427912
    % %    0.272808035650415   0.299349380953031
    % %    0.368800317722688   0.305277014484574
    % %    0.473686474141603   0.295955379128744
    % %    0.600733470684711   0.248637186068960
    % %    0.747001884020969   0.168362315002401
    % %    0.931784011091559   0.022052558159896
    % %    0.989643954790913   0.008104085275412];
    % 
    % % xyMidChord = [0.003   0.003
    % %    0.011900899226019   0.022010534473781
    % %    0.022099897892833   0.089780617294087
    % %    0.043846090146471   0.137752100197530
    % %    0.073184363539016   0.172706344247624
    % %    0.103314843692981   0.206302285554901
    % %    0.136214284225593   0.235150619721944
    % %    0.171832517935802   0.259337362483664
    % %    0.203472806334644   0.290344615105094
    % %    0.250885642315712   0.294308555695582
    % %    0.294730428985119   0.304390206210480
    % %    0.346638258603599   0.300647115652073
    % %    0.391564089454086   0.308875228242563
    % %    0.445012961862948   0.302489894836824
    % %    0.501235349883211   0.291349143974849
    % %    0.565068672054260   0.267158829241619
    % %    0.634686967796641   0.233049708141524
    % %    0.707440700461415   0.193564626144587
    % %    0.788313850923384   0.140158143825069
    % %    0.867798377388409   0.089132569969907
    % %    0.965128210899839   0.007509769875418
    % %    0.982   -4.8e-3];
    % 
    % xyMidChord = [-0.004850771437097   0.008104085834677
    %    0.011900899226019   0.022010534473781
    %    0.022099897892833   0.089780617294087
    %    0.043846090146471   0.137752100197530
    %    0.073184363539016   0.172706344247624
    %    0.103314843692981   0.206302285554901
    %    0.136214284225593   0.235150619721944
    %    0.171832517935802   0.259337362483664
    %    0.203472806334644   0.290344615105094
    %    0.250885642315712   0.294308555695582
    %    0.294730428985119   0.304390206210480
    %    0.346638258603599   0.300647115652073
    %    0.391564089454086   0.308875228242563
    %    0.445012961862948   0.302489894836824
    %    0.501235349883211   0.291349143974849
    %    0.565068672054260   0.267158829241619
    %    0.634686967796641   0.233049708141524
    %    0.707440700461415   0.193564626144587
    %    0.788313850923384   0.140158143825069
    %    0.867798377388409   0.089132569969907
    %    0.965128210899839   0.007509769875418
    %    0.989643954790913   0.008104085275412];
    % 
    % yMeanLine = interp1(xyMidChord(:,1),xyMidChord(:,2),xWithNoRotation);
    % 
    % % % Mid-chord (x,y) coordinates with rotation:
    % %   -0.000107288164668   0.009444293515160
    % %    0.042845768527915   0.021324391866838
    % %    0.128751881913082   0.104016436967154
    % %    0.214657995298249   0.130329055583308
    % %    0.300564108683416   0.142378736834755
    % %    0.386470222068583   0.121140594544775
    % %    0.472376335453749   0.077899471296392
    % %    0.558282448838916   0.017004985691412
    % %    0.644188562224083  -0.087876330128531
    % %    0.730094675609250  -0.230910044153280
    % %    0.816000788994417  -0.450389243885509
    % %    0.858953845687000  -0.491588368894000
    % %     
    % % % Mid-chord (x,y) coordinates with no rotation:
    % %     -0.004850771437097   0.008104085834677
    % %    0.026267527614243   0.040006331527312
    % %    0.058813891977280   0.154717136626771
    % %    0.119764563959287   0.220726452638423
    % %    0.187900992975752   0.274415205427912
    % %    0.272808035650415   0.299349380953031
    % %    0.368800317722688   0.305277014484574
    % %    0.473686474141603   0.295955379128744
    % %    0.600733470684711   0.248637186068960
    % %    0.747001884020969   0.168362315002401
    % %    0.931784011091559   0.022052558159896
    % %    0.989643954790913   0.008104085275412
    % 
    % end
    

% % setapplicationpath('FM/ns');
% % wid = 2;
% % tStart = 140;
% % tEnd = 250;
% % tFrequency = 10;
% % nproc = 128;
% % nd = app.nd;
% % if rem(tEnd - tStart,tFrequency) ~= 0; error('Start and end time-steps not consistent with time-step frequency'); end
% % nt = (tEnd - tStart)/tFrequency + 1;
% % if strcmp(mesh.hybrid,'hdg'); UH = reshape(UH, [app.nch size(mesh.perm,1) mesh.nf]); end
% % timeStepNo = 0;
% % Cp_avg = 0;
% % Cf_avg = 0;
% % Cp2d_avg = 0;
% % Cf2d_avg = 0;
% % Cp = cell(0);
% % Cf = cell(0);
% % x = cell(0);
% % if nd == 3
% %   Cp2d = cell(0);
% %   Cf2d = cell(0);
% %   x2d = cell(0);
% % end
% % for n = tStart:tFrequency:tEnd
% %     timeStepNo = timeStepNo + 1;
% %     disp(['Time-step No. ',num2str(timeStepNo),' / ',num2str(nt)]);
% %     for i = 1:nproc
% %         filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA45/3DLESIEDG2_t' num2str(n) '_np' num2str(i-1)];
% %         fileID = fopen([filename,'.bin'],'r');
% %         data = fread(fileID,'double');    
% %         fclose(fileID);
% %         
% %         nei = sum(dmd{i}.elempartpts(1:2));
% %         nfi = sum(dmd{i}.entpartpts(1:2)); 
% %         inde = dmd{i}.elempart(1:nei);               
% %         indf = dmd{i}.entpart(1:nfi);
% %         
% %         npf = size(mesh.perm,1);
% %         nc = app.nc;    
% %         nch = app.nch;    
% %         npe = size(mesh.plocal,1);        
% %         n1 = npe*nc*nei;        
% %         
% %         UDG(:,:,inde) = reshape(data(1:n1),[npe nc nei]);       
% %         if strcmp(app.hybrid,'hdg')
% %             UH(:,:,indf) = reshape(data(n1+1:end),[nch npf nfi]);    
% %         elseif strcmp(app.hybrid,'edg') || strcmp(app.hybrid,'iedg')
% %             UH(:,indf) = reshape(data(n1+1:end),[nch nfi]);
% %         else
% %             error('Invalid hybrid flag.');
% %         end
% %     end
% %     if nd ==2
% %         [Cp{timeStepNo},Cf{timeStepNo},x{timeStepNo}] = getsurfacedata(master,mesh,app,UDG,UH,wid);
% %     elseif nd == 3
% %         [Cp{timeStepNo},Cf{timeStepNo},x{timeStepNo},Cp2d{timeStepNo},Cf2d{timeStepNo},x2d{timeStepNo}] = getsurfacedata(master,mesh,app,UDG,UH,wid);
% %     end
% %     if timeStepNo == 1
% %         Cp_avg = zeros(size(Cp{1}));
% %         Cf_avg = zeros(size(Cf{1}));
% %         if nd == 3
% %             Cp2d_avg = zeros(size(Cp2d{1}));
% %             Cf2d_avg = zeros(size(Cf2d{1}));
% %         end
% %     end
% %     Cp_avg = Cp_avg + Cp{timeStepNo};
% %     Cf_avg = Cf_avg + Cf{timeStepNo};
% %     if nd == 3
% %         Cp2d_avg = Cp2d_avg + Cp2d{timeStepNo};
% %         Cf2d_avg = Cf2d_avg + Cf2d{timeStepNo};
% %     end
% % end
% % Cp_avg = Cp_avg / nt;
% % Cf_avg = Cf_avg / nt;
% % if nd == 3
% %     Cp2d_avg = Cp2d_avg / nt;
% %     Cf2d_avg = Cf2d_avg / nt;
% % end


% function [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata(master,mesh,app,UDG,UH,wid,elemAvg,deltaT)


% if nargin < 7; elemAvg = 0; end
% if nargin < 8; deltaT = 1; end
% % fhat   = str2func(app.fhat);


% gam = app.physicsparam(1);
% % Re = app.arg{3};
% % Pr = app.arg{4};
% Minf = app.physicsparam(4);
% pinf = 1/(gam*Minf^2);
% % pinf = 1.69061708; % For T106C with Min = 0.65
% % pinf = 1.269841269841270; % For T106C with Min = 0.75


% nd  = master.nd;
% ngf = master.ngf;
% npf = master.npf;
% ne  = size(mesh.t,1);
% nfe = size(master.perm,2);


% elcon  = mesh.elcon;
% perm   = master.perm(:,:,1);
% shapft = squeeze(master.shapft(:,:,1));
% dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);


% xc = [];
% yc = [];
% zc = [];
% cpData = [];
% cfData = [];
% chData = [];

% ftmp = mkt2f(mesh.t',1);
% mesh.f = setbndnbrs(mesh.p', ftmp, mesh.bndexpr);

% for i=1:ne
%     dg = mesh.dgnodes(:,:,i);       
%     uh = UH(:,elcon(:,i))';          % uh(3*nps,nc)
%     udg = UDG(:,:,i);
    
    
%     fc = mesh.f(:,i);
    
%     for is = 1:nfe
%         if fc(is) == wid
                
%                 pn    = dg(perm(:,is),:);
%                 pg    = shapft*pn;
%                 udgg  = shapft*udg(perm(:,is),:);
%                 uhg   = shapft*uh((is-1)*npf+1:is*npf,:);   
% %                 uhg = udgg(:,1:size(UH,1));
%                 dpg   = dshapft*pn(:,1:nd);
%                 dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  
                
%                 if nd==2
%                     jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
%                     nlg   = [dpg(:,2),-dpg(:,1)];
%                     nlg   = bsxfun(@rdivide, nlg, jac);
%                 elseif nd==3
%                     nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
%                     nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
%                     nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
%                     jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
%                     nlg   = bsxfun(@rdivide, nlg, jac);
%                 end             
                
%                 [fh] = fhat2(nlg,pg,udgg,uhg,app.arg,0);                                 
%                 % fh = zeros(size(nlg,1),3);
%                 pr = fh(:,2).*nlg(:,1) + fh(:,3).*nlg(:,2);
%                 if any(abs(nlg(:,2)) < 1.0e-2)
%                     1;
%                 end
% %                 pres = eulereval(udgg,'p',1.4);
% %                 prew = eulereval(uhg,'p',1.4);
% %                 [pr pres prew] 
% %                 [udgg(:,2:3) uhg(:,2:3)]
% %                 pause
                
%                 q = fh(:,nd+2);
%                 if nlg(:,2) > 0     % This is only an approximation for the boundary layer side
%                     tau = fh(:,2).*nlg(:,2) - fh(:,3).*nlg(:,1);
%                 else
%                     tau = - fh(:,2).*nlg(:,2) + fh(:,3).*nlg(:,1);
%                 end
                
%                 if elemAvg == 0
%                     xNew = pg(:,1);
%                     yNew = pg(:,2);
%                     if nd == 3
%                         zNew = pg(:,3);
%                     end
%                     cpNew = pr;
%                     cfNew = tau;
%                     chNew = q;
%                 elseif elemAvg == 1
%                     xNew = sum(abs(jac).*master.gwfc.*pg(:,1))/sum(abs(jac).*master.gwfc);
%                     yNew = sum(abs(jac).*master.gwfc.*pg(:,2))/sum(abs(jac).*master.gwfc);
%                     if nd == 3
%                         zNew = sum(abs(jac).*master.gwfc.*pg(:,3))/sum(abs(jac).*master.gwfc);
%                     end
%                     cpNew = sum(abs(jac).*master.gwfc.*pr)/sum(abs(jac).*master.gwfc);
%                     cfNew = sum(abs(jac).*master.gwfc.*tau)/sum(abs(jac).*master.gwfc);
%                     chNew = sum(abs(jac).*master.gwfc.*q)/sum(abs(jac).*master.gwfc);
%                 else
%                     error('elemAvg has invalid value');
%                 end
                
%                 xc = [xc; xNew];
%                 yc = [yc; yNew];
%                 if nd == 3
%                     zc = [zc; zNew];
%                 end
%                 cpData = [cpData; cpNew];
%                 cfData = [cfData; cfNew];
%                 chData = [chData; chNew];
%             end
%         end
%     end



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%  AIRFOIL DEPENDENT %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % xcmax = max(xc);
% % % % % % % ind = yc < -(1/xcmax^2)*0.1*xc.*xc+0.1;    % Trefftz foil dividing curve
% % % 
% % % NACA 65-(18)10:
% % % alpha = (pi/180)*16.73;
% % % theta = (pi/180)*45 - alpha;
% % theta = 0;      % This for Eppler 387
% % [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,0,0,1,18/10,1);
% % ind = yWithNoRotation < yMeanLine;       % % NACA 65-(18)10
% % % 
% % ind = yc < 0;                             % $ Symemtric foil
% %ind = (yc < 0.0015*xc.^2 & xc<0.99) | (xc>=0.99 & yc<=-0.1*xc+0.1);    
% xWithNoRotation = xc;
% yWithNoRotation = yc;
% ind = yWithNoRotation < 0;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % T106C:
% % theta = -0.527999463;
% % [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,0,0,1);
% % ind = yWithNoRotation < yMeanLine;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Cylinder
% % xWithNoRotation = xc;
% % yWithNoRotation = yc;
% % ind = yWithNoRotation < 0;
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% xclower = xWithNoRotation(ind);
% yclower = yWithNoRotation(ind);
% if nd == 3; zclower = zc(ind); end
% datalower = cpData(ind);
% data1lower = cfData(ind);
% data2lower = chData(ind);


% xcupper = xWithNoRotation(~ind);
% ycupper = yWithNoRotation(~ind);
% if nd == 3; zcupper = zc(~ind); end
% dataupper = cpData(~ind);
% data1upper = cfData(~ind);
% data2upper = chData(~ind);


% [xcl,ij] = sort(xclower);
% ycl = yclower(ij);
% if nd == 3; zcl = zclower(ij); end
% datal = datalower(ij);
% data1l = data1lower(ij);
% data2l = data2lower(ij);


% [xcu,ij] = sort(xcupper);
% ycu = ycupper(ij);
% if nd == 3; zcu = zcupper(ij); end
% datau = dataupper(ij);
% data1u = data1upper(ij);
% data2u = data2upper(ij);


% % x = [flipud(xcl); xcu];
% % xmin = min(x);
% % xmax = max(x);
% % x = (x-xmin)/(xmax-xmin);
% % y = [flipud(ycl); ycu]/(xmax-xmin);
% % z = [flipud(zcl); zcu]/(xmax-xmin);
% x = [flipud(xcl); xcu];
% y = [flipud(ycl); ycu];
% if nd == 3; z = [flipud(zcl); zcu]; end 
% if nd == 2
%     x = [x(:) y(:)];
% elseif nd == 3
%     x = [x(:) y(:) z(:)];
% end


% Cp = -2*([flipud(datal); datau]-pinf);
% Cf = -2*[flipud(data1l); data1u];
% Ch = [flipud(data2l); data2u] / (deltaT*gam);


% if nd == 3
%     [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch);
% else
%     Cp2d = [];
%     Cf2d = [];
%     x2d = [];
%     Ch2d = [];
% end


% % figure(1)
% % plot(x,Cp,'-'); 
% % grid on;
% % axis equal;
% % axis tight;
% % axis on;
% % 
% % figure(2) 
% % plot(x,Cf); grid on;
% % %axis square;
% % axis tight;
% % axis on;


% end


% function [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch)
% % z = x(:,3);
% % z = z(:);snap = 1.0e-8;
% snap = 1.0e-9;
% xyPoints = unique(round(squeeze(x(:,1:2))/snap)*snap,'rows');
% num_xyPoints = size(xyPoints,1);
% Cp2d = zeros(num_xyPoints,1);
% Cf2d = zeros(num_xyPoints,1);
% Ch2d = zeros(num_xyPoints,1);


% tol = 1.0e-5;
% for i=1:num_xyPoints
%     indices = find(sqrt((x(:,1)-xyPoints(i,1)).^2 + (x(:,2)-xyPoints(i,2)).^2) < tol);
%     Cp2d(i) = sum(Cp(indices))/length(indices);
%     Cf2d(i) = sum(Cf(indices))/length(indices);
%     Ch2d(i) = sum(Ch(indices))/length(indices);
% end
% % % NACA 65-(18)10
% % [yMeanLine,~,yWithNoRotation] = getNACA65_1810meanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1,18/10,1);
% % ind = yWithNoRotation < yMeanLine;
% % % T106C
% % [yMeanLine,~,yWithNoRotation] = getT106CmeanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1);
% % ind = yWithNoRotation < yMeanLine;
% % % Cylinder
% ind = xyPoints(:,2) < 0;


% xclower2d = xyPoints(ind,1);
% yclower2d = xyPoints(ind,2);
% Cp2dlower = Cp2d(ind);
% Cf2dlower = Cf2d(ind);
% Ch2dlower = Ch2d(ind);


% xcupper2d = xyPoints(~ind,1);
% ycupper2d = xyPoints(~ind,2);
% Cp2dupper = Cp2d(~ind);
% Cf2dupper = Cf2d(~ind);
% Ch2dupper = Ch2d(~ind);


% [xcl2d,ij] = sort(xclower2d);
% ycl2d = yclower2d(ij);
% Cp2dlower = Cp2dlower(ij);
% Cf2dlower = Cf2dlower(ij);
% Ch2dlower = Ch2dlower(ij);


% [xcu2d,ij] = sort(xcupper2d);
% ycu2d = ycupper2d(ij);
% Cp2dupper = Cp2dupper(ij);
% Cf2dupper = Cf2dupper(ij);
% Ch2dupper = Ch2dupper(ij);


% x2d = [flipud(xcl2d); xcu2d];
% y2d = [flipud(ycl2d); ycu2d];
% x2d = [x2d(:) y2d(:)];
% Cp2d = [flipud(Cp2dlower); Cp2dupper];
% Cf2d = [flipud(Cf2dlower); Cf2dupper];
% Ch2d = [flipud(Ch2dlower); Ch2dupper];


% end


% function [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,xTranslation,yTranslation,scaling,cl,a)


% if nargin < 5; a = 1; end
% if a == 1; a = 0.999999; end


% A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
% xyzWithNoRotation = [xc(:) yc(:)];
% xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
% xWithNoRotation = xyzWithNoRotation(:,1);
% yWithNoRotation = xyzWithNoRotation(:,2);
% xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
% yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
% g = -(1/(1-a)) * (a^2*(0.5*log(a)-0.25) + 0.25);
% h = (1/(1-a)) * (0.5*(1-a)^2*log(1-a) - 0.25*(1-a)^2) + g;
% yMeanLine = (cl/(2*pi*(a+1))) * ( (1/(1-a))* (0.5*(a-xWithNoRotation).^2.*log(abs(a-xWithNoRotation)) - 0.5*(1-xWithNoRotation).^2.*log(1-xWithNoRotation) + ...
%         0.25*(1-xWithNoRotation).^2 - 0.25*(a-xWithNoRotation).^2) - xWithNoRotation.*log(xWithNoRotation) + g - h*xWithNoRotation);


% end


% function [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,xTranslation,yTranslation,scaling)


% A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
% xyzWithNoRotation = [xc(:) yc(:)];
% xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
% xWithNoRotation = xyzWithNoRotation(:,1);
% yWithNoRotation = xyzWithNoRotation(:,2);
% xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
% yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;


% % xyMidChord = [-0.004850771437097   0.008104085834677
% %    0.026267527614243   0.040006331527312
% %    0.058813891977280   0.154717136626771
% %    0.119764563959287   0.220726452638423
% %    0.187900992975752   0.274415205427912
% %    0.272808035650415   0.299349380953031
% %    0.368800317722688   0.305277014484574
% %    0.473686474141603   0.295955379128744
% %    0.600733470684711   0.248637186068960
% %    0.747001884020969   0.168362315002401
% %    0.931784011091559   0.022052558159896
% %    0.989643954790913   0.008104085275412];


% % xyMidChord = [0.003   0.003
% %    0.011900899226019   0.022010534473781
% %    0.022099897892833   0.089780617294087
% %    0.043846090146471   0.137752100197530
% %    0.073184363539016   0.172706344247624
% %    0.103314843692981   0.206302285554901
% %    0.136214284225593   0.235150619721944
% %    0.171832517935802   0.259337362483664
% %    0.203472806334644   0.290344615105094
% %    0.250885642315712   0.294308555695582
% %    0.294730428985119   0.304390206210480
% %    0.346638258603599   0.300647115652073
% %    0.391564089454086   0.308875228242563
% %    0.445012961862948   0.302489894836824
% %    0.501235349883211   0.291349143974849
% %    0.565068672054260   0.267158829241619
% %    0.634686967796641   0.233049708141524
% %    0.707440700461415   0.193564626144587
% %    0.788313850923384   0.140158143825069
% %    0.867798377388409   0.089132569969907
% %    0.965128210899839   0.007509769875418
% %    0.982   -4.8e-3];


% xyMidChord = [-0.004850771437097   0.008104085834677
%    0.011900899226019   0.022010534473781
%    0.022099897892833   0.089780617294087
%    0.043846090146471   0.137752100197530
%    0.073184363539016   0.172706344247624
%    0.103314843692981   0.206302285554901
%    0.136214284225593   0.235150619721944
%    0.171832517935802   0.259337362483664
%    0.203472806334644   0.290344615105094
%    0.250885642315712   0.294308555695582
%    0.294730428985119   0.304390206210480
%    0.346638258603599   0.300647115652073
%    0.391564089454086   0.308875228242563
%    0.445012961862948   0.302489894836824
%    0.501235349883211   0.291349143974849
%    0.565068672054260   0.267158829241619
%    0.634686967796641   0.233049708141524
%    0.707440700461415   0.193564626144587
%    0.788313850923384   0.140158143825069
%    0.867798377388409   0.089132569969907
%    0.965128210899839   0.007509769875418
%    0.989643954790913   0.008104085275412];


% yMeanLine = interp1(xyMidChord(:,1),xyMidChord(:,2),xWithNoRotation);


% % % Mid-chord (x,y) coordinates with rotation:
% %   -0.000107288164668   0.009444293515160
% %    0.042845768527915   0.021324391866838
% %    0.128751881913082   0.104016436967154
% %    0.214657995298249   0.130329055583308
% %    0.300564108683416   0.142378736834755
% %    0.386470222068583   0.121140594544775
% %    0.472376335453749   0.077899471296392
% %    0.558282448838916   0.017004985691412
% %    0.644188562224083  -0.087876330128531
% %    0.730094675609250  -0.230910044153280
% %    0.816000788994417  -0.450389243885509
% %    0.858953845687000  -0.491588368894000
% %     
% % % Mid-chord (x,y) coordinates with no rotation:
% %     -0.004850771437097   0.008104085834677
% %    0.026267527614243   0.040006331527312
% %    0.058813891977280   0.154717136626771
% %    0.119764563959287   0.220726452638423
% %    0.187900992975752   0.274415205427912
% %    0.272808035650415   0.299349380953031
% %    0.368800317722688   0.305277014484574
% %    0.473686474141603   0.295955379128744
% %    0.600733470684711   0.248637186068960
% %    0.747001884020969   0.168362315002401
% %    0.931784011091559   0.022052558159896
% %    0.989643954790913   0.008104085275412
    
% end


% % 
% % % setapplicationpath('FM/ns');
% % % wid = 2;
% % % tStart = 140;
% % % tEnd = 250;
% % % tFrequency = 10;
% % % nproc = 128;
% % % nd = app.nd;
% % % if rem(tEnd - tStart,tFrequency) ~= 0; error('Start and end time-steps not consistent with time-step frequency'); end
% % % nt = (tEnd - tStart)/tFrequency + 1;
% % % if strcmp(mesh.hybrid,'hdg'); UH = reshape(UH, [app.nch size(mesh.perm,1) mesh.nf]); end
% % % timeStepNo = 0;
% % % Cp_avg = 0;
% % % Cf_avg = 0;
% % % Cp2d_avg = 0;
% % % Cf2d_avg = 0;
% % % Cp = cell(0);
% % % Cf = cell(0);
% % % x = cell(0);
% % % if nd == 3
% % %   Cp2d = cell(0);
% % %   Cf2d = cell(0);
% % %   x2d = cell(0);
% % % end
% % % for n = tStart:tFrequency:tEnd
% % %     timeStepNo = timeStepNo + 1;
% % %     disp(['Time-step No. ',num2str(timeStepNo),' / ',num2str(nt)]);
% % %     for i = 1:nproc
% % %         filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA45/3DLESIEDG2_t' num2str(n) '_np' num2str(i-1)];
% % %         fileID = fopen([filename,'.bin'],'r');
% % %         data = fread(fileID,'double');    
% % %         fclose(fileID);
% % %         
% % %         nei = sum(dmd{i}.elempartpts(1:2));
% % %         nfi = sum(dmd{i}.entpartpts(1:2)); 
% % %         inde = dmd{i}.elempart(1:nei);               
% % %         indf = dmd{i}.entpart(1:nfi);
% % %         
% % %         npf = size(mesh.perm,1);
% % %         nc = app.nc;    
% % %         nch = app.nch;    
% % %         npe = size(mesh.plocal,1);        
% % %         n1 = npe*nc*nei;        
% % %         
% % %         UDG(:,:,inde) = reshape(data(1:n1),[npe nc nei]);       
% % %         if strcmp(app.hybrid,'hdg')
% % %             UH(:,:,indf) = reshape(data(n1+1:end),[nch npf nfi]);    
% % %         elseif strcmp(app.hybrid,'edg') || strcmp(app.hybrid,'iedg')
% % %             UH(:,indf) = reshape(data(n1+1:end),[nch nfi]);
% % %         else
% % %             error('Invalid hybrid flag.');
% % %         end
% % %     end
% % %     if nd ==2
% % %         [Cp{timeStepNo},Cf{timeStepNo},x{timeStepNo}] = getsurfacedata(master,mesh,app,UDG,UH,wid);
% % %     elseif nd == 3
% % %         [Cp{timeStepNo},Cf{timeStepNo},x{timeStepNo},Cp2d{timeStepNo},Cf2d{timeStepNo},x2d{timeStepNo}] = getsurfacedata(master,mesh,app,UDG,UH,wid);
% % %     end
% % %     if timeStepNo == 1
% % %         Cp_avg = zeros(size(Cp{1}));
% % %         Cf_avg = zeros(size(Cf{1}));
% % %         if nd == 3
% % %             Cp2d_avg = zeros(size(Cp2d{1}));
% % %             Cf2d_avg = zeros(size(Cf2d{1}));
% % %         end
% % %     end
% % %     Cp_avg = Cp_avg + Cp{timeStepNo};
% % %     Cf_avg = Cf_avg + Cf{timeStepNo};
% % %     if nd == 3
% % %         Cp2d_avg = Cp2d_avg + Cp2d{timeStepNo};
% % %         Cf2d_avg = Cf2d_avg + Cf2d{timeStepNo};
% % %     end
% % % end
% % % Cp_avg = Cp_avg / nt;
% % % Cf_avg = Cf_avg / nt;
% % % if nd == 3
% % %     Cp2d_avg = Cp2d_avg / nt;
% % %     Cf2d_avg = Cf2d_avg / nt;
% % % end
% % 
% % function [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata(master,mesh,app,UDG,UH,wid,elemAvg,deltaT)
% % 
% % if nargin < 7; elemAvg = 0; end
% % if nargin < 8; deltaT = 1; end
% % fhat   = str2func(app.fhat);
% % 
% % gam = app.arg{1};
% % % Re = app.arg{3};
% % % Pr = app.arg{4};
% % Minf = app.arg{5};
% % pinf = 1/(gam*Minf^2);
% % % pinf = 1.69061708; % For T106C with Min = 0.65
% % % pinf = 1.269841269841270; % For T106C with Min = 0.75
% % 
% % nd  = master.nd;
% % ngf = master.ngf;
% % npf = master.npf;
% % ne  = size(mesh.t,1);
% % nfe = size(master.perm,2);
% % 
% % elcon  = mesh.elcon;
% % perm   = master.perm(:,:,1);
% % shapft = squeeze(master.shapft(:,:,1));
% % dshapft  = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
% % 
% % xc = [];
% % yc = [];
% % zc = [];
% % cpData = [];
% % cfData = [];
% % chData = [];
% % 
% % for i=1:ne
% %     dg = mesh.dgnodes(:,:,i);       
% %     uh = UH(:,elcon(:,i))';          % uh(3*nps,nc)
% %     udg = UDG(:,:,i);
% % 
% % 
% %     fc=mesh.f(abs(mesh.t2f(i,:)),end);  % fc=mesh.f(abs(mesh.t2f(i,:)),4);
% %     bf  = (fc<0);
% % 
% %     for is = 1:nfe
% %         if bf(is)==1 % boundary face
% % %             ib = app.bcm(-fc(is));
% %             if fc(is) == -wid %ib == wid
% % 
% %                 pn    = dg(perm(:,is),:);
% %                 pg    = shapft*pn;
% %                 udgg  = shapft*udg(perm(:,is),:);
% %                 % uhg   = shapft*uh((is-1)*npf+1:is*npf,:);   
% %                 uhg = udgg(:,1:size(UH,1));
% %                 dpg   = dshapft*pn(:,1:nd);
% %                 dpg   = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 2]);  
% % 
% %                 if nd==2
% %                     jac   = sqrt(dpg(:,1).^2+dpg(:,2).^2);
% %                     nlg   = [dpg(:,2),-dpg(:,1)];
% %                     nlg   = bsxfun(@rdivide, nlg, jac);
% %                 elseif nd==3
% %                     nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
% %                     nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
% %                     nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
% %                     jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
% %                     nlg   = bsxfun(@rdivide, nlg, jac);
% %                 end             
% % 
% %                 [fh] = fhat(nlg,pg,udgg,uhg,app.arg,0);                                 
% % 
% %                 pr = fh(:,2).*nlg(:,1) + fh(:,3).*nlg(:,2);
% %                 if any(abs(nlg(:,2)) < 1.0e-2)
% %                     1;
% %                 end
% % %                 pres = eulereval(udgg,'p',1.4);
% % %                 prew = eulereval(uhg,'p',1.4);
% % %                 [pr pres prew] 
% % %                 [udgg(:,2:3) uhg(:,2:3)]
% % %                 pause
% % 
% %                 q = fh(:,nd+2);
% %                 if nlg(:,2) > 0     % This is only an approximation for the boundary layer side
% %                     tau = fh(:,2).*nlg(:,2) - fh(:,3).*nlg(:,1);
% %                 else
% %                     tau = - fh(:,2).*nlg(:,2) + fh(:,3).*nlg(:,1);
% %                 end
% % 
% %                 if elemAvg == 0
% %                     xNew = pg(:,1);
% %                     yNew = pg(:,2);
% %                     if nd == 3
% %                         zNew = pg(:,3);
% %                     end
% %                     cpNew = pr;
% %                     cfNew = tau;
% %                     chNew = q;
% %                 elseif elemAvg == 1
% %                     xNew = sum(abs(jac).*master.gwfc.*pg(:,1))/sum(abs(jac).*master.gwfc);
% %                     yNew = sum(abs(jac).*master.gwfc.*pg(:,2))/sum(abs(jac).*master.gwfc);
% %                     if nd == 3
% %                         zNew = sum(abs(jac).*master.gwfc.*pg(:,3))/sum(abs(jac).*master.gwfc);
% %                     end
% %                     cpNew = sum(abs(jac).*master.gwfc.*pr)/sum(abs(jac).*master.gwfc);
% %                     cfNew = sum(abs(jac).*master.gwfc.*tau)/sum(abs(jac).*master.gwfc);
% %                     chNew = sum(abs(jac).*master.gwfc.*q)/sum(abs(jac).*master.gwfc);
% %                 else
% %                     error('elemAvg has invalid value');
% %                 end
% % 
% %                 xc = [xc; xNew];
% %                 yc = [yc; yNew];
% %                 if nd == 3
% %                     zc = [zc; zNew];
% %                 end
% %                 cpData = [cpData; cpNew];
% %                 cfData = [cfData; cfNew];
% %                 chData = [chData; chNew];
% %             end
% %         end
% %     end
% % end
% % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%  AIRFOIL DEPENDENT %%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % xcmax = max(xc);
% % % % % % % % ind = yc < -(1/xcmax^2)*0.1*xc.*xc+0.1;    % Trefftz foil dividing curve
% % % % 
% % % % NACA 65-(18)10:
% % % % alpha = (pi/180)*16.73;
% % % % theta = (pi/180)*45 - alpha;
% % % theta = 0;      % This for Eppler 387
% % % [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,0,0,1,18/10,1);
% % % ind = yWithNoRotation < yMeanLine;       % % NACA 65-(18)10
% % % % 
% % % ind = yc < 0;                             % $ Symemtric foil
% % %ind = (yc < 0.0015*xc.^2 & xc<0.99) | (xc>=0.99 & yc<=-0.1*xc+0.1);    
% % xWithNoRotation = xc;
% % yWithNoRotation = yc;
% % ind = yWithNoRotation < 0;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % T106C:
% % % theta = -0.527999463;
% % % [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,0,0,1);
% % % ind = yWithNoRotation < yMeanLine;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % Cylinder
% % % xWithNoRotation = xc;
% % % yWithNoRotation = yc;
% % % ind = yWithNoRotation < 0;
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % xclower = xWithNoRotation(ind);
% % yclower = yWithNoRotation(ind);
% % if nd == 3; zclower = zc(ind); end
% % datalower = cpData(ind);
% % data1lower = cfData(ind);
% % data2lower = chData(ind);
% % 
% % xcupper = xWithNoRotation(~ind);
% % ycupper = yWithNoRotation(~ind);
% % if nd == 3; zcupper = zc(~ind); end
% % dataupper = cpData(~ind);
% % data1upper = cfData(~ind);
% % data2upper = chData(~ind);
% % 
% % [xcl,ij] = sort(xclower);
% % ycl = yclower(ij);
% % if nd == 3; zcl = zclower(ij); end
% % datal = datalower(ij);
% % data1l = data1lower(ij);
% % data2l = data2lower(ij);
% % 
% % [xcu,ij] = sort(xcupper);
% % ycu = ycupper(ij);
% % if nd == 3; zcu = zcupper(ij); end
% % datau = dataupper(ij);
% % data1u = data1upper(ij);
% % data2u = data2upper(ij);
% % 
% % % x = [flipud(xcl); xcu];
% % % xmin = min(x);
% % % xmax = max(x);
% % % x = (x-xmin)/(xmax-xmin);
% % % y = [flipud(ycl); ycu]/(xmax-xmin);
% % % z = [flipud(zcl); zcu]/(xmax-xmin);
% % x = [flipud(xcl); xcu];
% % y = [flipud(ycl); ycu];
% % if nd == 3; z = [flipud(zcl); zcu]; end 
% % if nd == 2
% %     x = [x(:) y(:)];
% % elseif nd == 3
% %     x = [x(:) y(:) z(:)];
% % end
% % 
% % Cp = -2*([flipud(datal); datau]-0);
% % Cf = -2*[flipud(data1l); -data1u];
% % Ch = [flipud(data2l); data2u] / (deltaT*gam);
% % 
% % if nd == 3
% %     [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch);
% % else
% %     Cp2d = [];
% %     Cf2d = [];
% %     x2d = [];
% %     Ch2d = [];
% % end
% % 
% % % figure(1)
% % % plot(x,Cp,'-'); 
% % % grid on;
% % % axis equal;
% % % axis tight;
% % % axis on;
% % % 
% % % figure(2) 
% % % plot(x,Cf); grid on;
% % % %axis square;
% % % axis tight;
% % % axis on;
% % 
% % end
% % 
% % function [Cp2d,Cf2d,x2d,Ch2d] = getZaverage(Cp,Cf,x,Ch)
% % % z = x(:,3);
% % % z = z(:);snap = 1.0e-8;
% % snap = 1.0e-9;
% % xyPoints = unique(round(squeeze(x(:,1:2))/snap)*snap,'rows');
% % num_xyPoints = size(xyPoints,1);
% % Cp2d = zeros(num_xyPoints,1);
% % Cf2d = zeros(num_xyPoints,1);
% % Ch2d = zeros(num_xyPoints,1);
% % 
% % tol = 1.0e-5;
% % for i=1:num_xyPoints
% %     indices = find(sqrt((x(:,1)-xyPoints(i,1)).^2 + (x(:,2)-xyPoints(i,2)).^2) < tol);
% %     Cp2d(i) = sum(Cp(indices))/length(indices);
% %     Cf2d(i) = sum(Cf(indices))/length(indices);
% %     Ch2d(i) = sum(Ch(indices))/length(indices);
% % end
% % % % NACA 65-(18)10
% % % [yMeanLine,~,yWithNoRotation] = getNACA65_1810meanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1,18/10,1);
% % % ind = yWithNoRotation < yMeanLine;
% % % % T106C
% % % [yMeanLine,~,yWithNoRotation] = getT106CmeanLine(xyPoints(:,1),xyPoints(:,2),0,0,0,1);
% % % ind = yWithNoRotation < yMeanLine;
% % % % Cylinder
% % ind = xyPoints(:,2) < 0;
% % 
% % xclower2d = xyPoints(ind,1);
% % yclower2d = xyPoints(ind,2);
% % Cp2dlower = Cp2d(ind);
% % Cf2dlower = Cf2d(ind);
% % Ch2dlower = Ch2d(ind);
% % 
% % xcupper2d = xyPoints(~ind,1);
% % ycupper2d = xyPoints(~ind,2);
% % Cp2dupper = Cp2d(~ind);
% % Cf2dupper = Cf2d(~ind);
% % Ch2dupper = Ch2d(~ind);
% % 
% % [xcl2d,ij] = sort(xclower2d);
% % ycl2d = yclower2d(ij);
% % Cp2dlower = Cp2dlower(ij);
% % Cf2dlower = Cf2dlower(ij);
% % Ch2dlower = Ch2dlower(ij);
% % 
% % [xcu2d,ij] = sort(xcupper2d);
% % ycu2d = ycupper2d(ij);
% % Cp2dupper = Cp2dupper(ij);
% % Cf2dupper = Cf2dupper(ij);
% % Ch2dupper = Ch2dupper(ij);
% % 
% % x2d = [flipud(xcl2d); xcu2d];
% % y2d = [flipud(ycl2d); ycu2d];
% % x2d = [x2d(:) y2d(:)];
% % Cp2d = [flipud(Cp2dlower); Cp2dupper];
% % Cf2d = [flipud(Cf2dlower); Cf2dupper];
% % Ch2d = [flipud(Ch2dlower); Ch2dupper];
% % 
% % end
% % 
% % function [yMeanLine,xWithNoRotation,yWithNoRotation] = getNACA65_1810meanLine(xc,yc,theta,xTranslation,yTranslation,scaling,cl,a)
% % 
% % if nargin < 5; a = 1; end
% % if a == 1; a = 0.999999; end
% % 
% % A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
% % xyzWithNoRotation = [xc(:) yc(:)];
% % xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
% % xWithNoRotation = xyzWithNoRotation(:,1);
% % yWithNoRotation = xyzWithNoRotation(:,2);
% % xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
% % yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
% % g = -(1/(1-a)) * (a^2*(0.5*log(a)-0.25) + 0.25);
% % h = (1/(1-a)) * (0.5*(1-a)^2*log(1-a) - 0.25*(1-a)^2) + g;
% % yMeanLine = (cl/(2*pi*(a+1))) * ( (1/(1-a))* (0.5*(a-xWithNoRotation).^2.*log(abs(a-xWithNoRotation)) - 0.5*(1-xWithNoRotation).^2.*log(1-xWithNoRotation) + ...
% %         0.25*(1-xWithNoRotation).^2 - 0.25*(a-xWithNoRotation).^2) - xWithNoRotation.*log(xWithNoRotation) + g - h*xWithNoRotation);
% % 
% % end
% % 
% % function [yMeanLine,xWithNoRotation,yWithNoRotation] = getT106CmeanLine(xc,yc,theta,xTranslation,yTranslation,scaling)
% % 
% % A_ROT = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
% % xyzWithNoRotation = [xc(:) yc(:)];
% % xyzWithNoRotation = (A_ROT*xyzWithNoRotation')';
% % xWithNoRotation = xyzWithNoRotation(:,1);
% % yWithNoRotation = xyzWithNoRotation(:,2);
% % xWithNoRotation = (xWithNoRotation - xTranslation)/scaling;
% % yWithNoRotation = (yWithNoRotation - yTranslation)/scaling;
% % 
% % % xyMidChord = [-0.004850771437097   0.008104085834677
% % %    0.026267527614243   0.040006331527312
% % %    0.058813891977280   0.154717136626771
% % %    0.119764563959287   0.220726452638423
% % %    0.187900992975752   0.274415205427912
% % %    0.272808035650415   0.299349380953031
% % %    0.368800317722688   0.305277014484574
% % %    0.473686474141603   0.295955379128744
% % %    0.600733470684711   0.248637186068960
% % %    0.747001884020969   0.168362315002401
% % %    0.931784011091559   0.022052558159896
% % %    0.989643954790913   0.008104085275412];
% % 
% % % xyMidChord = [0.003   0.003
% % %    0.011900899226019   0.022010534473781
% % %    0.022099897892833   0.089780617294087
% % %    0.043846090146471   0.137752100197530
% % %    0.073184363539016   0.172706344247624
% % %    0.103314843692981   0.206302285554901
% % %    0.136214284225593   0.235150619721944
% % %    0.171832517935802   0.259337362483664
% % %    0.203472806334644   0.290344615105094
% % %    0.250885642315712   0.294308555695582
% % %    0.294730428985119   0.304390206210480
% % %    0.346638258603599   0.300647115652073
% % %    0.391564089454086   0.308875228242563
% % %    0.445012961862948   0.302489894836824
% % %    0.501235349883211   0.291349143974849
% % %    0.565068672054260   0.267158829241619
% % %    0.634686967796641   0.233049708141524
% % %    0.707440700461415   0.193564626144587
% % %    0.788313850923384   0.140158143825069
% % %    0.867798377388409   0.089132569969907
% % %    0.965128210899839   0.007509769875418
% % %    0.982   -4.8e-3];
% % 
% % xyMidChord = [-0.004850771437097   0.008104085834677
% %    0.011900899226019   0.022010534473781
% %    0.022099897892833   0.089780617294087
% %    0.043846090146471   0.137752100197530
% %    0.073184363539016   0.172706344247624
% %    0.103314843692981   0.206302285554901
% %    0.136214284225593   0.235150619721944
% %    0.171832517935802   0.259337362483664
% %    0.203472806334644   0.290344615105094
% %    0.250885642315712   0.294308555695582
% %    0.294730428985119   0.304390206210480
% %    0.346638258603599   0.300647115652073
% %    0.391564089454086   0.308875228242563
% %    0.445012961862948   0.302489894836824
% %    0.501235349883211   0.291349143974849
% %    0.565068672054260   0.267158829241619
% %    0.634686967796641   0.233049708141524
% %    0.707440700461415   0.193564626144587
% %    0.788313850923384   0.140158143825069
% %    0.867798377388409   0.089132569969907
% %    0.965128210899839   0.007509769875418
% %    0.989643954790913   0.008104085275412];
% % 
% % yMeanLine = interp1(xyMidChord(:,1),xyMidChord(:,2),xWithNoRotation);
% % 
% % % % Mid-chord (x,y) coordinates with rotation:
% % %   -0.000107288164668   0.009444293515160
% % %    0.042845768527915   0.021324391866838
% % %    0.128751881913082   0.104016436967154
% % %    0.214657995298249   0.130329055583308
% % %    0.300564108683416   0.142378736834755
% % %    0.386470222068583   0.121140594544775
% % %    0.472376335453749   0.077899471296392
% % %    0.558282448838916   0.017004985691412
% % %    0.644188562224083  -0.087876330128531
% % %    0.730094675609250  -0.230910044153280
% % %    0.816000788994417  -0.450389243885509
% % %    0.858953845687000  -0.491588368894000
% % %     
% % % % Mid-chord (x,y) coordinates with no rotation:
% % %     -0.004850771437097   0.008104085834677
% % %    0.026267527614243   0.040006331527312
% % %    0.058813891977280   0.154717136626771
% % %    0.119764563959287   0.220726452638423
% % %    0.187900992975752   0.274415205427912
% % %    0.272808035650415   0.299349380953031
% % %    0.368800317722688   0.305277014484574
% % %    0.473686474141603   0.295955379128744
% % %    0.600733470684711   0.248637186068960
% % %    0.747001884020969   0.168362315002401
% % %    0.931784011091559   0.022052558159896
% % %    0.989643954790913   0.008104085275412
% % 
% % end
