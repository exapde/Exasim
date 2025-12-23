pdeapp_ns; % solve the NS equations to obtain the flow field based on uniform av
mesh.dist = dist;

S0 = 0.2; lambda0 = 0.04; kappa0 = 4; gamma = 1e3;
eta = 0.9; m = 9;
lambda = ones(m,1)*lambda0;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa0;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end

for i = 1:length(lambda)
    if (i==1)
        pdeapp_hm;
    else
        pdehm.physicsparam = kappa(i)^2*3e-3;      
        div = divergence(sol, 1);
        meshhm.vdg = limiting(div,0,3,1e3,0);
        [pdehm,meshhm,masterhm,dmdhm] = preprocessing(pdehm,meshhm);
        runcode(pdehm, 1); % run C++ code
        solhm = fetchsolution(pdehm,masterhm,dmdhm,pdehm.buildpath + '/dataout');
        s = solhm(:,1,:);
        s = s/max(s(:));
        av = (s-S0).*(atan(gamma*(s-S0))/pi + 0.5) - atan(gamma)/pi + 0.5;    
        dist = tanh(mesh.dist*5);
        av = lambda(i)*(av.*dist);         
        figure(2); clf; scaplot(meshhm,av); axis on; axis equal; axis tight;
    end

    mesh.vdg(:,1,:) = av;
    mesh.udg = sol;
    [pde,mesh,master,dmd] = preprocessing(pde,mesh);
    runcode(pde, 1); % run C++ code
    sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
    figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[0 8]);
end

% fileID = fopen(pde.buildpath+"/dataout/outuhat_np0.bin",'r');
% UH = fread(fileID,'double');
% UH = UH(4:end);
% UH = reshape(UH, [pde.ncu mesh.nf*master.npf]);
% 
% wid=1;
% gam1=gam-1;
% Tinf = pinf/(gam-1);
% Ttinf = Tinf * (1 + (gam-1)/2 * Minf^2);
% TisoW = Twall/Tref * Tinf;
% deltaT = Ttinf - TisoW;
% elemAvg = 0;
% 
% mesh1 = hdgmesh(mesh,pde.porder);
% pde.arg = {gam, Minf, 0, Re, Pr, Tref, Twall, pde.tau}; 
% [Cp1,Cf1,x1,~,~,~,Ch1]=getsurfacedata2_hdgcode(master,mesh1,pde,sol,UH,wid,elemAvg,deltaT);
% theta1 = atan2(x1(:,2),-x1(:,1));
% ii = theta1>0; Cf1(ii) = -Cf1(ii);  
% 
% ExpCp=[-79.89803449878352, 0.1260397830018084
% -51.08143571799568, 0.4360019727108335
% -23.173505371641017, 0.8803550879500246
% -16.52005139279954, 0.9347197106690777
% -7.263606790410322, 1.015403583758014
% -4.895984254120992, 0.9964162419858622
% -2.651376397583445, 1.0018576360348512
% -0.068068122795987, 1.0024823277987833
% 2.3315382302288072, 1.0012987012987011
% 4.569585303846259, 1.04670392898241
% 7.04546075831716, 1.031530494821634
% 9.424564664716657, 0.9969587374650665
% 19.536371340313277, 0.884070360019727
% 26.368606653727344, 0.7807825086306099
% 33.379623301714005, 0.686059510110143
% 40.85727563489243, 0.5671214861088278
% 47.875673163664196, 0.4568633897747822
% 62.33727891528389, 0.2556962025316456];
% 
% ExpCh = [-39.48338121721711, 0.5721278800924818
% -36.33780542449722, 0.6621063541417523
% -33.402567279053216, 0.6916527146615641
% -30.272015603995886, 0.7323447341146456
% -27.242679037402144, 0.7949772781631188
% -24.13584965339097, 0.8322410906481702
% -21.227496771132607, 0.8350155465199712
% -14.874667228972823, 0.9798453320577213
% 0, 1.0031411942916366
% -4.832230686101373, 0.9652714661564218
% 1.3181686391312368, 0.9665151877541257
% 3.942644771870634, 0.8834250179382921
% 13.320066422414929, 0.8452204416806186
% 16.361264134531748, 0.8857689547955034
% 25.390758849732464, 0.7716176353344495
% 28.526054982999025, 0.72126285577613
% 34.67250059305727, 0.6612293709638841
% 37.803843011149475, 0.6634298014828988
% 40.72801075410527, 0.5387546838874272
% 43.93368301747542, 0.5104839352626963
% 46.887898995756345, 0.5009009009009009
% 49.92118927752445, 0.42978553775013956
% 53.01694825904742, 0.3975285019532807
% -0.8642821371148423, 0.9632942677190464
% -2.3129233769999207, 0.9927290121980387
% -5.58659954136904, 0.9211512397353104
% -11.97026806188882, 0.8921310691222195
% -8.991539049526871, 0.966276010523798
% -12.006642241492923, 0.9450689627680777
% -9.717441155539154, 0.9076297536474527
% -7.160178180763857, 0.9436976799808657];
% 
% figure(3); clf; hold on;
% plot(theta1*180/pi,-Cp1/max(-Cp1(:)),'-','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
% plot(ExpCp(:,1),ExpCp(:,2),'ko','LineWidth',2,'MarkerSize',8);
% set(gca,'FontSize',18); 
% axis tight; axis on; box on; 
% xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
% ylabel("$p/p_0$", 'interpreter', 'latex', 'FontSize', 20);
% leg = legend(["$\mbox{Exasim}$", "$\mbox{Experiment}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [30,10];
% axis([-90 90 0 1.05]);
% ax = gca;
% xticks([-90:30:90]);
% yticks([0:0.2:1.0]);
% %exportgraphics(ax,"nscyl8_adaptive1_pressurecoefficient.png",'Resolution',200); 
% 
% figure(4); clf; hold on;
% plot(theta1*180/pi,Ch1/max(Ch1(:)),'-','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
% plot(ExpCh(:,1),ExpCh(:,2),'ko','LineWidth',2,'MarkerSize',8);
% set(gca,'FontSize',18); 
% axis tight; axis on; box on; 
% xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
% ylabel("$q/q_0$", 'interpreter', 'latex', 'FontSize', 20);
% leg = legend(["$\mbox{Exasim}$", "$\mbox{Experiment}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [30,10];
% axis([-90 90 0 1.05]);
% ax = gca;
% xticks([-90:30:90]);
% yticks([0:0.2:1.0]);
%exportgraphics(ax,"nscyl8_adaptive1_heatflux.png",'Resolution',200); 
