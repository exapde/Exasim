function [UDG, UH, ACG, mine, minf, ming, resid_history] = avloop(master, master_mat, mesh, mesh_mat, app, UDG0, UH0, S0, lambda, kappa, save_history)
if nargin < 9
    save_history = 0;
end
alpha = 100;
href = 1;
nd = mesh.nd;
porder = mesh.porder;

% if nd == 1
%     np1 = porder;
%     A = tensorproduct(master.plocvl, mesh.porder);    
% elseif nd==2
%     if mesh.elemtype==0
%         np1 = porder*(porder+1)/2;
%         A = koornwinder(master.plocvl, mesh.porder);    
%     else
%         np1 = porder*porder;
%         A = tensorproduct(master.plocvl, mesh.porder);
%     end
% else
%     if mesh.elemtype==0
%         np1 = porder*(porder+1)*(porder+2)/6;
%         A = koornwinder(master.plocvl, mesh.porder);    
%     else
%         np1 = porder*porder*porder;
%         A = tensorproduct(master.plocvl, mesh.porder);
%     end
% end
%[~,cgelcon,rowent2elem,~,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:nd,:),1e-8);
% he = meshsize(mesh);
% disp(min(he(:)));

%dist = mesh.dgnodes(:,3,:);

UDG = cell(length(lambda),1);
UH = cell(length(lambda),1);
ACG = cell(length(lambda),1);
mine = zeros(length(lambda),1);
minf = zeros(length(lambda),1);
ming = zeros(length(lambda),1);
resid_history = cell(length(lambda),1);
% app.max_iter = 1;
for i = 1:length(lambda)
    % if i == length(lambda)
    %     app.max_iter = 20;
    % end
    % resid_history_curr = [];
    % if i==1
    %     div = divergence(UDG0, href);
    % else
    %     div = divergence(UDG{i-1}, href);
    % end
    %s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
    % sm = squeeze(mean(div/max(div(:)),1));
    % idx = (sm > S0);
    % if max(idx)==0
    %   hs = he;
    % else
    %   hs = he(:,idx);
    % end
    % hk =  (kappa(i)*min(hs(:)))^2;
    % mesh.hk = hk;
    % %div = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
    % %s = (div).*(atan(alpha*(div))/pi + 0.5) - atan(alpha)/pi + 0.5; 
    % divmax = max(div(:));
    % s = limiting(div,0,divmax/2,alpha,0);
    % s = cgpoisson(mesh, master, s, [hk 1.0]);        
    % s = s/max(s(:));
    % %figure(4); clf; scaplot(mesh, reshape(s(mesh.t2'), master.npe, 1, mesh.ne),[],2,0); axis off;
    % a = (s-S0).*(atan(alpha*(s-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
    % a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);    
    %a = a/max(a(:));
    
%     s = (div-S0).*(atan(alpha*(div-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
%     a = cgpoisson(mesh, master, s, [hk 1.0]);        
%     a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);

    app.lambda = lambda(i);
    app.S0 = S0;
    app.kappa = kappa(i);
    if i>1
      mesh.vdg(:,1,:) = avf(mesh_mat, master_mat, app, UDG{i-1});    
    else
      mesh.vdg(:,1,:) = avf(mesh_mat, master_mat, app, UDG0);    
    end
    %mesh.dgnodes(:,nd+1,:) = lambda(i)*(dist + a);       
    %mesh.dgnodes(:,3,:) = 0.015*tanh(meshdist(mesh,2)*5);
    ACG{i} = mesh.vdg(:,1,:); 
    % figure(2); clf; scaplot(mesh, ACG{i}); axis([-0.2 1.2 -0.2 0.8]); drawnow

    if i == 1
        mesh.udg = real(UDG0);
        % [UDG{i}, UH{i}, resid_history{i}] = hdg_solve(master,mesh,app,real(UDG0),real(UH0),[],save_history);        
    else
        mesh.udg = UDG{i-1};
        % [UDG{i}, UH{i},resid_history{i}] = hdg_solve(master,mesh,app,real(UDG{i-1}),real(UH{i-1}),[],save_history);
    end
     figure(1); clf; scaplot(mesh, eulereval(mesh.udg, 'M',1.4,8),[]);
     figure(2); clf; scaplot(mesh, ACG{i});

    [app,mesh,master,dmd] = preprocessing(app,mesh);
    runcode(app, 1); % run C++ code
    UDG_i = fetchsolution(app,master,dmd, app.buildpath + '/dataout');
    fileID = fopen(app.buildpath+"/dataout/out_uhat_np0.bin",'r');
    UH_i = fread(fileID,'double');
    UH_i = reshape(UH_i, [app.ncu mesh.nf*master.npf]);
    UDG{i} = real(UDG_i);
    UH{i} = real(UH_i);
    % UH{i} = real(UH{i}); 
    % usave = UDG{i}; uhsave = UH{i}; asave=ACG{i};
    % save("av_iter_curr.mat", "usave", "uhsave","asave");

    % figure(1); clf; scaplot(mesh, UDG{i}(:,1,:),[],2); axis off;
    
%     if size(UDG{i},2)==1
%         pres = UDG{i}(:,1,:);
%     else
%         pres = eulereval(UDG{i},'p',app.arg{1},app.arg{2});
%     end
    % pres = UDG{i}(:,1,:);
    % u = A\squeeze(pres);
    % u((np1+1):end,:) = 0;
    % U = reshape(A*u,[master.npe 1 mesh.ne]);    
    % err = calerror(mesh,master,pres,U,1);
    % mine(i) = max(err(:));
    
    % pres = eulereval(UDG{i},'p',app.arg{1},app.arg{2});
    % u = A\squeeze(pres);
    % u((np1+1):end,:) = 0;
    % U = reshape(A*u,[master.npe 1 mesh.ne]);    
    % err = calerror(mesh,master,pres,U,1);
    % minf(i) = max(err(:));    
    
    % pres = eulereval(UDG{i},'M',app.arg{1},app.arg{2});
    % u = A\squeeze(pres);
    % u((np1+1):end,:) = 0;
    % U = reshape(A*u,[master.npe 1 mesh.ne]);    
    % err = calerror(mesh,master,pres,U,1);
    % ming(i) = max(err(:));    
    [i lambda(i) kappa(i)]       
    
%     UDGt = UDG; UHt=UH;
%     save tmp.mat UDGt UHt
    
    %u = A\squeeze(eulereval(UDG{i},'p',app.arg{1},app.arg{2}));
    %u = A\squeeze(UDG{i}(:,1,:));
%     u1 = u(1:np1,:);
%     e1 = sum(abs(u1),1);
%     e = sum(abs(u),1);
%     idx = e>1e-6;
%     [mine(i),jj] = min(e1(idx)./e(idx));
%     u(:,idx(jj))'
%     u1(:,idx(jj))'
       
%     figure(1); clf; scaplot(mesh, UDG{i}(:,1,:),[],2,0); axis off;
%     figure(2); clf; scaplot(mesh, U,[],2,0); axis off;    
%     figure(3); clf; scaplot(mesh, abs(UDG{1}(:,1,:)-U),[],2,0); axis off;    
%     U1 = 0*U;
%     U2 = 0*U;
%     for j = 1:mesh.ne
%         U1(:,1,j) = sum(abs(UDG{1}(:,1,j)-U(:,1,j)))/master.npe;
%         U2(:,1,j) = err(1,j);
%     end
%     figure(4); clf; scaplot(mesh, U1,[],2,0); axis off;
%     figure(5); clf; scaplot(mesh, U2,[],2,0); axis off;
%     [i lambda(i) max(err) max(U1(:))]       
end

function err = calerror(mesh,master,UDG, U0, p)

[npv, nc, ne] = size(UDG);
ngv = master.ngv;
nd  = master.nd;

shapvt    = squeeze(master.shapvt(:,:,1));
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

err = zeros(nc,ne);
for i = 1:ne
    dg = mesh.dgnodes(:,:,i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    jac = volgeom(Jg);                   
    udgg = shapvt*UDG(:,:,i);    
    udge = shapvt*U0(:,:,i);        
    for j = 1:nc
        err(j,i) = (master.gwvl.*jac)'*(abs(udgg(:,j)-udge(:,j)).^p)/sum(master.gwvl.*jac);  
    end    
end
%err  = sqrt(err);

function [jac] = volgeom(Jg)

nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;        
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);                    
    otherwise
        error('Dimension is not implemented');
end





