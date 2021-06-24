% ileID = fopen(['nsout1_np0','_Ru.bin'],'r');
% Ru = fread(fileID,'double');
% fclose(fileID);
% fileID = fopen(['nsout_np0','_Ru.bin'],'r');
% Ru1 = fread(fileID,'double');
% fclose(fileID);
% fileID = fopen(['nsout_np1','_Ru.bin'],'r');
% Ru2 = fread(fileID,'double');
% fclose(fileID);
% Ru = reshape(Ru,64,5,[]);
% Ru1 = reshape(Ru1,64,5,[]);
% Ru2 = reshape(Ru2,64,5,[]);
% 
fileID = fopen(['nsout1_np0','_udg.bin'],'r');
udg = fread(fileID,'double'); 
fclose(fileID);

QDG = getq(master, mesh, UDG, UH);
udg = reshape(udg,size(UDG));
for i=1:mesh.ne
    e = udg(:,6:end,i)-QDG(:,:,i);
    [i max(abs(e(:)))]
    pause
end
% 
% fileID = fopen(['nsout_np0','_udg.bin'],'r');
% udg1 = fread(fileID,'double');
% fclose(fileID);
% fileID = fopen(['nsout_np1','_udg.bin'],'r');
% udg2 = fread(fileID,'double');
% fclose(fileID);
% udg = reshape(udg,64,20,[]);
% udg1 = reshape(udg1,64,20,[]);
% udg2 = reshape(udg2,64,20,[]);
% 
% m=5;
% for j=1:length(elempart)
% for i=1:length(elempart{j})
%     if j==1
%         e = udg1(:,1:m,i)-udg(:,1:m,elempart{j}(i));
%     elseif j==2
%         e = udg2(:,1:m,i)-udg(:,1:m,elempart{j}(i));
%     end
%     [j i max(abs(e(:)))]
% end
% end

master  = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);
UH=inituhat(master,mesh.elcon,UDG,app.ncu);
%QDG = getq(master, mesh, UDG, UH);

npe = master.npv;
npf = master.npf;
ngf = master.ngf;
nd = app.nd;
ncu = app.ncu;
ncq = app.ncq;

fileID = fopen('uh.bin','r');
suh = fread(fileID,'double'); suh = reshape(suh,npf,ncu,[]);
fclose(fileID);    

fileID = fopen('unf.bin','r');
sunf = fread(fileID,'double'); sunf = reshape(sunf,npf,[],ncu);
fclose(fileID);    
e=permute(sunf,[1 3 2])-suh;
max(abs(e(:)))
fileID = fopen('uhg.bin','r');
suhg = fread(fileID,'double'); suhg = reshape(suhg,ngf,[],ncu);
fclose(fileID);    
fileID = fopen('xgf.bin','r');
sxgf = fread(fileID,'double'); sxgf = reshape(sxgf,ngf,[],nd);
fclose(fileID);    
fileID = fopen('nlgf.bin','r');
snlgf = fread(fileID,'double');snlgf = reshape(snlgf,ngf,[],nd);
fclose(fileID);    
fileID = fopen('jacf.bin','r'); 
sjacf = fread(fileID,'double'); sjacf = reshape(sjacf,ngf,[]);
fclose(fileID);    
fileID = fopen('fhg.bin','r'); 
sfhg = fread(fileID,'double');  sfhg = reshape(sfhg,ngf,[],ncq);
fclose(fileID);
fileID = fopen('Rqnf.bin','r'); 
sRqf = fread(fileID,'double');  sRqf = reshape(sRqf,npf,[],ncq);
fclose(fileID);
fileID = fopen('Rqf.bin','r'); 
sRq = fread(fileID,'double');  sRq = reshape(sRq,npe,ncq,[]);
fclose(fileID);    

master.shapft(:,:,2:end)
nf = mesh.nf;
f2e = reshape(permute(reshape(mesh.facecon,npf,2,nf),[2 1 3]),[2 npf*nf]);
qdg = zeros(npe,ncq,mesh.ne);
for i=1:nf
    ei = mesh.f(i,end-1:end);
    e1 = ei(1);
    e2 = ei(2);
    j = mesh.t2f(e1,:) == i;    
    u1 = UDG(mesh.perm(:,j),1:ncu,e1);
    ug = master.shapft(:,:,1)*u1;
    pb = mesh.dgnodes(mesh.perm(:,j),1:nd,e1);
    xg = master.shapft(:,:,1)*pb;
    for k = 1:nd-1
        dpg(:,:,k) = reshape(master.shapft(:,:,k+1)*pb,[npf nd]);
    end
    nmg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
    nmg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
    nmg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
    jac = sqrt(nmg(:,1).^2+nmg(:,2).^2+nmg(:,3).^2);
    nmg(:,1) = nmg(:,1)./jac;
    nmg(:,2) = nmg(:,2)./jac;
    nmg(:,3) = nmg(:,3)./jac;            
    eu = ug-reshape(suhg(:,i,:),[ngf ncu]);
    en = nmg-reshape(snlgf(:,i,:),[ngf nd]);
    ej = jac-reshape(sjacf(:,i,:),[ngf 1]);
    [max(abs(eu(:))) max(abs(en(:))) max(abs(ej(:)))]
    fg = [ug.*nmg(:,1) ug.*nmg(:,2) ug.*nmg(:,3)].*jac;
    Rqf = master.shapfg(:,:,1)*fg;
    ef = fg-reshape(sfhg(:,i,:),[npf ncq]);
    eq = Rqf-reshape(sRqf(:,i,:),[npf ncq]);
    [max(abs(ef(:))) max(abs(eq(:)))]
    if e2>0
        qdg = cuda_putfacenodes(qdg,Rqf,f2e,npf,ncq,npe,ncq,i,i,0);
        %pause
    end
    %[u1-suh(:,:,i)]
    %pause
%     pb = mesh.dgnodes(perm(:,j),1:2,e1);
%     dpg = reshape(dshapft*pb,[npf nd]);
%     jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
%     nl = [dpg(:,2)./jac,-dpg(:,1)./jac];
%     if e2>0
%         k = mesh.t2f(e2,:)==i;
%         u2 = UDG(perm(end:-1:1,k),5:6,e2);
%         qn =sum((u1+u2).*nl,2)./2;
%         Tot2 = Tot2 + (master.gwfc.*jac)'*(shapft*qn);
%     end
end
