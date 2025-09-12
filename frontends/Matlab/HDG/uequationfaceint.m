function [Ru, Rh, B, D, F, G, K, H] = uequationfaceint(master,pde,bf,dgnodes,UDG,UH,Ru,BD)
% FACEINTND compute face integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, nlg, jac] = facegeom(master.shapft,dgnodes,master.perm);

ne   = size(UDG,2);
nd   = master.nd;
npe  = master.npe;
npf  = master.npf;
ngf  = master.ngf;
nfe  = size(master.perm,2);
nc   = pde.nc;
ncu  = pde.ncu;
arg  = pde.arg;
time = pde.time;
bcm  = pde.bcm;        
bcs  = pde.bcs;
fbou   = str2func(pde.fbou);
fhat   = str2func(pde.fhat);

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

% DG solution at Gauss points
udgn = reshape(UDG(perm,:,:),[npf nfe*ne*nc]);
udgg = shapft*udgn;
udgg = reshape(udgg,[ngf*nfe*ne nc]);

uh  = reshape(UH,[npf nfe*ne*ncu]);
uhg = shapft*uh;
uhg = reshape(uhg,[ngf*nfe*ne ncu]);

% Compute numerical fluxes 
[FH, FH_udg, FH_uh] = fhat(nlg, pg, udgg, uhg, arg, time);     

% dataout = "/Users/ngoccuongnguyen/Exasim/build/dataout/";
% filename = dataout + "outuEquationElemFace_fh.bin";
% tmp = readbin(filename);
% fh = reshape(tmp, size(FH));
% max(abs(fh(:)-FH(:)))
% pause

tm     = find(bf<0);
nbf    = length(tm(:));
BH     = zeros(ngf*nbf,ncu);
BH_udg = zeros(ngf*nbf,ncu,nc);
BH_uh  = zeros(ngf*nbf,ncu,ncu);
an     = zeros(nbf*npf,1);
bn     = zeros(nbf*ngf,1);
j = 1;
for i = 1:length(bcm)    
    [I,J] = find(bf==-i);
    nfb = length(I);
    pgb  = zeros(ngf*nfb, size(pg,2));
    nlgb = zeros(ngf*nfb, nd);
    udgb = zeros(ngf*nfb, nc);
    uhgb = zeros(ngf*nfb, ncu);
    im = zeros(ngf*nfb,1);
    for k = 1:nfb
        i0 = (j-1)*ngf+1:j*ngf;
        i1 = (k-1)*ngf+1:k*ngf; 
        i2 = (J(k)-1)*nfe*ngf+(I(k)-1)*ngf+1:(J(k)-1)*nfe*ngf+I(k)*ngf; 
        bn(i0) = i2;
        im(i1)  = i0;
        pgb(i1,:)  = pg(i2,:);
        nlgb(i1,:) = nlg(i2,:);
        udgb(i1,:) = udgg(i2,:);
        uhgb(i1,:) = uhg(i2,:);
        i0 = (j-1)*npf+1:j*npf;
        i2 = (J(k)-1)*nfe*npf+(I(k)-1)*npf+1:(J(k)-1)*nfe*npf+I(k)*npf; 
        an(i0) = i2;
        j = j + 1;
    end  
    
    ib = bcm(i);
    bv = repmat(bcs(i,:),[ngf*nfb 1]);
    [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nlgb,pgb,udgb,uhgb,arg,time);                                                 
    BH(im,:) = fh;
    BH_udg(im,:,:) = dfh_dudg;
    BH_uh(im,:,:) = dfh_duh;                  
end

wrk = reshape(bsxfun(@times,FH,jac),[ngf nfe*ne*ncu]);
Run = reshape(shapfg*wrk,[npf nfe ne ncu]);

wrk = reshape(bsxfun(@times,FH_udg,jac),[ngf nfe*ne*ncu*nc]);
BDn = reshape(shapfgdotshapfc*wrk,[npf npf nfe ne ncu nc]);

wrk = reshape(bsxfun(@times,FH_uh,jac),[ngf nfe*ne*ncu*ncu]);
Fn = reshape(shapfgdotshapfc*wrk,[npf npf nfe ne ncu ncu]);
   
Rut = zeros(npe,1,ne,ncu);
BDt = zeros(npe,npe,1,ne,ncu,nc);
F = zeros(npe,npf,nfe,ne,ncu,ncu); 
for is=1:nfe  % Adding face contributions - vector dependencies avoided
    IP = perm(:,is);
    Rut(IP,1,:,:) = Rut(IP,1,:,:) + Run(:,is,:,:);
    BDt(IP,IP,1,:,:,:) = BDt(IP,IP,1,:,:,:) + BDn(:,:,is,:,:,:);
    F(IP,:,is,:,:,:) = F(IP,:,is,:,:,:) + Fn(:,:,is,:,:,:);
end
Ru = Ru-reshape(Rut,[npe ne ncu]);
BD = BD+reshape(BDt,[npe npe ne ncu nc]);
F = reshape(F,[npe npf*nfe ne ncu ncu]);

wrk = reshape(bsxfun(@times,BH,jac(bn)),[ngf nbf*ncu]);
Run = reshape(Run,[npf*nfe*ne ncu]);
Run(an,:) = reshape(shapfg*wrk,[npf*nbf ncu]);

Rh = -reshape(Run,[npf*nfe ne ncu]);

wrk = reshape(bsxfun(@times,BH_udg,jac(bn)),[ngf nbf*ncu*nc]);
Gt = reshape(BDn,[npf npf*nfe*ne ncu nc]);
Gt(:,an,:,:) = reshape(shapfgdotshapfc*wrk,[npf npf*nbf ncu nc]);
Gt = reshape(Gt,[npf npf nfe ne ncu nc]);

wrk = reshape(bsxfun(@times,BH_uh,jac(bn)),[ngf nbf*ncu*ncu]);
Ht = reshape(Fn,[npf npf*nfe*ne ncu ncu]);
Ht(:,an,:,:) = reshape(shapfgdotshapfc*wrk,[npf npf*nbf ncu ncu]);
Ht  = reshape(Ht,[npf npf nfe ne ncu ncu]);

Gt = permute(Gt,[1 3 2 4 5 6]); %Gt = npf npf nfe ne ncu nc
G   = zeros(npf,nfe,npe,ne,ncu,nc); 
H   = zeros(npf*nfe,npf*nfe,1,ne,ncu,ncu);
for is=1:nfe  
    IP = (is-1)*npf+1:is*npf;
    G(:,is,perm(:,is),:,:,:) = G(:,is,perm(:,is),:,:,:) + Gt(:,is,:,:,:,:);
    H(IP,IP,1,:,:,:) = H(IP,IP,1,:,:,:) + Ht(:,:,is,:,:,:);
end
G  = reshape(G,[npf*nfe npe ne ncu nc]);
H  = reshape(H,[npf*nfe npf*nfe ne ncu ncu]);

D = BD(:,:,:,:,1:ncu);
K = G(:,:,:,:,1:ncu);

if pde.ncq > 1
  B = BD(:,:,:,:,(ncu+1):nc);
  G = G(:,:,:,:,(ncu+1):nc);
else
  B = [];
  G = [];
end




