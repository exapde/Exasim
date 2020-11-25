function meshplot3D(p,t)

clf;

n = [0 0 1];
p0 =[0, 0, -100];
D  = n*p0';
vis = (n*p - D > 0.0);
%vis = sum(p.^2,1) > 0;

face = getelemface(3,1);    
f2t = mkf2e(t,1,3);
nf = size(f2t,2);
ib = find(f2t(3,:) == 0,1);

fi = zeros(4,ib-1);
fb = zeros(4,nf-ib+1);

for i = 1:ib-1
    fi(:,i) = t(face(:,f2t(2,i)),f2t(1,i));
end
ind = 1:ib-1;
fi = fi(:,ind(prod(vis(fi),1)>0));

for i = 1:nf-ib+1
    fb(:,i) = t(face(:,f2t(2,i+ib-1)),f2t(1,i+ib-1));
end
ind = 1:nf-ib+1;
fb = fb(:,ind(prod(vis(fb),1)>0));

pars = {'facecolor',[0.8,1.0,0.8],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',fi','vertices',p',pars{:}); hold on;

pars = {'facecolor',[0.8,1.0,1.0],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',fb','vertices',p',pars{:}); hold on;
axis equal;