function meshplot3D(p,t)

% clf;

n = [1, 0, 0];
p0 =[0, 0, 0];
D  = n*p0';
vis1 = (1- (p(1,:)>1e-3).*(p(2,:)>1e-3).*(p(3,:)>1e-3));
vis2 = (p(1,:).^2 + p(2,:).^2 + p(3,:).^2) < 636^2 - 1;
vis2 = min(vis1+vis2,1);

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
fi = fi(:,ind(prod(vis1(fi),1)>0));

for i = 1:nf-ib+1
    fb(:,i) = t(face(:,f2t(2,i+ib-1)),f2t(1,i+ib-1));
end
fb = fb(:,ind(prod(vis2(fb),1)>0));

figure
pars = {'facecolor',[0.8,1.0,0.8],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',fi','vertices',p',pars{:}); hold on;

pars = {'facecolor',[0.8,1.0,1.0],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',fb','vertices',p',pars{:}); hold on;
axis equal;