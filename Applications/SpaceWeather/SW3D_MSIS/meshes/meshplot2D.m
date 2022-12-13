function meshplot2D(p,t)

figure; clf;
pars = {'facecolor',[.8,1,.8],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',t','vertices',p',pars{:});
axis equal
