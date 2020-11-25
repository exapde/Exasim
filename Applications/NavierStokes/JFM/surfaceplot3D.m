function surfaceplot3D(p,t,t_flag)

clf;

ind = ismember(t_flag,[2,3,14,25,33,34,35,36]);

fp = t(:,~ind);

pars = {'facecolor',[0.8,1.0,0.8],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',fp','vertices',p',pars{:}); hold on;

axis equal; grid on;