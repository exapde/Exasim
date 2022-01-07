function ptplot(p,t)

pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1.5,'FaceAlpha',1.0,'EdgeAlpha',1};
face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]]';
figure(2); clf;
hold on;
for i = 1:size(t,1)       
    f = reshape(t(i,face),[4 6]);     
    patch('faces',f','vertices',p,pars{:});
end
