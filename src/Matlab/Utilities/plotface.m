function plotface(p, t, f)

pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
patch('faces',t,'vertices',p,pars{:});  

for it=1:size(f,1)
    pmid=mean(p(f(it,1:2),:),1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
    text(pmid(1),pmid(2),num2str(it),txtpars{:});
end        

