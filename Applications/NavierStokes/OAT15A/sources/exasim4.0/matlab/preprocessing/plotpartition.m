function plotpartition(p,t,f,dmd)

bcol = [1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
       ];
      
np = length(dmd);
t = t';

% plot nonoverlapping subdomains
figure(1); clf;
hold on;        
for i=1:np
    ind = dmd{i}.elempart(1:sum(dmd{i}.elempartpts(1:2)));
    ti = t(:,ind);
    simpplot(p,ti',[],bcol(i,:));                       
end
for it=1:size(t,2)
    pmid=mean(p(t(:,it),:),1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
    text(pmid(1),pmid(2),num2str(it),txtpars{:});
end
% for it=1:size(f,1)
%     pmid=mean(p(f(it,1:2),:),1);
%     txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%     text(pmid(1),pmid(2),num2str(it),txtpars{:});
% end

hold off;
axis equal;      
axis tight;
axis on;  

% plot overlapping subdomains
for n=1:np               
    figure(n+1); clf;    
    hold on;        
    simpplot(p,t',[],'none');   
    
    tn = t(:,dmd{n}.elempart(1:sum(dmd{n}.elempartpts(1:2))));
    simpplot(p,tn',[],bcol(n,:));                           
        
    tn = t(:,dmd{n}.elempart(sum(dmd{n}.elempartpts(1:2))+1:end));
    simpplot(p,tn',[],[0 0 0]);                           
    
    nelem = length(dmd{n}.elempart);
    for it=1:nelem
        pmid=mean(p(t(:,dmd{n}.elempart(it)),:),1);
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
        text(pmid(1),pmid(2),num2str(dmd{n}.elempart(it)),txtpars{:});
    end
    
%     for it=1:length(dmd{n}.facepart)
%         pmid=mean(p(f(dmd{n}.facepart(it),1:2),:),1);
%         txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%         text(pmid(1),pmid(2),num2str(dmd{n}.facepart(it)),txtpars{:});
%     end
    
    hold off;
    axis equal;      
    axis tight;
    axis on;  
end

figure(4); clf;
hold on;        
for i=1:np
    ind = dmd{i}.elempart(1:sum(dmd{i}.elempartpts(1:2)));
    ti = t(:,ind);
    simpplot(p,ti',[],bcol(i,:));                       
end
for it=1:size(f,1)
    pmid=mean(p(f(it,1:2),:),1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
    text(pmid(1),pmid(2),num2str(it),txtpars{:});
end

% plot overlapping subdomains
for n=1:np               
    figure(n+4); clf;    
    hold on;        
    simpplot(p,t',[],'none');   
    
    tn = t(:,dmd{n}.elempart(1:sum(dmd{n}.elempartpts(1:2))));
    simpplot(p,tn',[],bcol(n,:));                           
        
    tn = t(:,dmd{n}.elempart(sum(dmd{n}.elempartpts(1:2))+1:end));
    simpplot(p,tn',[],[0 0 0]);                           
        
    for it=1:length(dmd{n}.facepart)
        pmid=mean(p(f(dmd{n}.facepart(it),1:2),:),1);
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
        text(pmid(1),pmid(2),num2str(dmd{n}.facepart(it)),txtpars{:});
    end
    
    hold off;
    axis equal;      
    axis tight;
    axis on;  
end

% % plot overlapping subdomains
% for n=1:np               
%     figure(n+1); clf;    
%     hold on;        
%     simpplot(p,t',[],'none');   
%     
%     tn = t(:,extintelem{n}(1:sum(extintelempts(n,1:2)))+1);
%     simpplot(p,tn',[],bcol(n,:));                           
%         
%     tn = t(:,extintelem{n}(sum(extintelempts(n,1:2))+1:end)+1);
%     simpplot(p,tn',[],[0 0 0]);                           
% 
% %     for it=1:size(fedg,1)
% %         pmid=mean(p(fedg(it,1:2),:),1);
% %         if ismember(it,extintent{n}+1)
% %             txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-n,:)};
% %         else
% %             txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1 1 1]};   
% %         end
% %         text(pmid(1),pmid(2),num2str(it),txtpars{:});
% %     end
% 
%     hold off;
%     axis equal;      
%     axis tight;
%     axis on;  
%     
% %     intent = find(ent2cpu==n-1);
% %     extintent{n}'+1
% %     intent'
% %     pause    
% end
% 
