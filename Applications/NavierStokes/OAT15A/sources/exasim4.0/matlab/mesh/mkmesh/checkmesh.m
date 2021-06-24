function checkmesh(f, t2f, facecon, dgnodes, perm)

ne = size(t2f,1);
nf = size(f,1);

figure(1); clf;
hold on;
for i = 1:ne
    x = dgnodes(perm(:),:,i);
    plot(x(:,1),x(:,2),'-b');
    pmid=mean(x,1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
    text(pmid(1),pmid(2),num2str(i),txtpars{:});    
end

% figure(2); clf;
% hold on;
x = dgnodes(:,1,:);
y = dgnodes(:,2,:);
for i = 1:nf
    ei = f(i,end-1:end);
    f1 = find(t2f(ei(1),:)==i);
    x1 = dgnodes(perm(:,f1),:,ei(1));
    if ei(2)>0
        f2 = find(t2f(ei(2),:)==i);
        x2 = dgnodes(perm(:,f2),:,ei(2));        
    end    
    plot(x1(:,1),x1(:,2),'-r');        
    pmid=mean(x1,1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[0,1,0]};
    text(pmid(1),pmid(2),num2str(i),txtpars{:});    
    
    fi = facecon(:,:,i);    
    [x(fi(:,1)) y(fi(:,1)) x(fi(:,2)) y(fi(:,2))]
end
