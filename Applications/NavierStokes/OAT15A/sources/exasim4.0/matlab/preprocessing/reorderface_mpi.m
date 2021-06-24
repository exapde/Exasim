function [f,t2f,facecon,mf,bcn] = reorderface_mpi(f,t2f,facecon,bcm,perm)

disp('run reorderface...');  

% npe = max(perm(:));
% ind = f(:,end)<0; 
% a = -unique(f(ind,end));
% a = sort(a);
% b = (1:length(bcm));
% if length(a)~=length(b)
%     error('Something wrong. Check bcm or f');
% end
% if max(abs(a(:)-b(:)))~=0
%     error('Something wrong. Check bcm or f');
% end

%bcm = bcm(a);      % reorder bcm
bcn = unique(bcm); % remove duplications
nbc = length(bcn); % number of boundary conditions
inb = cell(nbc,1);
for i = 1:nbc
    inb{i} = find(bcm==bcn(i));
end

nf = size(f,1);
mp = zeros(nf,1);
mf = zeros(2,1);

% interior faces
ind = find(f(:,end)>0); 
n = length(ind);

% reorder the list
mp(1:length(ind))=ind;
mf(2) = n;

% boundary faces
for i = 1:nbc % for each boundary condition
    for j = 1:length(inb{i}) % for each boundary with boundary condition i 
        ind = find(f(:,end)==-inb{i}(j));
        if ~isempty(ind)
        mp((n+1):(n+length(ind)))=ind;
        n = n + length(ind);
        mf(i+2) = n;
        end
    end
end


f = f(mp,:);
facecon = facecon(:,:,mp);

[nt,nfv] = size(t2f);
a = zeros(nf,1);
a(mp) = (1:nf)';
t2f = reshape(t2f,nfv*nt,1);
ii = find(t2f > 0);
t2f(ii) = a(t2f(ii));
ii = find(t2f < 0);
t2f(ii) = -a(-t2f(ii));
t2f = reshape(t2f,nt,nfv);

% % break faces into groups so that all faces in one group do not share the same dofs
% fcolor = [];
% newmf = 0;
% newbcn = [];
% 
% for i = 1:length(mf)-1 % for each block of faces
%     % break the block into smaller blocks
%     %[acolor,ncolor] = mkcolorface(facecon,npe,mf(i)+1,mf(i+1));
%     [acolor,ncolor] = colorface(f,t2f,perm,(mf(i)+1):1:mf(i+1));
%     fcolor = [fcolor acolor];
%     newmf = [newmf ncolor(2:end)+mf(i)];
%     if i==1
%         newbcn = [newbcn zeros(1,length(ncolor(2:end)))];
%     else
%         newbcn = [newbcn bcn(i-1)*ones(1,length(ncolor(2:end)))];
%     end
% end
% if (length(unique(fcolor))~=nf) || (newmf(end)~=nf)
%     error('something wrong');
% end
% mf = newmf;
% bcn = newbcn(2:end);
% 
% mp = fcolor;
% f = f(mp,:);
% facecon = facecon(:,:,mp);
% [nt,nfv] = size(t2f);
% a = zeros(nf,1);
% a(mp) = (1:nf)';
% t2f = reshape(t2f,nfv*nt,1);
% ii = find(t2f > 0);
% t2f(ii) = a(t2f(ii));
% ii = find(t2f < 0);
% t2f(ii) = -a(-t2f(ii));
% t2f = reshape(t2f,nt,nfv);
