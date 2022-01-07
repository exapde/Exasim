n = size(a,1);
[c,ia,ic] = unique(a,'rows');
ia = ia';
ic = ic';
ii = setdiff(1:n,ia);

nc = length(ia);
ib = zeros(1,n);
ib(ia) = 1:nc;

while (1)
    size(ii)
    [b,ij] = unique(a(ii,:),'rows');

    nb = size(b,1);
    ind = zeros(1,nb);
    for i = 1:nc
        if max(abs(b(1,:)-c(i,:)))<1e-10
            ind(1) = i;
            break;
        end
    end
    m = ind(1);
    for k = 2:nb
        for q = m:nc
            if max(abs(b(k,:)-c(q,:)))<1e-10  %b(k)==c(q)
                ind(k) = q;    
                m = q;
                break;
            end
        end            
    end
    %a(ii(ij(k)),:)==b(k,:)
    ib(ii(ij)) = ind;
    
    i1 = setdiff(1:length(ii),ij);
    if isempty(i1)
        break;
    else
        ii = ii(i1);
    end
end


%i1 = setdiff(1:length(ii),ij);

%a(ii(i1),:)
%a(ii(ij),:) = b;

% for k = 1:length(ij)
%     %a(ii(ij(k)),:)==b(k,:)
%     ib(ii(ij(k))) = ind(k);
% end


%[~,tm]=ismember(b,c,'rows');

% [b,ij] = sortrows(a(ii,:));
% 
% for i = 1:nc
%     if max(abs(b(1,:)-c(i,:)))<1e-10
%         m = i;
%         break;
%     end
% end
% ib(ii(ij(1))) = m;
% for k = 2:length(ii)
%     if max(abs(b(k,:)-b(k-1,:)))<1e-10 
%         ib(ii(ij(k))) = ib(ii(ij(k-1)));
%     else
%         for q = m:nc
%             if max(abs(b(k,:)-c(q,:)))<1e-10  %b(k)==c(q)
%                 ib(ii(ij(k))) = q;    
%                 m = q;
%                 break;
%             end
%         end        
%     end
% end
% 







% 
% %a = [1 2 3 5 5 5 6 7 7 7 8 8 9 10 10 10 10 11];
% % n = length(a);
% %a = a(randperm(n));
% % [c,ia,ic] = unique(a);
% % ia = ia'; ic=ic';
% % id = 0*ic;
% % id(ia) = 1:length(c);
% % ib = setdiff(1:n,ia);
% %id(ib) = ;
% 
% [c,ia,ic] = unique(a);
% ia=ia';ic=ic';
% ii = setdiff(1:n,ia);
% 
% ib = 0*ic;
% ib(ia) = 1:length(ia);
% 
% [b,ij] = sort(a(ii));
% %m = find(c==b(1));
% for i = 1:length(c)
%     if b(1)==c(i)
%         m = i;
%         break;
%     end
% end
% ib(ii(ij(1))) = m;
% for k = 2:length(ii)
%     if b(k)==b(k-1)
%         ib(ii(ij(k))) = ib(ii(ij(k-1)));
%     else
%         for q = m:length(c)
%             if b(k)==c(q)
%                 ib(ii(ij(k))) = q;    
%                 m = q;
%                 break;
%             end
%         end
%     end
% end
% 
% % b = unique(a(ii));
% % ij = 0*b;
% % for k = 1:length(b)
% %     ij(k) = find(c==b(k));    
% % end
% % for k = 1:length(ii)
% %     m = find(b==a(ii(k)));
% %     ib(ii(k)) = ij(m);    
% % end
% 
% % for k = 1:length(ii)
% %     ib(ii(k)) = find(c==a(ii(k)));    
% % end
% 
% % [b,i] = sort(a);
% % 
% % c = b;
% % ib = ones(1,n);
% % ic = ones(1,n);
% % m = 1;
% % for k = 2:length(b)
% %     if b(k)~=b(k-1)
% %         m = m + 1;
% %         c(m) = b(k);
% %         ib(m) = k;            
% %     end
% %     ic(k) = m;
% % end
% % c = c(1:m);
% % ib = ib(1:m);
% % 
% % [c1,i1,i2]=unique(b);
% % i1 = i1'; i2 = i2';
% % [i1; ib]
% % [i2; ic]
% % 
% % id = ic(i);
% % ia = i(ib);
% % 
% % [c1,i1,i2]=unique(a);
% % i1 = i1'; i2 = i2';
% % [i1; ia]
% % [i2; id]

% a = [1 2 3 5 5 5 6 7 7 7 8 8 9 10 10 10 10 11];
% n = length(a);
% a = a(randperm(n));

% [b,i] = sort(a);
% 
% d = b(2:end)-b(1:end-1);
% i1 = find(d==0);
% i2 = i1 + 1;
% i3 = unique([i1 i2]);  % indices for nonunique members of b
% i4 = setdiff(1:n, i3); % indices for unique members of b
% 
% c = b(i3);
% e = c(2:end) - c(1:end-1);
% i5 = find(e~=0);
% i5 = [1 i5+1]; % indices 
% 
% i6 = sort([i4 i3(i5)]);
% 
% c = b(i6);
% ia = i(i6);
% 
% 
% 
% 
% 
