function print_results(a,type)
if type==1
a=round(a/0.01)*0.01;
c=[];
for i=1:size(a,1)
    b=' ';
    for j=1:size(a,2)
        if j==1
            b=[b   num2str(a(i,j)) '  &  '];
        elseif j < size(a,2)
            b=[b num2str(a(i,j)) '  &  '];
        else
            b=[b num2str(a(i,j)) '  \\  '];
        end
    end
    b=strrep(b, 'e', 'E');
    disp(b)
end
else
ind=3:2:size(a,2);    
a(:,ind)=round(a(:,ind)/0.01)*0.01;
ind=2:2:size(a,2)-1;    
%ind=1:1:size(a,2);    
a(:,ind)=a(:,ind)/1e10;
c=[];
for i=1:size(a,1)
    for j=1:length(ind)
        for k=0:30
            if a(i,ind(j))*10^k>100
                break;
            end
        end
        a(i,ind(j))=round(a(i,ind(j))*(10^(k)))/10^(k);
    end
end
for i=1:size(a,1)
    b=' & ';
    for j=1:size(a,2)
        if j==1
            b=[b   num2str(a(i,j)) '  &  '];
        elseif j < size(a,2)
            b=[b num2str(a(i,j)) '  &  '];
        else
            b=[b num2str(a(i,j)) '  \\  '];
        end
    end
    %b=strrep(b, 'e', 'E');
    b=strrep(b, 'e-0', '\mbox{e-}');
    b=strrep(b, 'e-1', '\mbox{e-}');
    b=strrep(b, 'e-2', '\mbox{e-}1');
    b=strrep(b, '&  0', '&  --');
    b=strrep(b, '--.', ' 0.');
    disp(b)
    if mod(i,5)==0
        disp('\hline');
    end
end    
end
