function checkelementpartition

load dmdelem.mat

f = ["t2t", "elem2cpu"];

 for i = 1:length(f)
  s = f(i);  
  a = readbin(s + ".bin",'int32');
  t = eval(s);
  if numel(a) == numel(t)
    fprintf("Check " + s + ":  (%g, %g)\n", max(abs(a(:)-t(:))), max(abs(1+a(:)-t(:))));
  else
    fprintf(s + ": sizes do not match (%d, %d)\n", numel(a), numel(t));
  end
 end

f = ["elempart", "elempartpts", "elem2cpu", "nbsd", "elemrecv", "elemsend", "elemrecvpts", "elemsendpts"];

for j = 1:length(dmd)
 fprintf("Subdomain: %d\n", j); 
 for i = 1:length(f)
  s = f(i);  
  a = readbin(s + num2str(j-1) + ".bin",'int32');
  t = eval("dmd{j}." + s);
  if numel(a) == numel(t)
    fprintf("Check " + s + ":  (%g, %g)\n", max(abs(a(:)-t(:))), max(abs(1+a(:)-t(:))));
    %[t reshape(a, size(t))]
  else
    fprintf(s + ": sizes do not match (%d, %d)\n", numel(a), numel(t));
  end
 end
end

% load tmp.mat
% % save tmp.mat dmd eblks fblks rowe2f1 cole2f1 ent2ind1 rowe2f2 cole2f2 ent2ind2 cgelcon rowent2elem colent2elem cgent2dgent
%     
% f = ["elemcon" "facecon" "bf" "f2t" "facepartpts" "facepartbnd" "eblks" "fblks" ];
% f = [f "facecon1" "facecon2" "rowe2f1" "cole2f1" "ent2ind1" "rowe2f2" "cole2f2"];
% f = [f "ent2ind2" "cgelcon" "rowent2elem" "colent2elem" "cgent2dgent"]; %  
% 
% i = 1;
% elemcon = dmd{i}.elemcon;
% facecon = dmd{i}.facecon;
% bf = dmd{i}.bf;
% f2t = dmd{i}.f2t;
% facepartpts = dmd{i}.facepartpts;
% facepartbnd = dmd{i}.facepartbnd;
% facecon1 = reshape(facecon(:,1,:),[size(facecon,1) size(facecon,3)]);
% facecon2 = reshape(facecon(:,2,:),[size(facecon,1) size(facecon,3)]);      
% ind = [];        
% for ii = 1:size(fblks,2)
%     if fblks(3,ii)>0
%         ind = [ind fblks(1,ii):fblks(2,ii)];
%     end
% end            
% facecon2(:,ind)=[];        
% facecon = permute(facecon, [2 1 3]);
% 
% for i = 1:length(f)
%   s = f(i);  
%   a = readbin(s + ".bin",'int32');
%   t = eval(s);
%   if numel(a) == numel(t)
%     fprintf("Check " + s + ":  (%g, %g)\n", max(abs(a(:)-t(:))), max(abs(1+a(:)-t(:))));
%   else
%     fprintf(s + ": sizes do not match (%d, %d)\n", numel(a), numel(t));
%   end
% end
%  
