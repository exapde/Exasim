function checkpreprocessing(path1, path2, mpiprocs)

disp("********************************");
match = checkpdebin(path1 + "/app.bin", path2 + "/app.bin");
if match==0
  disp("app.bin files DO NOT match.");
else
  disp("app.bin files match!!!");
end
disp("********************************");

disp(" ");

disp("********************************");
match = checkmasterbin(path1 + "/master.bin", path2 + "/master.bin");
if match==0
  disp("master.bin files DO NOT match.");
else
  disp("master.bin files match!!!");
end
disp("********************************");

disp(" ");

disp("********************************");

if mpiprocs==1
    match = checkmeshbin(path1 + "/mesh.bin", path2 + "/mesh.bin");
    if match==0
      disp("mesh.bin files DO NOT match.");
    else
      disp("mesh.bin files match!!!");
    end

    disp("********************************");

    disp(" ");

    disp("********************************");
    
    match = checksolbin(path1 + "/sol.bin", path2 + "/sol.bin");
    if match==0
      disp("sol.bin files DO NOT match.");
    else
      disp("sol.bin files match!!!");
    end    
    disp("********************************");
else
    for i = 1:mpiprocs
        match = checkmeshbin(path1 + "/mesh" + num2str(i) + ".bin", path2 + "/mesh" + num2str(i) + ".bin");
        if match==0
          disp("mesh" + num2str(i) + ".bin files DO NOT match.");
        else
          disp("mesh" + num2str(i) + ".bin files match!!!");
        end

        disp("********************************");

        disp(" ");

        disp("********************************");
        
        match = checksolbin(path1 + "/sol" + num2str(i) + ".bin", path2 + "/sol" + num2str(i) + ".bin");
        if match==0
          disp("sol" + num2str(i) + ".bin files DO NOT match.");
        else
          disp("sol" + num2str(i) + ".bin files match!!!");
        end          
        disp("********************************");

        if i<mpiprocs
            disp(" ");

            disp("********************************");        
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
