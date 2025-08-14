function checkuhat(dmd,UH)

% number of subdomains
nproc = length(dmd);

for i = 1:nproc
  for j = 1:size(dmd{i}.interfaces,1)
    n1 = i;
    f1 = dmd{i}.interfaces(j,2);
    n2 = dmd{i}.interfaces(j,9);    
    f2 = dmd{i}.interfaces(j,10);    
    eUH = UH{n1}(:,:,f1) - UH{n2}(:,:,f2);
    if max(abs(eUH(:))) > 1e-12
      disp(dmd{i}.interfaces(j,:));
      error("UH does not match");
    end
  end
end

disp("UH matches on the interfaces between subdomains!");


