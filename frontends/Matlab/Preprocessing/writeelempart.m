function writeelempart(dmd)

for j = 1:length(dmd)
     fname = sprintf("elempartpts_np%d.bin", j-1);
     fname= fullfile("elempart", fname);
     writebin(fname, dmd{j}.elempartpts)     

     fname = sprintf("elempart_np%d.bin", j-1);
     fname= fullfile("elempart", fname);
     writebin(fname, dmd{j}.elempart)
end

end

