function partition(t2f,ne,np,metis)

# Generate a temporary file to be used in METIS
print("Writing input files for METIS...\n");
open("temp.txt", "w") do io
     Main.DelimitedFiles.writedlm(io, ne)
end
open("temp.txt", "a") do io
     Main.DelimitedFiles.writedlm(io, t2f')
end

metisstatus0 = Sys.which(metis);
metisstatus1 = Sys.which("mpmetis");
metisstatus2 = Sys.which("/usr/bin/mpmetis");
metisstatus3 = Sys.which("/usr/local/bin/mpmetis");
metisstatus4 = Sys.which("/opt/local/bin/mpmetis");

if metisstatus0 != nothing
elseif metisstatus1 != nothing
    metis = "mpmetis"
elseif metisstatus2 != nothing
    metis = "/usr/bin/mpmetis";
elseif metisstatus3 != nothing
    metis = "/usr/local/bin/mpmetis";
elseif metisstatus4 != nothing
    metis = "/opt/local/bin/mpmetis";
else
    error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Metis. Please see the documentation to install it. After installation, please set its path to app.metis");
end

print("Calling METIS and reading output files...\n");
# call mpmetis
nps = string(np);
str = `$metis $"temp.txt" $nps`;
run(str);

# get mesh partitioning data
estr = string("temp.txt.epart.", nps);
#epart = Main.DelimitedFiles.readdlm(estr, ',', Int);
epart = Main.DelimitedFiles.readdlm(estr, Int);

# get node partitioning data
nstr = string("temp.txt.npart.", nps);
#npart = Main.DelimitedFiles.readdlm(nstr, ',', Int);
npart = Main.DelimitedFiles.readdlm(nstr, Int);

# remove files
rm("temp.txt");
rm(estr);
rm(nstr);

return epart, npart

end

# #if size(t2f,2) == ne; t2f = t2f'; end
#
# # current directory
# cdir = pwd();
# ii = findlast("Exasim", cdir);
#
# # move to directory that contains metis programs
# if Sys.isapple()
#     #cd([cdir "/metis/mac"]);
#     metisdir = string(cdir[1:ii[end]],"/metis/mac");
# elseif Sys.isunix()
#     #cd([cdir "/metis/linux"]);
#     metisdir = string(cdir[1:ii[end]],"/metis/linux");
# elseif Sys.iswindows()
#     #cd([cdir "/metis/linux"]);
#     metisdir = string(cdir[1:ii[end]],"/metis/windows");
# end
# cd(metisdir);
#
# # Generate a temporary file to be used in METIS
# display("Writing input files for METIS...");
# open("temp.txt", "w") do io
#      Main.DelimitedFiles.writedlm(io, ne)
# end
# open("temp.txt", "a") do io
#      Main.DelimitedFiles.writedlm(io, t2f')
# end
#
# #open(io -> writedlm(io, ne, ','), "temp.txt", "w") # write
# #open(io -> writedlm(io,t2f, ','), "temp.txt", "a") # append
# # writedlm('temp.txt', ne, 'delimiter', ' ','precision',10);
# # writedlm('temp.txt', t2f, '-append', 'delimiter', ' ','precision',10);
#
# display("Calling METIS and reading output files...");
# # call mpmetis
# nps = string(np);
# str = `$"./mpmetis" $"temp.txt" $nps`;
# run(str);
#
# # get mesh partitioning data
# estr = string("temp.txt.epart.", nps);
# epart = Main.DelimitedFiles.readdlm(estr, ',', Int);
#
# # get node partitioning data
# nstr = string("temp.txt.npart.", nps);
# npart = Main.DelimitedFiles.readdlm(nstr, ',', Int);
#
# # remove files
# rm("temp.txt");
# rm(estr);
# rm(nstr);
#
# # move back to current directory
# cd(cdir);
#
# return epart, npart
#
# end
