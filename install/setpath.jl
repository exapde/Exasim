
# Add Exasim to Julia search path
src = "frontends"; 
srcdir = cdir[1:ii[end]] * "/" * src * "/Julia";
push!(LOAD_PATH, cdir[1:ii[end]] * "/install");
push!(LOAD_PATH, srcdir * "/Gencode");
push!(LOAD_PATH, srcdir * "/Mesh");
push!(LOAD_PATH, srcdir * "/Preprocessing");
push!(LOAD_PATH, srcdir * "/Postprocessing");
push!(LOAD_PATH, srcdir * "/Utilities");

include(cdir[1:ii[end]] * "/install/cmakecompile.jl");

# Set Julia's PATH enviroment variable so that Exasim can call external programs
ENV["PATH"] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";
# Add more paths if neccesary
ENV["PATH"] =  ENV["PATH"] * ":/Applications/ParaView-5.8.1.app/Contents/MacOS";

print("==> Exasim ...\n");
