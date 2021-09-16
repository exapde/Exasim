
# Add Exasim to Julia search path
versiondir = cdir[1:ii[end]] * "/" * version * "/Julia";
push!(LOAD_PATH, cdir[1:ii[end]] * "/Installation");
push!(LOAD_PATH, versiondir * "/Gencode");
push!(LOAD_PATH, versiondir * "/Mesh");
push!(LOAD_PATH, versiondir * "/Preprocessing");
push!(LOAD_PATH, versiondir * "/Postprocessing");
push!(LOAD_PATH, versiondir * "/Utilities");

include(cdir[1:ii[end]] * "/Installation/cmakecompile.jl");

# Set Julia's PATH enviroment variable so that Exasim can call external programs
ENV["PATH"] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";
# Add more paths if neccesary
ENV["PATH"] =  ENV["PATH"] * ":/Applications/ParaView-5.8.1.app/Contents/MacOS";

print("==> Exasim " * version * " ...\n");
