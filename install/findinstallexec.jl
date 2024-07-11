function string2cmd(str::String)

ind = findall(" ", str);
n = length(ind);
cmdstr = Array{String,1}(undef,n+1);
cmdstr[1] = str[1:(ind[1][1]-1)];
for i = 2:n
    cmdstr[i] = str[(ind[i-1][1]+1):(ind[i][1]-1)];
end
i = n+1;
cmdstr[i] = str[(ind[i-1][1]+1):end];

return Cmd(cmdstr)

end

function findinstallexec(filename, appname, brew, searchopt)

print("Finding " * appname * "...\n");

dirname = "";
appinstall = 0;

status = Sys.which(filename);
if status != nothing
    dirname = filename;
    print("Exasim found " * dirname * "\n");
    return dirname, appinstall;
end
status = Sys.which("/usr/bin/" * filename);
if status != nothing
    dirname = "/usr/bin/" * filename;
    print("Exasim found " * dirname * "\n");
    return dirname, appinstall;
end
status = Sys.which("/usr/local/bin/" * filename);
if status != nothing
    dirname = "/usr/local/bin/" * filename;
    print("Exasim found " * dirname * "\n");
    return dirname, appinstall;
end
status = Sys.which("/opt/local/bin/" * filename);
if status != nothing
    dirname = "/opt/local/bin/" * filename;
    print("Exasim found " * dirname * "\n");
    return dirname, appinstall;
end

if searchopt==1
    if Sys.isapple()
        a = read(string2cmd("find /Applications -name " * filename * " -type f"), String);
        if length(a)>0
            ii = findlast("/MacOS/" * filename, a);
            dirname = a[1:ii[end]];
            print("Exasim found " * dirname * "\n");
            return dirname, appinstall;
        end
    end
end

print("Exasim could not find " * appname * " on your computer.\n");
if searchopt==10
    dirname = filename;
    print("CUDA Toolkit is not found on your system.\n");
    print("If you have Nividia GPUs on your system, please visit https://docs.nvidia.com/cuda/ to install it.\n");
    return dirname, appinstall;
else
    appinstall = 1;
    if Sys.isapple()
        print("Installing " * appname * " via brew.\n");
        if searchopt==1
            run(string2cmd(brew * " cask install " * appname));
        else
            run(string2cmd(brew * " install " * appname));
        end
    elseif Sys.isunix()
        print("Installing " * appname * " via apt.\n");
        run(string2cmd("sudo apt install " * appname));
    elseif Sys.windows()
        print("Installing " * appname * " via apt.\n");
        run(string2cmd("sudo apt install " * appname));
    end

    return dirname, appinstall;
end

return dirname, appinstall;

end
