function [dirname,appinstall] = findinstallexec(filename, appname, brew, searchopt)

disp("Finding " + appname + "...");

dirname = "";
appinstall = 0;

[status,~] = system("which " + filename);
if status==0    
    dirname = filename;    
    disp("Exasim found " + dirname);
    return;
end
[status,~] = system("which /usr/bin/" + filename);
if status==0
    dirname = "/usr/bin/" + filename;    
    disp("Exasim found " + dirname);
    return;
end
[status,~] = system("which /usr/local/bin/" + filename);
if status==0
    dirname = "/usr/local/bin/" + filename;
    disp("Exasim found " + dirname);
    return;
end
[status,~] = system("which /opt/local/bin/" + filename);
if status==0
    dirname = "/opt/local/bin/" + filename;
    disp("Exasim found " + dirname);
    return;
end

if searchopt==1
    if ismac                    
        [~, b] = system("find /Applications -name " + filename + " -type f");
        if ~isempty(b)
            ii = strfind(b, "/MacOS/" + filename);
            dirname = b(1:ii(1)) + "MacOS/" + filename;           
            disp("Exasim found " + dirname);
            return;
        end        
    end        
end

disp("Exasim could not find " + appname + " on your computer.");
if searchopt==10
    dirname = filename;
    disp("CUDA Toolkit is not found on your system.");
    disp("If you have Nividia GPUs on your system, please visit https://docs.nvidia.com/cuda/ to install it.");
    return;
else
    appinstall = 1;
    if ismac
        disp("Installing " + appname + " via brew.");
        if searchopt==1
            system(brew + " cask install " + appname);     
        else
            system(brew + " install " + appname);     
        end
    elseif isunix
        disp("Installing " + appname + " via apt.");
        system("sudo apt install " + appname);
    elseif ispc
        disp("Installing " + appname + " via apt.");
        system("sudo apt install " + appname);    
    end
end

end

