function filename = findexec(filename, version)

[status,~] = system("which " + filename);
if status==0
    return;
end
[status,~] = system("which /usr/bin/" + filename);
if status==0
    filename = "/usr/bin/" + filename;
    return;
end
[status,~] = system("which /usr/local/bin/" + filename);
if status==0
    filename = "/usr/local/bin/" + filename;
    return;
end
[status,~] = system("which /opt/local/bin/" + filename);
if status==0
    filename = "/opt/local/bin/" + filename;
    return;
end
if status~=0
    disp("Exasim could not find " + filename + " in /usr/bin, /usr/local/bin, /opt/local/bin");
    if ismac
        disp("Exasim tries to find " + filename + " in /Applications. It may take a while.");
        [~, b] = system("find /Applications -name " + filename + " -type f");
        if ~isempty(b)
            ii = strfind(b, "/MacOS/" + filename);
            newfilename = b(1:ii(1)) + "MacOS/" + filename;
            disp("Exasim found " + filename + " at " + newfilename);
            disp("Please open initializepde.m in the folder Exasim/" + version + "/Matlab/Preprocessing");
            disp("Replace pde." + filename + " = " + """" + filename + """" + " with pde." + filename + " = " + """" + newfilename + """");
            disp("Doing so will prevent Exasim from searching " + filename + " again");
            disp("Read the above instructions carefully and press any key to continue...");
            open("initializepde.m");
            pause();
            filename = newfilename;
            return;
        end
%     elseif isunix
%     elseif ispc
    end
    disp("Exasim could not find " + filename + " on your system.");
    disp("Please see the documentation to install " + filename + ".");
    disp("After installing " + filename + ", open initializepde.m in the folder Exasim/Version0.1/Matlab/Preprocessing");
    disp("and replace pde." + filename + " = " + """" + filename + """" + " with pde." + filename + " = " + """" + "/path/to/executable/" + filename + """");
    error("Exit Exasim!");
end

end
