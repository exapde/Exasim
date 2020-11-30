using Gencode

function findexec(filename, version)

status = Sys.which(filename);
if status != nothing
    return filename;
end
status = Sys.which("/usr/bin/" * filename);
if status != nothing
    filename = "/usr/bin/" * filename;
    return filename;
end
status = Sys.which("/usr/local/bin/" * filename);
if status != nothing
    filename = "/usr/local/bin/" * filename;
    return filename;
end
status = Sys.which("/opt/local/bin/" * filename);
if status != nothing
    filename = "/opt/local/bin/" * filename;
    return filename;
end
if status == nothing
    print("Exasim could not find " * filename * " in /usr/bin, /usr/local/bin, /opt/local/bin\n");
    if Sys.isapple()
        print("Exasim tries to find " * filename * " in /Applications. It may take a while.\n");
        a = read(Gencode.string2cmd("find /Applications -name " * filename * " -type f"), String);
        if length(a)>0
            ii = findlast("/MacOS/" * filename, a);
            newfilename = a[1:ii[end]];
            print("Exasim found " * filename * " at " * newfilename * "\n");
            print("Please open initializepde.jl in the folder Exasim/" * version * "/Julia/Preprocessing\n");
            qstr1 = """ "$filename" """;
            qstr2 = """ "$newfilename" """;
            print("Replace pde." * filename * " =" * qstr1 * " with pde." * filename * " =" * qstr2 * "\n");
            print("Doing so will prevent Exasim from searching " * filename * " again to save time.\n");
            print("Read the above instructions carefully and press any key to continue...\n");
            sleep(20);
            filename = newfilename;
            return filename;
        end
    # elseif Sys.isunix()
    # elseif Sys.windows()
    end
    mystr = "Exasim could not find " * filename * " on your system.\n";
    mystr = mystr * "Please see the documentation to install " * filename * ".\n"
    mystr = mystr * "After installing " * filename * ", open initializepde.jl in the folder Exasim/" * version * "/Julia/Preprocessing\n";
    qstr1 = """ "$filename" """;
    newfilename = "/path/to/executable/" * filename;
    qstr2 = """ "$newfilename" """;
    mystr = mystr * "and replace pde." * filename * " =" * qstr1 * "with pde." * filename * " =" * qstr2 * "\n";
    print(mystr);
    error("Exit Exasim!")
end

end
