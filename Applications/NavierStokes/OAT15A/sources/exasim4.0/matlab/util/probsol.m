cd('app');

if ismac
    mpirun = '!/opt/local/bin/mpirun-openmpi-mp';
else
    mpirun = '!mpirun';
end

if (gpu==0)
    if mpiprocs>1
        str = [mpirun ' -np ' num2str(mpiprocs) ' ./mpi'];    
    else
        str = '!./cpu';
    end
else
    if mpiprocs>1
        str = [mpirun ' -np ' num2str(mpiprocs) ' ./gpumpi'];          
    else
        str = '!./gpu';
    end
end
str = [str app.appname ' ../data/' app.appname ' ../data/' app.appname 'out'];
eval(char(str));    

cd('..');
