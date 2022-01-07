% create app folder
mkdir('app');
cd('app');

% reshape symbolic arrays
flux{d} = reshape(flux{d},ncu,d);
source{d} = reshape(source{d},ncu,1);
tdfunc{d} = reshape(tdfunc{d},ncu,1);
if ~isempty(ubou{d})
    ubou{d} = reshape(ubou{d},ncu,[]);
end
if ~isempty(fbou{d})
    fbou{d} = reshape(fbou{d},ncu,[]);    
end
if ~isempty(uhat{d})
    uhat{d} = reshape(uhat{d},ncu,[]);
end
if ~isempty(fhat{d})
    fhat{d} = reshape(fhat{d},ncu,[]);    
end

gencppcode(app.appname,xdg,udg,odg,uhg,nlg,tau,uinf,param,time,...
    flux,source,tdfunc,ubou,fbou,uhat,fhat,avfd);

gencppcode2(app.appname,xdg,uhg,udg1,udg2,odg1,odg2,nlg,param,time,stab);

% gencppcode3(app.appname,xdg,odg,udg,udg0,udg1,udg2,udg3,uhg,uhg0,uhg1,uhg2,uhg3,nlg,...
%             tau,uinf,param,time,dt,app.tdepbc,app.nstage,uboutdep);

% function gencppcode3(appname,xdg,odg,udg,udg0,udg1,udg2,udg3,uhg,uhg0,uhg1,uhg2,uhg3,nlg,tau,...
%     uinf,param,time,dt,tdepbc,nstage,uboutdep)        
%         function gencppcode3(appname,xdg,odg,udg,udg0,udg1,udg2,udg3,uhg,uhg0,uhg1,uhg2,uhg3,nlg,tau,...
%     uinf,param,time,dt,tdepbc,tstage,uboutdep)

if ismac
% cpucompiler = '/opt/local/bin/g++-mp-7';
% mpicompiler = '/opt/local/bin/mpicxx-openmpi-mp';
cpucompiler = 'g++';
mpicompiler = 'mpicxx';
gpucompiler = [];
gpuflags    = [];
% blaslapack  = 'mkl';
ompblaslapack  = [];
blaslapack  = '-lblas -llapack';
elseif isunix
cpucompiler = 'g++';
mpicompiler = 'mpicxx';
gpucompiler = 'nvcc';    
gpuflags    = '-lcudart -lcublas';
blaslapack  = 'mkl';
ompblaslapack  = 'mkl';  
end
%cpuflags = '-pg -g';
cpuflags = ' ';
if isfield(app,'compilerflags') == 1
    cpuflags = [cpuflags ' ' app.compilerflags];
end

deletesourcefiles = 0;
compilerstr = compilecppcode(app.version,app.appname,cpucompiler,mpicompiler,gpucompiler,cpuflags,gpuflags,blaslapack,ompblaslapack,deletesourcefiles);

cd('..');
