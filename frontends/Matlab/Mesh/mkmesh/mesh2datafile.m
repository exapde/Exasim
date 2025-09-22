function mesh2datafile(mesh,master,app,UDG,UH,filename)

% open binary file to write
fid = fopen(filename,'wb');

% mesh structure
fwritearray(fid, int64(mesh.nd));
fwritearray(fid, int64(mesh.porder));
fwritearray(fid, int64(mesh.elemtype));
fwritearray(fid, int64(mesh.nodetype));
fwritearray(fid, mesh.dgnodes);
fwritearray(fid, int64(mesh.bf));
fwritearray(fid, int64(mesh.elcon));
fwritearray(fid, int64(mesh.f));
fwritearray(fid, int64(mesh.t2f));

% master structure
fwritearray(fid, master.shapmv);
fwritearray(fid, master.shapvt);
fwritearray(fid, master.shapvg);
fwritearray(fid, master.shapvgdotshapvl);
fwritearray(fid, master.shapmf);
fwritearray(fid, master.shapft);
fwritearray(fid, master.shapfg);
fwritearray(fid, master.shapfgdotshapfc);
fwritearray(fid, int64(master.perm));
fwritearray(fid, int64(master.permgeom));

% solution structure
fwritearray(fid, UDG);
fwritearray(fid, UH);

% app structure
fwritearray(fid, cell2mat(app.arg));
fwritearray(fid, int64(app.bcm));
fwritearray(fid, app.bcs);
fwritearray(fid, app.time);
fwritearray(fid, app.fc_q);
fwritearray(fid, app.fc_u);
fwritearray(fid, app.fc_p);
fwritearray(fid, int64(app.localsolve));
fwritearray(fid, int64(app.tdep));
fwritearray(fid, int64(app.wave));
fwritearray(fid, int64(app.flg_p));
fwritearray(fid, int64(app.adjoint));
fwritearray(fid, int64(app.bcd));
fwritearray(fid, app.bcv);

time = app.time;
tdep = app.tdep;
fc_q = app.fc_q;
fc_u = app.fc_u;
flg_p = app.flg_p;
localsolve = app.localsolve;
adjoint = app.adjoint;
source = str2func(app.source);
arg  = app.arg;
time = app.time;
bcm  = app.bcm;        
bcs  = app.bcs;
bcd  = app.bcd;        
bcv  = app.bcv;
localsolve = app.localsolve;
adjoint = app.adjoint;
fbou   = str2func(app.fbou);


% close the file
fclose(fid);

% check
fid = fopen(filename,'r');
[dimp,szp,clnp,p] = freadarray(fid);
[dimt,szt,clnt,t] = freadarray(fid);
fclose(fid);

max(abs(p(:)-mesh.p(:)))
max(abs(t(:)-mesh.t(:)))
[dimp,szp,clnp]
[dimt,szt,clnt]

function fwritearray(fid, a)

sz  = uint64(size(a));
dim = uint64(length(sz));

if isa(a,'double')
    cln = 0;
elseif isa(a,'int64')
    cln = 1;
elseif isa(a,'uint64')
    cln = 1;
elseif isa(a,'logical')
    cln = 2;    
else
    str = ['does not support this datatype' cl];
    error(str);
end
cln = uint64(cln);

fwrite(fid,dim,'uint64');
fwrite(fid,sz,'uint64');
fwrite(fid,cln,'uint64');
fwrite(fid,a,class(a));

function [dim,sz,cln,a] = freadarray(fid)

dim = fread(fid,1,'uint64');
sz = fread(fid,dim,'uint64');
cln = fread(fid,1,'uint64');
sz = reshape(sz,1,dim);

if cln == 0 % double
    cl = class(1.0);
elseif cln == 1 % 64-bit integer
    cl = class(int64(1));
elseif cln == 2 % boolean
    cl = class(false);
else
    str = 'does not support this datatype';
    error(str);
end

a = fread(fid,sz,cl);
