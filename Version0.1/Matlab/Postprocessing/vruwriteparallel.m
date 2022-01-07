function vruwriteparallel(app, mesh, dmd)

[~,~,endian] = computer;
if strcmp(endian,'L')
    byte_order = 'LittleEndian';
else
    byte_order = 'BigEndian';
end

if isempty(app.viselem)
    ne = size(mesh.t,2);
    app.viselem = 1:ne;
end

% if app.porder>1
%     visorder = min(2*app.porder,8);
% else
%     visorder = app.porder;
% end
% if isfield(mesh, 'dgnodes')
%     visorder = app.porder;
% end
visorder = app.porder;
[xpe,telem] = masternodes(visorder,app.nd,app.elemtype);
% shape = mkshape(app.porder,mesh.xpe,xpe,app.elemtype);
% shape = shape(:,:,1)';

if isfield(mesh, 'dgnodes')
    dgnodes = mesh.dgnodes;
else
    dgnodes = createdgnodes(mesh.p,mesh.t(:,app.viselem),mesh.f(:,app.viselem),mesh.curvedboundary,mesh.curvedboundaryexpr,visorder);    
end
npe = size(dgnodes,1);
nc = app.nc;
size(dgnodes)

nproc = length(dmd);
nei = zeros(1,nproc);
for i = 1:nproc
    nei(i) = sum(dmd{i}.elempartpts(1:2));
end

% Write Header for Paraview PVD File
filename = app.visfilename;
pvdfile = filename + ".pvd";
fid = fopen(pvdfile, 'w');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<VTKFile type="Collection" version="0.1"\n');
fprintf(fid,'         byte_order="%s"\n',byte_order);
fprintf(fid,'         compressor="vtkZLibDataCompressor">\n');
fprintf(fid,'  <Collection>\n');

for i = 1:nproc
    fileID = fopen(['dataout/out' '_np' num2str(i-1) '.bin'],'r');       
    tm = fread(fileID,'double');
    udg = reshape(tm,[npe nc nei(i)]);
    
    elempart = dmd{i}.elempart(1:nei(i));   
    xdg = dgnodes(:,:,elempart);
    [cgnodes, cgelcon, cgcells, celltype] = createcggrid(xdg,telem);
       
    %tm = shape*reshape(udg(:,:,app.viselem),[size(shape,2) size(udg,2)*length(app.viselem)]);
    %tm = reshape(tm,[size(shape,1) size(udg,2) length(app.viselem)]);
    
    vtufile = filename + num2str(i);
    vtuwrite(vtufile, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, udg);       

    ind = strfind(vtufile, "/");
    ch = char(vtufile);
    outfile = string(ch((ind(end)+1):end));
    
    % Write pointer to .vtu file in .pvd file
    fprintf(fid,'    <DataSet subdomain="%e" group="" part="0"\n',i);
    fprintf(fid,'             file="%s"/>\n', outfile + ".vtu"); 
end

% Finalize and close .pvd file
fprintf(fid,'  </Collection>\n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);


