function dgnodes = vis(visfields,app,mesh)

nt = size(visfields,4);

if isempty(app.viselem)
    ne = size(mesh.t,2);
    app.viselem = 1:ne;
end

dgnodes = createdgnodes(mesh.p,mesh.t(:,app.viselem),mesh.f(:,app.viselem),mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);    
[cgnodes, cgelcon, cgcells, celltype] = createcggrid(dgnodes,mesh.telem);

% find paraview executable
app.paraview = findexec(app.paraview, app.version);

% paraview = app.paraview;
% [paraviewstatus0,~] = system("which " + paraview);
% [paraviewstatus1,~] = system("which paraview");
% [paraviewstatus2,~] = system("which /usr/bin/paraview");
% [paraviewstatus3,~] = system("which /usr/local/bin/paraview");
% [paraviewstatus4,~] = system("which /opt/local/bin/paraview");        
% 
% if paraviewstatus0==0
% elseif paraviewstatus1==0        
%     paraview = "paraview";        
% elseif paraviewstatus2==0        
%     paraview = "/usr/bin/paraview";    
% elseif paraviewstatus3==0        
%     paraview = "/usr/local/bin/paraview";    
% elseif paraviewstatus4==0        
%     paraview = "/opt/local/bin/paraview";    
% else            
%     error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Paraview. Please see the documentation to install it. After installation, please set its path to app.paraview"); 
% end
% app.paraview = paraview;

if nt==1
    vtuwrite(app.visfilename, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, visfields(:,:,app.viselem));    
    if ~isempty(app.paraview)        
        str = app.paraview + " --data='" + app.visfilename + ".vtu'" + " &";           
        eval(char("!" + str));
    end    
else
    if isempty(app.visdt)    
        app.visdt = 1;
    end    
    pvdwrite(app.visfilename, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, visfields(:,:,app.viselem,:),app.visdt);
    if ~isempty(app.paraview) 
        str = app.paraview + " --data='" + app.visfilename + ".pvd'" + " &";
        eval(char("!" + str));
    end    
end
