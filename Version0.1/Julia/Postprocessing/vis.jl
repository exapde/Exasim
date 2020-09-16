using Preprocessing, Gencode

function vis(visfields,app,mesh)

nt = size(visfields,4);

if app.viselem == []
    ne = size(mesh.t,2);
    app.viselem = 1:ne;
end

dgnodes = Preprocessing.createdgnodes(mesh.p,mesh.t[:,app.viselem],mesh.f[:,app.viselem],mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);
cgnodes, cgelcon, cgcells, celltype = createcggrid(dgnodes,mesh.telem);

# find paraview executable
app.paraview = Preprocessing.findexec(app.paraview, app.version);

# paraview = app.paraview;
# paraviewstatus0 = Sys.which(paraview);
# paraviewstatus1 = Sys.which("paraview");
# paraviewstatus2 = Sys.which("usr/bin/paraview");
# paraviewstatus3 = Sys.which("/usr/local/bin/paraview");
# paraviewstatus4 = Sys.which("/opt/local/bin/paraview");
#
# if paraviewstatus0 != nothing
# elseif paraviewstatus1 != nothing
#     paraview = "paraview"
# elseif paraviewstatus2 != nothing
#     paraview = "/usr/bin/paraview";
# elseif paraviewstatus3 != nothing
#     paraview = "/usr/local/bin/paraview";
# elseif paraviewstatus4 != nothing
#     paraview = "/opt/local/bin/paraview";
# else
#     error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Paraview. Please see the documentation to install it. After installation, please set its path to app.paraview");
# end
# app.paraview = paraview;

if nt==1
    vtuwrite(app.visfilename, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, visfields[:,:,app.viselem]);
    if length(app.paraview)>0
        str = app.paraview * " --data="* app.visfilename * ".vtu";
        run(Gencode.string2cmd(str), wait=false);
    end
else
    if length(app.visdt)==0
        app.visdt = 1;
    end
    pvdwrite(app.visfilename, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, visfields[:,:,app.viselem,:],app.visdt);
    if length(app.paraview)>0
        str = app.paraview * " --data=" * app.visfilename * ".pvd";
        run(Gencode.string2cmd(str), wait=false);
    end
end

return dgnodes

end
