using Preprocessing, Gencode

function vis(visfields,app,mesh)

nt = size(visfields,4);

if app.viselem == []
    ne = size(mesh.t,2);
    app.viselem = 1:ne;
end

if app.porder>1
    visorder = min(2*app.porder,8);
else
    visorder = app.porder;
end

xpe,telem,~,~,~ = Preprocessing.masternodes(visorder,app.nd,app.elemtype);
shape = Preprocessing.mkshape(app.porder,mesh.xpe,xpe,app.elemtype);
shape = shape[:,:,1]';

dgnodes = Preprocessing.createdgnodes(mesh.p,mesh.t[:,app.viselem],mesh.f[:,app.viselem],mesh.curvedboundary,mesh.curvedboundaryexpr,visorder);
cgnodes, cgelcon, cgcells, celltype = createcggrid(dgnodes,telem);
dgnodes = Preprocessing.createdgnodes(mesh.p,mesh.t[:,app.viselem],mesh.f[:,app.viselem],mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);

# find paraview executable
app.paraview = Preprocessing.findexec(app.paraview, app.version);

if nt==1
    tm = shape*reshape(visfields[:,:,app.viselem],(size(shape,2), size(visfields,2)*length(app.viselem)));
    tm = reshape(tm,(size(shape,1), size(visfields,2), length(app.viselem)));
    vtuwrite(app.visfilename, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, tm);
    if length(app.paraview)>0
        str = app.paraview * " --data="* app.visfilename * ".vtu";
        run(Gencode.string2cmd(str), wait=false);
    end
else
    if length(app.visdt)==0
        app.visdt = 1;
    end
    pvdwrite(app.visfilename, cgnodes, cgelcon, cgcells, celltype, app.visscalars, app.visvectors, visfields[:,:,app.viselem,:],app.visdt,shape);
    if length(app.paraview)>0
        str = app.paraview * " --data=" * app.visfilename * ".pvd";
        run(Gencode.string2cmd(str), wait=false);
    end
end

return dgnodes

end
