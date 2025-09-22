mutable struct MESHStruct
    p::Array{FloatP,2};       # points of a linear mesh
    t::Array{IntP,2};         # elements of a linear mesh
    f::Array{IntP,2};         # faces of a linear mesh
    dgnodes::Array{FloatP,3}; # spatial nodes of a high-order mesh
    udg::Array{FloatP,3}; # spatial nodes of (u,q)
    odg::Array{FloatP,3}; # spatial nodes of o
    wdg::Array{FloatP,3}; # spatial nodes of w
    tprd::Array{IntP,2};      # elements for periodic conditions
    boundaryexpr;             # expressions to determine boundaries
    boundarycondition::Array{IntP,2}; # a list of boundary conditions
    curvedboundary::Array{IntP,2}; # flags for curved boundaries
    curvedboundaryexpr;       # expressions to determine curved boundaries
    periodicboundary::Array{IntP,2}; # flags for curved boundaries
    periodicexpr;             # expressions to determine periodic boundaries
    xpe::Array{FloatP,2};     # nodes on the master element
    xpf::Array{FloatP,2};     # nodes on the master face
    telem::Array{IntP,2};     # master element connectivity
    tface::Array{IntP,2};     # master face connectivity
    MESHStruct() = new();
end

function initializemesh(version)
    mesh = MESHStruct();
    mesh.boundaryexpr = [];
    mesh.periodicexpr = [];
    mesh.curvedboundaryexpr = [];
    mesh.curvedboundary = [0 0];
    mesh.periodicboundary = [0 0];
    mesh.boundarycondition = [0 0];
    mesh.dgnodes = zeros(0,0,0);
    mesh.udg = zeros(0,0,0);
    mesh.odg = zeros(0,0,0);
    mesh.wdg = zeros(0,0,0);
    return mesh;
end
