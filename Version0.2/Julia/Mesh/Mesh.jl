__precompile__()

module Mesh

export squaremesh, cubemesh, gmshcall
export readmesh, writemesh

include("squaremesh.jl");
include("cubemesh.jl");
include("gmshcall.jl");
include("readmesh.jl");
include("writemesh.jl");

end
