__precompile__()

module Postprocessing

using Gencode, Preprocessing

export fetchsolution, vis, exasim, producecode

include("createcgcells.jl");
include("createcggrid.jl");
include("getcelltype.jl");
include("getsolution.jl");
include("getmeansolution.jl");
include("fetchsolution.jl");
include("vtuwrite.jl");
include("pvdwrite.jl");
include("vis.jl");
include("exasim.jl");
include("producecode.jl");
include("generatecode.jl");

end
