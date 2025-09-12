__precompile__()

module Postprocessing

using Gencode, Preprocessing

export fetchsolution, getsolutions, vis, exasim, producecode

include("createcgcells.jl");
include("createcggrid.jl");
include("getcelltype.jl");
#include("readheader10.jl");
#include("readdoubles.jl");
#include("getsolution.jl");
include("read_rank.jl");
include("getsolutions.jl");
include("getmeansolution.jl");
include("fetchsolution.jl");
include("fetchresidual.jl");
include("vtuwrite.jl");
include("pvdwrite.jl");
include("vis.jl");
include("exasim.jl");
include("producecode.jl");
include("generatecode.jl");

end
