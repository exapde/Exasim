__precompile__()

module Types

export IntP, FloatP
export DMDStruct

const IntP = Int64
const FloatP = Float64

mutable struct DMDStruct
    facecon::Array{IntP,3}; # faces-to-DG nodes connectivities (used to obtain DG nodes on faces)
    eblks::Array{IntP,2};  # blocks of elements for parallel computation
    fblks::Array{IntP,2};  # blocks of faces for parallel computation
    nbsd::Array{IntP,2};   # neighboring subdomains (neighbors)
    elemsend::Array{IntP,2};  # elements to be sent to neighbors
    elemrecv::Array{IntP,2};  # elements to be received from neighbors
    elemsendpts::Array{IntP,2}; # markers for elements to be sent to neighbors
    elemrecvpts::Array{IntP,2}; # markers for elements to be received from neighbors
    elempart::Array{IntP,2};    # element partitions
    elempartpts::Array{IntP,2}; # markers for element partitions
    facepart::Array{IntP,2};    # element partitions
    facepartpts::Array{IntP,2}; # markers for element partitions
    facepartbnd::Array{IntP,2};
    t2f::Array{IntP,2};
    f::Array{IntP,2};
    bf::Array{IntP,2};
    f2t::Array{IntP,2};
    elemcon::Array{IntP,3};
    #elem2cpu::Array{IntP,2};    # element partitions
    DMDStruct() = new();
end

end
