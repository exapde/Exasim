# Read only the first 10 Float64 values (header).
function readheader10(path::AbstractString)::NTuple{10,Float64}
    open(path, "r") do io
        hdr = Vector{Float64}(undef, 10)
        read!(io, hdr)
        return Tuple(hdr)  # small, stack-friendly
    end
end
