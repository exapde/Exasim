# Read entire file as Float64 vector.
function readdoubles(path::AbstractString)::Vector{Float64}
    nbytes = filesize(path)
    @assert nbytes % 8 == 0 "File $path size is not a multiple of 8."
    n = Int(nbytes ÷ 8)
    open(path, "r") do io
        read!(io, Vector{Float64}(undef, n))
    end
end
