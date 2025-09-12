function read_rank(fname::AbstractString)
    # Open and ensure the stream is closed on error
    io = open(fname, "r")
    try
        # Header: three Float64 values [n1, n2, n3]
        hdr = Vector{Float64}(undef, 3)
        read!(io, hdr)
        
        if length(hdr) != 3
            error("Header read failed for $fname")
        end
        n1 = Int(hdr[1]); n2 = Int(hdr[2]); n3 = Int(hdr[3])
        N  = n1 * n2 * n3

        # Read the rest of the file and reinterpret bytes as Float64 (no copy)
        bytes   = read(io)                       # Vector{UInt8}
        payload = reinterpret(Float64, bytes)    # Vector{Float64} view

        # Sanity check: payload length must be a multiple of N
        if length(payload) % N != 0
            error("Payload size is not a multiple of n1*n2*n3 in $fname")
        end
        nsteps = length(payload) ÷ N

        # Reshape into [n1, n2, n3, nsteps]
        payload = reshape(payload, n1, n2, n3, nsteps)

        return n1, n2, n3, nsteps, payload
    finally
        close(io)
    end
end

# function read_rank(fname::AbstractString)
#     io = open(fname, "r")
#     try
#         hdr = Vector{Float64}(undef, 3)
#         read!(io, hdr)
#         n1 = Int(hdr[1]); n2 = Int(hdr[2]); n3 = Int(hdr[3])
#         N  = n1 * n2 * n3

#         data = Vector{Float64}(undef, N)
#         read!(io, data)
#         nd   = length(data)
#         if nd % N != 0
#             error("Payload size not a multiple of n1*n2*n3 in $fname (got $nd, block $N).")
#         end
#         nsteps = nd ÷ N
#         data4d = reshape(data, n1, n2, n3, nsteps)
#         return n1, n2, n3, nsteps, data4d
#     finally
#         close(io)
#     end
# end