function getsolutions(basename::AbstractString, dmd::AbstractVector)
    nproc = length(dmd)

    if nproc == 1
        n1, n2, n3, nsteps, block =
            read_rank(string(basename, "_np", 0, ".bin"))
        display([n1, n2, n3, nsteps])    
        return block  # size: (n1, n2, n3, nsteps)
    else
        # how many elements each rank contributes
        nei = zeros(Int, nproc)
        for r in 1:nproc
            # assume dmd[r].elempartpts is an array where sum of first two is local ne
            nei[r] = sum(dmd[r].elempartpts[1:2])
        end
        net = sum(nei)

        # read rank 0 for sizing and initial placement
        n1, n2, _, nsteps, block0 =
            read_rank(string(basename, "_np", 0, ".bin"))
        sol = zeros(Float64, n1, n2, net, nsteps)

        # place rank 0
        elempart = dmd[1].elempart[1:nei[1]]
        @views sol[:,:,elempart,:] .= block0

        # place ranks 1..nproc-1
        for r in 2:nproc
            _, _, _, _, blockr =
                read_rank(string(basename, "_np", r-1, ".bin"))
            elempart = dmd[r].elempart[1:nei[r]]
            @views sol[:,:,elempart,:] .= blockr
        end

        return sol
    end
end

# function getsolutions(pde, dmd)
#     nproc = length(dmd)    
#     if nproc === 1
#         # ---------- Single-process ----------
#         path0 = joinpath(pde.buildpath, "dataout", "outsol_np0.bin")
#         tmp   = readdoubles(path0)

#         npe  = Int(tmp[1])
#         nc   = Int(tmp[2])
#         ne   = Int(tmp[3])
#         ncw  = Int(tmp[5])
#         ncu  = Int(tmp[7])
#         npf  = Int(tmp[8])
#         nf   = Int(tmp[9])
#         N    = Int(tmp[10])

#         payload = length(tmp) - 10
#         @assert payload % N == 0 "File payload not a multiple of N."
#         nsteps = payload ÷ N

#         UDG  = zeros(Float64, npe, nc,  ne,  nsteps)
#         WDG  = zeros(Float64, npe, ncw, ne,  nsteps)
#         UHAT = zeros(Float64, ncu, npf, nf,  nsteps)

#         n1 = npe*nc*ne
#         n2 = npe*ncw*ne
#         n3 = ncu*npf*nf

#         m = 11
#         for i in 1:nsteps
#             UDG[:,:,:,i]  = reshape(view(tmp, m:(m+n1-1)),  npe,  nc,  ne);  m += n1
#             WDG[:,:,:,i]  = reshape(view(tmp, m:(m+n2-1)),  npe,  ncw, ne);  m += n2
#             UHAT[:,:,:,i] = reshape(view(tmp, m:(m+n3-1)),  ncu,  npf, nf);  m += n3
#         end

#         return UDG, WDG, UHAT
#     else
#         # ---------- Multi-process ----------
        
#         # (1) Rank 0: read full file once to discover sizes, and nf0 for faces
#         path0  = joinpath(pde.buildpath, "dataout", "outsol_np0.bin")
#         tmp0   = readdoubles(path0)
#         npe    = Int(tmp0[1])
#         nc     = Int(tmp0[2])
#         ne0    = Int(tmp0[3])
#         ncw    = Int(tmp0[5])
#         ncu    = Int(tmp0[7])
#         npf    = Int(tmp0[8])
#         nf0    = Int(tmp0[9])
#         N0     = Int(tmp0[10])

#         payload0 = length(tmp0) - 10
#         @assert payload0 % N0 == 0 "Rank 0: payload not a multiple of N0."
#         nsteps   = payload0 ÷ N0

#         # (2) Elements per rank from dmd
#         nei = [sum(dmd[r].elempartpts[1:2]) for r in 1:nproc]  # local ne (rank view)
#         net = sum(nei)                                         # total elements

#         # (3) Faces per rank: read just headers for ranks 2..nproc to sum nf
#         #     (we already have nf0 from rank 0)
#         nf_others = 0
#         for r in 2:nproc
#             path = joinpath(pde.buildpath, "dataout", "outsol_np$(r-1).bin")
#             hdr  = readheader10(path)
#             # sanity checks for shared header fields (optional but recommended)
#             @assert Int(hdr[1]) == npe  "Rank $(r-1): npe mismatch."
#             @assert Int(hdr[2]) == nc   "Rank $(r-1): nc mismatch."
#             @assert Int(hdr[5]) == ncw  "Rank $(r-1): ncw mismatch."
#             @assert Int(hdr[7]) == ncu  "Rank $(r-1): ncu mismatch."
#             @assert Int(hdr[8]) == npf  "Rank $(r-1): npf mismatch."
#             # accumulate faces
#             nf_others += Int(hdr[9])
#             # also ensure each rank has same number of steps:
#             Nr = Int(hdr[10])
#             bytes = filesize(path)
#             count = Int(bytes ÷ 8)         # doubles
#             payload = count - 10
#             @assert payload % Nr == 0      "Rank $(r-1): payload not multiple of N."
#             @assert (payload ÷ Nr) == nsteps "Rank $(r-1): step count differs from rank 0."
#         end
#         nft = nf0 + nf_others  # total faces across all ranks

#         # (4) Allocate global arrays
#         UDG  = zeros(Float64, npe, nc,  net, nsteps)
#         WDG  = zeros(Float64, npe, ncw, net, nsteps)
#         UHAT = zeros(Float64, ncu, npf, nft, nsteps)

#         # (5) Place rank 0 data
#         elempart0 = dmd[1].elempart[1:nei[1]]
#         n1_0 = npe*nc*ne0
#         n2_0 = npe*ncw*ne0
#         n3_0 = ncu*npf*nf0
#         m = 11
#         for i in 1:nsteps
#             UDG[:,:,elempart0,i]  = reshape(view(tmp0, m:(m+n1_0-1)), npe,  nc,  ne0);   m += n1_0
#             WDG[:,:,elempart0,i]  = reshape(view(tmp0, m:(m+n2_0-1)), npe,  ncw, ne0);   m += n2_0
#             UHAT[:,:,(1:nf0),i]   = reshape(view(tmp0, m:(m+n3_0-1)), ncu,  npf, nf0);   m += n3_0
#         end

#         # (6) Subsequent ranks: append face blocks and place element blocks by mapping
#         k = nf0
#         for r in 2:nproc
#             path = joinpath(pde.buildpath, "dataout", "outsol_np$(r-1).bin")
#             tmp  = readdoubles(path)

#             ne_r = Int(tmp[3])
#             nf_r = Int(tmp[9])
#             N_r  = Int(tmp[10])

#             # rank-local block sizes
#             n1 = npe*nc*ne_r
#             n2 = npe*ncw*ne_r
#             n3 = ncu*npf*nf_r

#             elems = dmd[r].elempart[1:nei[r]]

#             m = 11
#             for i in 1:nsteps
#                 UDG[:,:,elems,i]        = reshape(view(tmp, m:(m+n1-1)), npe,  nc,  ne_r);  m += n1
#                 WDG[:,:,elems,i]        = reshape(view(tmp, m:(m+n2-1)), npe,  ncw, ne_r);  m += n2
#                 UHAT[:,:,(k+1):(k+nf_r),i] = reshape(view(tmp, m:(m+n3-1)), ncu,  npf, nf_r); m += n3
#             end
#             k += nf_r
#         end

#         return UDG, WDG, UHAT
#     end
# end


