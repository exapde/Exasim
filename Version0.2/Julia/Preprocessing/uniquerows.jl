function uniquerows(a)

ind = findall(abs.(a[:]).<1e-10);
a[ind] .= 0;

n = size(a,1);
ia = unique(i -> a[i,:], 1:n);
c = a[ia,:];
ic = [findfirst(i->all(j->abs(a[k,j] - c[i,j])<1e-8,1:size(c,2)),1:size(c,1)) for k = 1:size(a,1)];

# d,tm = sortrows(c);
# e = d[2:end,:]-d[1:end-1,:];
# e = sum(e.*e,dims=2);
# i = argmin(e[:]);
# display(d[i-1:i+1,:])
# display(minimum(e[:]))

# display(size(c))
# display(size(a))
# display(size(ic))
#
# ent = unique(ic[:]);
# ndof = length(ent);
# display([minimum(ic) maximum(ic) ndof])

# ii = setdiff(1:n,ia);
# nc = length(ia);
# ic = zeros(Int,n);
# ic[ia] = 1:nc;
#
# #ic[ii] = xiny(a[ii,:],c,1);
#
# # nb = length(ii);
# # for k = 1:nb
# #     for q = 1:nc
# #         if maximum(abs.(a[ii[k],:]-c[q,:]))<1e-7
# #             ic[ii[k]] = q;
# #             break;
# #         end
# #     end
# # end
#
# while (true)
#     tm = a[ii,:];
#     display(size(tm))
#     if length(ii)==1
#         ij = 1;
#         b = reshape(tm[1,:],1,length(tm));
#     else
#         ij = unique(ix -> tm[ix,:], 1:size(tm,1));
#         b = tm[ij,:];
#     end
#     # dsiplay(ii)
#     # display(size(a))
#     # display(size(tm))
#
#     nb = size(b,1);
#     ind = zeros(Int,nb);
#     for i = 1:nc
#         if maximum(abs.(b[1,:]-c[i,:]))<1e-7
#             ind[1] = i;
#             break;
#         end
#     end
#     m = ind[1];
#     for k = 2:nb
#         found = 0;
#         for q = m:nc
#             if maximum(abs.(b[k,:]-c[q,:]))<1e-7
#                 ind[k] = q;
#                 m = q;
#                 found = 1;
#                 break;
#             end
#         end
#         if found==0
#             for q = 1:m-1
#                 if maximum(abs.(b[k,:]-c[q,:]))<1e-7
#                     ind[k] = q;
#                     m = q;
#                     break;
#                 end
#             end
#         end
#     end
#     if length(ii)==1
#         ic[ii] = ind;
#     else
#         ic[ii[ij]] = ind;
#     end
#
#     i1 = setdiff(1:length(ii),ij);
#     if length(i1)==0
#         break;
#     else
#         ii = ii[i1];
#     end
# end

# ii = setdiff(1:n,ia);
# nc = length(ia);
# ic = zeros(Int,n);
# ic[ia] = 1:nc;
#
# b,ij = sortrows(a[ii,:]);
# m = 0;
# for i = 1:nc
#     if maximum(abs.(b[1,:]-c[i,:]))<1e-7
#         m = i;
#         break;
#     end
# end
# ic[ii[ij[1]]] = m;
#
# for k = 2:length(ii)
#     if maximum(abs.(b[k,:]-b[k-1,:]))<1e-7
#         ic[ii[ij[k]]] = ic[ii[ij[k-1]]];
#     else
#         for q = 1:nc
#             if maximum(abs.(b[k,:]-c[q,:]))<1e-7
#                 ic[ii[ij[k]]] = q;
#                 m = q;
#                 break;
#             end
#         end
#     end
# end

return c,ia,ic

end
