function mkshape(porder,plocal,pts,elemtype)
#MKSHAPE Calculates the nodal shape functions and its derivatives for
#        various master elements
#
#   NFS=MKSHAPE(PORDER,PLOCAL,PTS,ELEMTYPE)
#
#      PORDER:    Polynomial order
#      PLOCAL:    Node positions (np,dim)
#      PTS:       Coordinates of the points where the shape fucntions
#                 and derivatives are to be evaluated (npoints,dim)
#      ELEMTYPE:  Flag determining element type
#                 Flag = 0 simplex elements (default)
#                 Flag = 1 tensor elements
#
#      NFS:       shape function and derivatives (np,npoints,dim+1)
#                 nsf(:,:,1) shape functions
#                 nsf(:,:,2) shape fucntions derivatives w.r.t. x
#                 nsf(:,:,3) shape fucntions derivatives w.r.t. y
#                 nsf(:,:,4) shape fucntions derivatives w.r.t. z
#


if elemtype==0 # simplex elements
    nf,nfx,nfy,nfz=simplexmonomial(pts,porder);   # orthogonal shape functions
    A,~,~,~ = simplexmonomial(plocal,porder);               # Vandermonde matrix
elseif elemtype==1           # tensor product elements
    nf,nfx,nfy,nfz=tensorproduct(pts,porder); # orthogonal shape functions
    A,~,~,~ = tensorproduct(plocal,porder);             # Vandermonde matrix
else
    error("This element type is not implemented");
end

# Divide orthogonal shape functions by the Vandemonde matrix to
# obtain nodal shape functions
dim = size(plocal,2);
if dim== 1
        nfs=[nf/A nfx/A];
        nfs=reshape(nfs,size(nf,1),size(nf,2),2);
        nfs=permutedims(nfs,[2,1,3]);
elseif dim==2
        nfs=[nf/A nfx/A nfy/A];
        nfs=reshape(nfs,size(nf,1),size(nf,2),3);
        nfs=permutedims(nfs,[2,1,3]);
elseif dim==3
        nfs=[nf/A nfx/A nfy/A nfz/A];
        nfs=reshape(nfs,size(nf,1),size(nf,2),4);
        nfs=permutedims(nfs,[2,1,3]);
else
        error("Only can handle dim=1, dim=2 or dim=3");
end

return nfs;

end
