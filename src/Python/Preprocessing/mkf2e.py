from numpy import *
from getelemface import *
from sortrows import *

def mkf2e(t,elemtype,nd):

    ne = t.shape[1]
    face = getelemface(nd,elemtype)
    nvf = face.shape[0]
    nfe = face.shape[1]

    # sort faces for all elements to make two identical faces have the same numbering order
    N = nfe*ne;
    tf = reshape(t[face.flatten(order='F'),:], (nvf, N), order='F');
    tf = sort(tf, axis=0);

    tf, jx = sortrows(tf.T); # do this so that two identical faces are next each other
    tf = tf.T;
    dx = tf[:,1:N]-tf[:,0:(N-1)]; # do this to find two identical faces
    dx = sum(dx*dx,axis=0);           # identical faces corresponds to zero entries
    in1 = nonzero(dx.flatten('F')==0)[0];           # interior global faces connected to 1st elements (1st identical faces)
    in1 = in1.flatten('F');
    in2 = in1 + 1;                 # interior global faces connected to 2nd elements (2nd identical faces)

    for i in range(0,len(in1)):
        if jx[in1[i]] > jx[in2[i]]:
            tm = jx[in1[i]];
            jx[in1[i]] = jx[in2[i]];
            jx[in2[i]] = tm;

    in3 = concatenate([in1, in2]);
    in0 = setdiff1d(arange(0,N), unique(in3)); # boundary global faces connected to 1st elements

    nf = len(in0)+len(in1);
    f2e = zeros((4,nf)).astype(int);   # allocate memory
    e2e = -ones((nfe,ne)).astype(int); # allocate memory

    # interior faces
    e1 = int64(ceil(float64(jx[in1]+1)/float(nfe))) - 1;      # 1st elements
    l1 = jx[in1] - e1*nfe;   # 1st local faces
    e2 = int64(ceil(float64(jx[in2]+1)/float(nfe))) - 1;      # 2nd elements
    l2 = jx[in2] - e2*nfe;   # 2nd local faces
    g = arange(0,len(in1));           # indices for interior faces

    f2e[0,g] = e1+1;
    f2e[1,g] = l1+1;
    f2e[2,g] = e2+1;
    f2e[3,g] = l2+1;
    for i in range(0,len(e1)):
        e2e[l1[i],e1[i]] = e2[i];
        e2e[l2[i],e2[i]] = e1[i];

    # boundary faces
    e1 = int64(ceil(float64(jx[in0]+1)/float(nfe))) - 1;      # 1st elements
    l1 = jx[in0] - e1*nfe;   # 1st local faces
    g = arange(len(in1),nf);      # indices for boundary faces
    f2e[0,g] = e1+1;
    f2e[1,g] = l1+1;

    # print(f2e.T)
    # print(e2e.T)
    # error("here")

    return f2e, e2e
