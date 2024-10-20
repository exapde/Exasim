function [C] = mapContractK(A,B,ia,ka,la,kb,jb,lb)
% A is a multi-array and depends on multi-indices ia,ka and la.
% B is a multi-array and depends on multi-indices kb,jb and lb.
% C is the tensor contraction that depends only on ia,jb,la (ka,kb are
%   contracted).
% ia determine the position of indices that are in A and not in B.
% jb determine the position of indices that are in B and not in A.
% ka and kb are the position of indices shared by A and B. These indices
%   are contracted.
% la and lb are the position of the rest of indices shared by A and B.
%   These indices are not contracted.
%
% C(ia,jb,lb) = Sum_k A(ia,ka,la)B(kb,jb,lb)
%   for la1 = 1,...,ra1; ... ; lah = 1,...,rah, where la = (la1,...,lah).
%
% Where:
% multi-index ia = (ia1,...,iae), e is the number of components.
% multi-index jb = (jb1,...,jbf), f is the number of components.
% multi-index ka = (ka1,...,kag), g is the number of components.
% multi-index kb = (kb1,...,kbg), g is the number of components.
% multi-index la = (la1,...,lah), h is the number of components.
% multi-index lb = (lb1,...,lbh), h is the number of components.
%
% The multi-index ia covers ia1 = 1,...,ma1; ... ; iae = 1,...,mae.
% The multi-index jb covers jb1 = 1,...,nb1; ... ; ibf = 1,...,nbf.
% The multi-index ka covers ka1 = 1,...,qa1; ... ; kag = 1,...,qag.
% The multi-index kb covers kb1 = 1,...,qb1; ... ; kbg = 1,...,qbg.
% The multi-index la covers la1 = 1,...,ra1; ... ; lah = 1,...,rah.
% The multi-index lb covers lb1 = 1,...,rb1; ... ; lbh = 1,...,rbh.
%
% The order of (ia1,...,iae) determines the order of ia in C.
% The order of (jb1,...,jbf) determines the order of jb in C.
% The order of (ka1,...,kag) have to be the same order
%   of (kb1,...,kbg). Thus, qa1 = qb1, ... , qag = qbg.
% The order of (la1,...,lah) determines the order of la in C. Moreover,
%   it has to be the same order of (lb1,...,lbh). Thus,
%   ra1 = rb1, ... , rah = rbh.
% The size of A is ( ma1 X...X mae X qa1 X...X qag ) X ra1 X ... X rah.
% The size of B is ( qb1 X...X qbg X nb1 X...X nbf ) X rb1 X ... X rbh.
% The size of C is ( ma1 X...X mae X nb1 X...X nbf ) X ra1 X ... X rah.
% Note that qa1 X...X qag = qb1 X...X qbg and
%   la1 X ... X lah = lb1 X ... X lbh.
%
% ia, ka and la should determine the positions of all indices of A
% kb, jb and lb should determine the positions of all indices of B

sA = size(A);
sB = size(B);

Ap = my_permute(A,[ia ka la]);
Bp = my_permute(B,[kb jb lb]);

lenA = max([ia ka la]);
lenB = max([kb jb lb]);

if numel(sA) < lenA
    sA = [sA ones(1,lenA-numel(sA))];
end

if numel(sB) < lenB
    sB = [sB ones(1,lenB-numel(sB))];
end

ma = sA(ia);
qa = sA(ka);
ra = sA(la);

qb = sB(kb);
nb = sB(jb);
rb = sB(lb);

Ac = reshape(Ap,[prod(ma) prod(qa) prod(ra)]);
Bc = reshape(Bp,[prod(qb) prod(nb) prod(rb)]);

Cc = mapContract_ikl_kjl_ijl(Ac,Bc); 
%Cc = arraymul(Ac,Bc);

C = reshape(Cc,[ma nb ra 1 1]);


function Ap = my_permute(A,pA)
nA = numel(pA);
if min(1:nA == pA) == 0
    Ap = permute(A,pA);
else
    Ap = A;
end


function C = mapContract_ikl_kjl_ijl(A,B)
% tic;

r = size(A,3);

C = zeros(size(A,1),size(B,2),size(A,3));

for l = 1:r
    C(:,:,l) = A(:,:,l)*B(:,:,l);
end

% t = toc;
% 
%  m = size(A,1);
%  q = size(A,2);
%  n = size(B,2);
%  
%  s = [m q n r]
% 
% numOps = r*m*n*2*q;
% 
% GIGAFLOPS = numOps/(10^9*t)
