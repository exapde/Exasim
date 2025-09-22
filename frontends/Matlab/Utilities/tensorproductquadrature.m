
nd = 3;
porder = 3;
pgauss = 3*porder;
elemtype = 1;
nodetype = 1;

master=mkmasterelement(nd,porder,porder,pgauss,pgauss,elemtype,nodetype);
p1d = mkmasternodes(porder,1,elemtype,nodetype);
g1d = gaussquad(pgauss,1,elemtype);
shap1d = mkshape(porder,p1d,g1d,elemtype);
a=shap1d(:,:,1)'; % node-2-gauss
at = a';

if nd==2    
    b=master.shapvt(:,:,1); % node-2-gauss
    c=master.shapvl(:,:,1); % gauss-2-node
    e=b-kron(a,a);
    max(abs(e(:)))
    e=c-kron(at,at);
    max(abs(e(:)))
    
    l = size(a,1); % # gauss points in 1D
    k = size(a,2); % # shape functions in 1D
    n = k*k;       % # shape functions in 2D
    p = l*l;       % # gauss points in 2D
    
    % node-2-gauss
    m = 100;
    u = rand(n,m);
    v = b*u;
    
    w = a*reshape(u,[k k*m]); 
    w = reshape(w, [l k m]);
    w = permute(w, [2 1 3]);
    w = a*reshape(w, [k l*m]);    
    w = reshape(w, [l l m]);
    w = permute(w, [2 1 3]);
    e = v(:)-w(:);
    max(abs(e(:)))    
    
    % gauss-2-node
    u = rand(p,m);
    v = c*u;
    
    w = at*reshape(u,[l l*m]); 
    w = reshape(w, [k l m]);
    w = permute(w, [2 1 3]);
    w = at*reshape(w, [l k*m]);    
    w = reshape(w, [k k m]);
    w = permute(w, [2 1 3]);
    e = v(:)-w(:);
    max(abs(e(:)))    
else
    b=master.shapvt(:,:,1); % node-2-gauss
    c=master.shapvl(:,:,1); % gauss-2-node
    e=b-kron(kron(a,a),a);
    max(abs(e(:)))
    e=c-kron(kron(at,at),at);
    max(abs(e(:)))
    
    l = size(a,1); % # gauss points in 1D
    k = size(a,2); % # shape functions in 1D
    n = k*k*k;     % # shape functions in 3D
    p = l*l*l;     % # gauss points in 3D
    
    % node-2-gauss
    m = 100;
    u = rand(n,m);
    v = b*u;
    
    w = a*reshape(u,[k k*k*m]); 
    w = reshape(w, [l k k*m]);
    w = permute(w, [2 1 3]);
    w = a*reshape(w, [k l*k*m]);    
    w = reshape(w, [l*l k m]);
    w = permute(w, [2 1 3]);
    w = a*reshape(w, [k l*l*m]);    
    w = reshape(w, [l l l m]);
    w = permute(w, [3 2 1 4]);
    e = v(:)-w(:);
    max(abs(e(:)))    
    
    % gauss-2-node
    u = rand(p,m);
    v = c*u;
    
    w = at*reshape(u,[l l*l*m]); 
    w = reshape(w, [k l l*m]);
    w = permute(w, [2 1 3]);
    w = at*reshape(w, [l k*l*m]);    
    w = reshape(w, [k*k l m]);
    w = permute(w, [2 1 3]);
    w = at*reshape(w, [l k*k*m]);    
    w = reshape(w, [k k k m]);
    w = permute(w, [3 2 1 4]);
    e = v(:)-w(:);
    max(abs(e(:)))        
end



