%load dmd.mat;
%load solM07.mat;
%ind = 200:200:148600;

nn = 20000;
for i = 1:7
    inda{i} = ((i-1)*nn+200):200:i*nn;
end
inda{i} = 120200:200:148600;


mpiprocs = length(elempart);
npe = size(UDG,1);
ne = size(UDG,3);

for j=2:length(inda)
    ind = inda{j};    
    ns = length(ind);
    u2d = zeros(9,5,mesh2d.ne,ns);
    uav = u2d;
    i = 0;
    for itime=ind
        i = i + 1;  
        fn = ['/data/data1/oat15a3d/data1/nsout_t' num2str(itime)]
        tmp = getsolution(fn, mpiprocs, npe, 5, elempart);    
        u = reshape(tmp,[9 3 5 ne/32 32]);
        u2d(:,:,:,i) = squeeze(u(:,2,:,:,16)); 
        uav(:,:,:,i) = squeeze(sum(sum(u,5),2))/96;
        [i itime]
    end
    filen = ['u2d' num2str(j) '.mat'];
    save(filen, 'u2d', 'uav');
end





