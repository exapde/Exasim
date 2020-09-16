function Uout = fetchsolution(app,master,dmd)

nt = length(app.dt);
if nt==1
    Uout = getsolution('dataout/out',dmd,master.npe);
else    
    if isempty(app.soltime)    
        app.soltime = 1:nt;
    end    
    tmp = getsolution(['dataout/out_t' num2str(app.soltime(1))],dmd,master.npe);
    Uout = zeros([size(tmp) length(app.soltime)]);
    Uout(:,:,:,1) = tmp;
    clear tmp;
    for i = 2:length(app.soltime)
        Uout(:,:,:,i) = getsolution(['dataout/out_t' num2str(app.soltime(i))],dmd,npe);
    end
end

