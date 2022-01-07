function Uout = fetchsolution(app,master,dmd)

nt = length(app.dt);
if nt==1
    Uout = getsolution('dataout/out',dmd,master.npe);
else    
    if isempty(app.soltime)    
        app.soltime = 1:nt;
    end        
    tmp = getsolution(['dataout/out_t' num2str(app.soltime(1))],dmd,master.npe);
    if app.wave==1
        w = getsolution(['dataout/out_wdg_t' num2str(app.soltime(1))],dmd,master.npe);
        tmp = cat(2,tmp,w);
    end
    Uout = zeros([size(tmp) length(app.soltime)]);
    Uout(:,:,:,1) = tmp;
    clear tmp;
    for i = 2:length(app.soltime)
        tmp = getsolution(['dataout/out_t' num2str(app.soltime(i))],dmd,master.npe);
        if app.wave==1
            w = getsolution(['dataout/out_wdg_t' num2str(app.soltime(i))],dmd,master.npe);
            tmp = cat(2,tmp,w);
        end
        Uout(:,:,:,i) = tmp;
    end
end


