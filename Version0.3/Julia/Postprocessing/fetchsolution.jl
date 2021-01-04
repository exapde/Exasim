function fetchsolution(app,master,dmd,dirname)

nt = length(app.dt);
if nt==1
    Uout = getsolution(dirname * "/out",dmd,master.npe);
else
    if app.soltime != []
        app.soltime = 1:nt;
    end
    fn = dirname * "/out_t" * string(app.soltime[1]);
    tmp = getsolution(fn,dmd,master.npe);
    if app.wave==1
        w = getsolution(dirname * "/out_wdg_t" * string(app.soltime[1]),dmd,master.npe);
        tmp = hcat(tmp,w);
    end
    Uout = zeros((size(tmp,1), size(tmp,2), size(tmp,3), length(app.soltime)));
    Uout[:,:,:,1] = tmp;
    for i = 2:length(app.soltime)
        fn = dirname * "/out_t" * string(app.soltime[i]);
        tmp = getsolution(fn,dmd,master.npe);
        if app.wave==1
            w = getsolution(dirname * "/out_wdg_t" * string(app.soltime[i]),dmd,master.npe);
            tmp = hcat(tmp,w);
        end
        Uout[:,:,:,i] = tmp;
    end
end

return Uout

end
