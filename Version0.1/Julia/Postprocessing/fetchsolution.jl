function fetchsolution(app,master,dmd)

nt = length(app.dt);
if nt==1
    Uout = getsolution("dataout/out",dmd,master.npe);
else
    if app.vistime != []
        app.vistime = 1:nt;
    end
    fn = "dataout/out_t" * string(app.vistime[1]);
    tmp = getsolution(fn,dmd,master.npe);
    Uout = zeros((size(tmp,1), size(tmp,2), size(tmp,3), length(app.vistime)));
    Uout[:,:,:,1] = tmp;
    tmp = [];
    for i = 2:length(app.vistime)
        fn = "dataout/out_t" * string(app.vistime[i]);
        Uout[:,:,:,i] = getsolution(fn,dmd,master.npe);
    end
end

return Uout

end
