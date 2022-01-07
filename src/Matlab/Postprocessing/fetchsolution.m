function Uout = fetchsolution(app,master,dmd,dirname)

if nargin<4
    dirname = 'dataout';
end

nt = length(app.dt);
if nt==1
    Uout = getsolution([dirname '/out'],dmd,master.npe);
else    
    if isempty(app.soltime)    
        app.soltime = 1:nt;
    end        
    tmp = getsolution([dirname '/out_t' num2str(app.soltime(1))],dmd,master.npe);
    if app.wave==1
        w = getsolution([dirname '/out_wdg_t' num2str(app.soltime(1))],dmd,master.npe);
        tmp = cat(2,tmp,w);
    end
    Uout = zeros([size(tmp) length(app.soltime)]);
    Uout(:,:,:,1) = tmp;
    clear tmp;
    for i = 2:length(app.soltime)
        tmp = getsolution([dirname '/out_t' num2str(app.soltime(i))],dmd,master.npe);
        if app.wave==1
            w = getsolution([dirname '/out_wdg_t' num2str(app.soltime(i))],dmd,master.npe);
            tmp = cat(2,tmp,w);
        end
        Uout(:,:,:,i) = tmp;
    end
end


