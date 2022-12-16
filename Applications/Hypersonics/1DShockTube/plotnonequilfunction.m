function [] = plotnonequilfunction(pde, master, mesh, dmd, str)
    dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
    for ti = pde.saveSolFreq:pde.saveSolFreq:2e5
        %     sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
        if strcmp(str,'P')
            sol = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
            u = sol(:,1,:);
            figure(1);
            ay = 1.6e4; by = 1.7e4;
        elseif strcmp(str,'T')
            sol = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
            u = sol(:,2,:);
            figure(2)
            ay = 3500; by = 6500;
        else
            sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
            rho = sum(sol(:,1:5,:),2);
            if strcmp(str,'Y1')
                u = sol(:,1,:);
                u = u./rho;
                figure(1); 
                ay = 0; by = 1.8e-3;
            end
            if strcmp(str,'Y2')
                u = sol(:,2,:);
                u = u./rho;
                figure(2); 
                ay = 0; by = 0.25;
            end
            if strcmp(str,'Y3')
                u = sol(:,3,:);
                u = u./rho;
                figure(3); 
                ay = 0; by = 8e-2;
            end
            if strcmp(str,'Y4')
                u = sol(:,4,:);
                u = u./rho;
                figure(4); 
                ay = 0.73; by = 0.78;
            end
            if strcmp(str,'Y5')
                u = sol(:,5,:);
                u = u./rho;
                figure(5); 
                ay = 0; b = 0.25;
            end
            if strcmp(str,'r1')
                u = sol(:,1,:);
                figure(1); 
            end
            if strcmp(str,'r2')
                u = sol(:,2,:);
                figure(2); 
            end
            if strcmp(str,'r3')
                u = sol(:,3,:);
                figure(3); 
            end
            if strcmp(str,'r4')
                u = sol(:,4,:);
                figure(4); 
            end
            if strcmp(str,'r5')
                u = sol(:,1,:);
                figure(5); 
            end
            if strcmp(str,'r')
                u = rho;
                figure(6); 
            end
            if strcmp(str,'u')
                u = sol(:,6,:);
                u = u./rho;
                figure(7); 
            end
            if strcmp(str,'ru')
                u = sol(:,6,:);
                figure(8); 
            end
        end
        clf; hold on
        plot(dgnodes(:),u(:),'LineWidth',1.3); 
        set(gca,'Xscale','log');
        ylim([ay by])
        xlim([1e-6 1.5])
        drawnow
%         waitforbuttonpress
    end
end
    
            
    