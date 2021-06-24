function gencppcode(appname,xdg,udg,odg,uhg,nlg,tau,uinf,param,time,...
    flux,source,tdfunc,ubou,fbou,uhat,fhat,avfd)

gencppflux(flux,xdg,udg,odg,param,time,appname);
gencppfbou(fbou,xdg,udg,uhg,odg,nlg,tau,uinf,param,time,appname);
gencppubou(ubou,xdg,udg,uhg,odg,nlg,tau,uinf,param,time,appname);
gencppsource(source,xdg,udg,odg,param,time,appname);
gencpptdfunc(tdfunc,xdg,udg,odg,param,time,appname);
gencppavfield(avfd,xdg,udg,odg,param,time,appname);

% if ~isempty(uhat)
%     gencppuhat(uhat,xdg,udg,uhg,odg,nlg,tau,uinf,param,time,appname);
% end
% if ~isempty(fhat)
%     gencppfhat(fhat,xdg,udg,uhg,odg,nlg,tau,uinf,param,time,appname);
% end

