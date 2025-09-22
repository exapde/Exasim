function ttp = distribute(m,spx,spy,n)
ttp = 0:(n-1)/m:n-1;
for iter = 1:50
 xp = fnval(spx,ttp);
 yp = fnval(spy,ttp);
 dis = sqrt((xp(2:end)-xp(1:end-1)).^2+(yp(2:end)-yp(1:end-1)).^2);
 ttp(2:end-1) = ttp(2:end-1) + 0.4*(ttp(3:end)-ttp(1:end-2)).*(dis(2:end)-dis(1:end-1))./(dis(2:end)+dis(1:end-1));
end
end
