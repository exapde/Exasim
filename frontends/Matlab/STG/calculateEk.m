function Ek = calculateEk(k, ke, kcut, keta) 

% calculate feta
feta = exp(-1.0 * (12.0*k./keta) .* (12.0*k./keta));

% calculate fcut
base = 4.0 * max(k - 0.9*kcut, 0.0)./kcut;
fcut = exp( -1.0 * base .* base .* base);

% calculate Ek
kke = k./ke;
kke2 = kke.*kke;
Ek = feta.*fcut.*(kke2 .* kke2)./((1.0 + 2.4 * kke2).^(17.0/6.0));



