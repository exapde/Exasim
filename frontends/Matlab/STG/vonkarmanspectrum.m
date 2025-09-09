function Ek = vonkarmanspectrum(k, ke, keta, Eamp) 

% Energy spectrum as a function of the wave number k
kke = k/ke;
kke2 = kke.*kke;
Ek = Eamp.*((kke2 .* kke2)./((1.0 + kke2).^(17/6))).*exp(-2.0*(k/keta).^2);

