function [kcute, lemax] = calculatek(hx, hyz, dwall, lt, nu, eps)

%hyz = max(hy(:),hz(:));
hmax = max(hx, hyz);

% spectral maximum wavenumber
le = min(2.0 * dwall(:), 3.0 * lt(:));
lemax = max(le);
kcute(:,1) = reshape(2.0*(pi./le), size(dwall));

% cut-off wavenumber
lcut = 2.0 * min( max(0.3 * hmax, hyz) + 0.1 * dwall, hmax);
kcute(:,2) = reshape(2.0*(pi./lcut), size(dwall));

% Kolmogorov wavenumber
leta = (nu.*nu.*nu./eps).^(0.25);
kcute(:,3) = 2*pi./leta;


