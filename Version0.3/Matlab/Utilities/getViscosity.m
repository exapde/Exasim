function [mu] = getViscosity(muRef,Tref,Tdim,vModel)
%GETVISCOSITY Summary of this function goes here
%   Detailed explanation goes here

Ts = 110.4;

if vModel==0
    mu = muRef;
elseif vModel==1
    mu = muRef*(Tdim./Tref).^(3/2) .* (Tref + Ts)./(Tdim + Ts);
else
    error('Viscosity Model Unknown');
end

