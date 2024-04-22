function [x,y] = read_foil(name)
%ht13foil - sd7003foil
[x,y]=textread(name, '%f %f','headerlines',1); 
x(end)=x(1);
y(end)=y(1);