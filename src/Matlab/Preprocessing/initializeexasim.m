function [pde,mesh] = initializeexasim()

version = "src";
pde = initializepde(version);  % initialize pde struct
mesh = initializemesh(version);  % initialize mesh struct

