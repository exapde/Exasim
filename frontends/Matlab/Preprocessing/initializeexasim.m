function [pde,mesh] = initializeexasim()

version = "backend";
pde = initializepde(version);  % initialize pde struct
mesh = initializemesh(version);  % initialize mesh struct

