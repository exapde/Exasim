fprintf('  --> Initializing EXASIM ...')
d0=fileparts([pwd,filesep]);
addpath([d0,'/matlab/codegen']);
addpath([d0,'/matlab/mesh/distmesh']);
addpath([d0,'/matlab/mesh/mkmesh']);
addpath([d0,'/matlab/mesh/cmesh']);
addpath([d0,'/matlab/mesh/foilmesh']);
addpath([d0,'/matlab/mesh/lesmesh']);
addpath([d0,'/matlab/mesh/airfoilTools']);
addpath([d0,'/matlab/mesh']);
addpath([d0,'/matlab/master']);
addpath([d0,'/matlab/util']);
addpath([d0,'/matlab/plot']);
addpath([d0,'/matlab/preprocessing']);
fprintf(' \n Done. \n\n');
clear d0

