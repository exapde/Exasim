## Chen test subproblem notes
This is a living document that serves as a log of the simulation setup for the Chen ionic wind case


## Prior to 2.7.23
- Workekd on minimum example for the Poisson case -> inhomogeneous forcing function on the interior, dirichet BCs on the exterior
- Proved for both Sam's and Cuong's files
- Started implementing the subproblem but ran into errors when trying to run the coupled case
- Issue was something to do with the dimension of the vector w in the pdemodel files.   

File breakdown in the research/chen_test_subproblem directory:
- Cuong's minimal example for the poisson problem
1. pdepoi_cuong.m
2. pdemodel.m

- Sam's minimal example for the poisson problem
1. pde_poisson.m
2. pdemodel_poisson.m

- Full-fledged coupled case
1. pdeapp.m
2. pdeapp1.m
3. pdeapp2.m
4. pdemodel1.m
5. pdemodel2.m (shared with Cuong's example)


