function cmp = compilepdemodel(pde)

cdir = pwd();
cmp = pde.exasimpath + "/backend/Model/build";
cd(cmp);

!cmake ..
!make 

cd(cdir);

