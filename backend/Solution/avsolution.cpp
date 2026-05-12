void avdistfunc(CSolution** pdemodel, ofstream* out, Int nummodels, Int backend)
{  
  for (int i=0; i<nummodels; i++) 
    pdemodel[i]->InitSolution(backend); 
  
  Int aviter = pdemodel[0]->disc.app.szavparam/2;
  for (Int n=0; n<aviter; n++) {
    if (pdemodel[0]->disc.common.mpiRank==0)
      printf("AV continuation iteration: %d\n", n+1);
    
    for (Int i=0; i<nummodels; i++) {
      Int m = pdemodel[i]->disc.app.szphysicsparam;      
      ArrayCopy(&pdemodel[i]->disc.app.physicsparam[m-2], &pdemodel[i]->disc.app.avparam[2*n], 2);
      pdemodel[i]->SteadyProblem(out[i], backend);        
    }
  }
  
  for (int i=0; i<nummodels; i++) {        
    pdemodel[i]->SaveSolutions(backend);    
    pdemodel[i]->SaveSolutionsOnBoundary(backend);         
    if (pdemodel[i]->disc.common.nce>0)
      pdemodel[i]->SaveOutputCG(backend);            
  }
}



