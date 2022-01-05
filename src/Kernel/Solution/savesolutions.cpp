#ifndef __SAVESOLUTIONS
#define __SAVESOLUTIONS

void SaveSolutions(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{           
   if (common.tdep==1) { 
        if (((common.currentstep+1) % common.saveSolFreq) == 0)             
        {        
            string filename = common.fileout + "_t" + NumberToString(common.currentstep+1) + "_np" + NumberToString(common.mpiRank) + ".bin";     
            if (common.saveSolOpt==0)
                writearray2file(filename, sys.u, common.ndof1, backend);
            else
                writearray2file(filename, sol.udg, common.ndofudg1, backend);
            
            if (common.wave==1) {
                string fn = common.fileout + "_wdg_t" + NumberToString(common.currentstep+1) + "_np" + NumberToString(common.mpiRank) + ".bin";                    
                writearray2file(fn, sys.wtmp, common.ndof1, backend);
            }
            
            
        }                                
   }
   else {
        string filename = common.fileout + "_np" + NumberToString(common.mpiRank) + ".bin";                    
        if (common.saveSolOpt==0)
            writearray2file(filename, sys.u, common.ndof1, backend);
        else
            writearray2file(filename, sol.udg, common.ndofudg1, backend);       
   }
}

#endif


