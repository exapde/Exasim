#include <iostream>
#include <vector>
#include <string>

#include "readbinfiles.cpp"

// ./readsolution /Users/cuongnguyen/Documents/GitHub/Exasim/apps/navierstokes/orion/dataout/outudg /Users/cuongnguyen/Documents/GitHub/Exasim/apps/navierstokes/orion/datain/mesh 4

int main(int argc, char** argv)
{
    try {
        // ./readsol dataout/outudg datain/mesh 4 
        if (argc < 4) {
            std::cerr
              << "Usage:\n  " << argv[0]
              << " <sol_base> <elempart_base> <nprocs>"              
              << " [nsteps] [stepoffsets]\n";
            return 1;
        }

        const std::string sol_base         = argv[1];
        const std::string elempart_base    = argv[2];        
        const int nprocs   = std::atoi(argv[3]);
        const int nsteps      = (argc >= 5) ? std::atoi(argv[4]) : 1;
        const int stepoffsets = (argc >= 6) ? std::atoi(argv[5]) : 0;

        if (nprocs <= 0)  throw std::runtime_error("nprocs must be positive");

        // ---------- Read elempartpts / elempart ----------
        std::vector<std::vector<int>> elempartpts(nprocs);
        std::vector<std::vector<int>> elempart(nprocs);
        readelempart(elempart_base, elempart, elempartpts, nprocs);
       
        for (int n=0; n < nprocs; n++) {
          int n1 = elempart[n].size();
          int n2 = elempartpts[n].size();
          printf("n = %d, elempart size = %d, elempartpts size = %d\n", n, n1, n2);
        }

        // ---------- Read solution ----------
        std::vector<double> sol;        

        for (int step = 0; step < nsteps; step++) {

            int n1 = 0, n2 = 0, ne = 0;

            // Keep sol's capacity across steps to reduce churn
            sol.clear();

            readsolution(sol_base,
                        elempartpts,
                        elempart,
                        sol,
                        1,
                        stepoffsets + step,
                        n1, n2, ne);
    
            const std::string out = "sol_step_" + std::to_string(stepoffsets + step) + ".bin";
            writearray2file(out, sol.data(), sol.size());
            std::cout << "Wrote " << out<< " (sol.size()=" << sol.size() << ")\n";
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}

