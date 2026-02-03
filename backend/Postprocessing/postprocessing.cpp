#include <iostream>
#include <vector>
#include <string>

#include "readbinfiles.cpp"

// ------------------------- main -------------------------
int main(int argc, char** argv)
{
    try {
        // Command line:
        //   argv[1] = sol_base          (solution file prefix, without "_npX.bin")
        //   argv[2] = elempart_base     (prefix for elempart rank files, e.g. ".../elempart_np")
        //   argv[3] = elempartpts_base  (prefix for elempartpts rank files, e.g. ".../elempartpts_np")
        //   argv[4] = nprocs
        //   argv[5] = nsteps (optional, default 1)
        //   argv[6] = stepoffsets (optional, default 0)
        //
        // Example:
        //   ./reader outsol ./elempart_np ./elempartpts_np 8 1 0
        if (argc < 5) {
            std::cerr
                << "Usage:\n  " << argv[0]
                << " <sol_base> <elempart_base> <elempartpts_base> <nprocs> [nsteps] [stepoffsets]\n\n"
                << "Where files are assumed as:\n"
                << "  Solution:     sol_base + \"_np\" + rank + \".bin\"\n"
                << "  elempart:     elempart_base + rank + \".bin\"\n"
                << "  elempartpts:  elempartpts_base + rank + \".bin\"\n";
            return 1;
        }

        const std::string sol_base        = argv[1];
        const std::string elempart_base   = argv[2];
        const std::string elempartpts_base= argv[3];
        const int nprocs = std::atoi(argv[4]);

        const int nsteps      = (argc >= 6) ? std::atoi(argv[5]) : 1;
        const int stepoffsets = (argc >= 7) ? std::atoi(argv[6]) : 0;

        if (nprocs <= 0) throw std::runtime_error("nprocs must be positive.");
        if (nsteps <= 0) throw std::runtime_error("nsteps must be positive.");
        if (stepoffsets < 0) throw std::runtime_error("stepoffsets must be nonnegative.");

        // --------- Read elempartpts and elempart from rank files ---------
        std::vector<std::vector<int>> elempartpts(nprocs);
        std::vector<std::vector<int>> elempart(nprocs);

        for (int r = 0; r < nprocs; ++r) {
            const std::string f_pts = elempartpts_base + std::to_string(r) + ".bin";
            const std::string f_ep  = elempart_base    + std::to_string(r) + ".bin";

            elempartpts[r] = readBinInts(f_pts);
            elempart[r]    = readBinInts(f_ep);

            if (elempartpts[r].size() < 2) {
                throw std::runtime_error("elempartpts[" + std::to_string(r) + "] has <2 ints in file: " + f_pts);
            }
            if (elempart[r].empty()) {
                throw std::runtime_error("elempart[" + std::to_string(r) + "] is empty in file: " + f_ep);
            }
        }

        // --------- Read solution and assemble only sol3dGlobal ---------
        std::vector<double> sol3dGlobal;
        int n1 = 0, n2 = 0, ne = 0;

        readsol_cpp(sol_base,
                    elempartpts,
                    elempart,
                    sol3dGlobal,
                    nsteps,
                    stepoffsets,
                    n1, n2, ne);

        // --------- Report ---------
        std::cout << "Read sol3dGlobal successfully.\n";
        std::cout << "n1=" << n1 << " n2=" << n2 << " ne=" << ne << " nsteps=" << nsteps << "\n";
        std::cout << "sol3dGlobal.size() = " << sol3dGlobal.size() << "\n";

        // Quick sanity check: print first entry of step 0, elem 0, (0,0)
        if (!sol3dGlobal.empty()) {
            //writearray2file("sol.bin", sol3dGlobal.data(), sol3dGlobal.size());

            int nc = n2;
            int p1 = n1/npe2d;
            int ne2 = ne/ne_z;
            std::vector<double> sol2d = extractSol2D(sol3dGlobal, npe2d, p1, nc, 
                                          ne2, ne_z, i_matlab, j_matlab);

            writearray2file("sol2d.bin", sol2d.data(), sol2d.size());
          
            // column-major index for sol(n1,n2,ne,nsteps):
            auto idx4 = [&](int i1, int i2, int e, int s) -> std::size_t {
                return static_cast<std::size_t>(i1)
                     + static_cast<std::size_t>(n1) * (
                       static_cast<std::size_t>(i2)
                     + static_cast<std::size_t>(n2) * (
                       static_cast<std::size_t>(e)
                     + static_cast<std::size_t>(ne) * static_cast<std::size_t>(s)));
            };

            std::cout << "sol(0,0,0,0) = " << sol3dGlobal[idx4(0,0,0,0)] << "\n";
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}


// // Declaration of readsol_cpp (must match the definition exactly)
// void readsol_cpp(const std::string& base,
//                  const std::vector<std::vector<int>>& elempartpts,
//                  const std::vector<std::vector<int>>& elempart,
//                  std::vector<double>& sol3dGlobal,   // OUT
//                  int nsteps,
//                  int stepoffsets,
//                  int& n1_out,                         // OUT
//                  int& n2_out,                         // OUT
//                  int& ne_out);                        // OUT
// 
// int main(int argc, char** argv)
// {
//     try {
//         // ---------------- User inputs ----------------
//         // Base name of MPI solution files:
//         //   base_np0.bin, base_np1.bin, ...
//         std::string base = "outsol";   // example
// 
//         int nsteps = 1;
//         int stepoffsets = 0;
// 
//         // ---------------- Read elempart data ----------------
//         // In practice, read these from binary files.
//         // Here we just illustrate the structure.
// 
//         // Number of MPI ranks
//         int nprocs = 2;
// 
//         // elempartpts{r} = [nInterior, nBoundary, ...]
//         std::vector<std::vector<int>> elempartpts(nprocs);
//         elempartpts[0] = {10, 0};
//         elempartpts[1] = {12, 0};
// 
//         // elempart{r} = global element indices (1-based, MATLAB style)
//         std::vector<std::vector<int>> elempart(nprocs);
//         elempart[0].resize(10);
//         elempart[1].resize(12);
// 
//         // Example: fill with dummy global indices
//         for (int i = 0; i < 10; ++i) elempart[0][i] = i + 1;
//         for (int i = 0; i < 12; ++i) elempart[1][i] = 10 + i + 1;
// 
//         // ---------------- Call readsol_cpp ----------------
//         std::vector<double> sol3dGlobal;
//         int n1 = 0, n2 = 0, ne = 0;
// 
//         readsol_cpp(base,
//                     elempartpts,
//                     elempart,
//                     sol3dGlobal,
//                     nsteps,
//                     stepoffsets,
//                     n1,
//                     n2,
//                     ne);
// 
//         // ---------------- Post-processing ----------------
//         std::cout << "Read solution successfully\n";
//         std::cout << "n1 = " << n1 << ", n2 = " << n2 << ", ne = " << ne << "\n";
//         std::cout << "Total entries in sol3dGlobal = "
//                   << sol3dGlobal.size() << "\n";
// 
//         // Access example: sol(i1,i2,elem,step)
//         auto idx4 = [&](int i1, int i2, int e, int s) -> std::size_t {
//             return static_cast<std::size_t>(i1)
//                  + static_cast<std::size_t>(n1) * (
//                    static_cast<std::size_t>(i2)
//                  + static_cast<std::size_t>(n2) * (
//                    static_cast<std::size_t>(e)
//                  + static_cast<std::size_t>(ne) * static_cast<std::size_t>(s)));
//         };
// 
//         // Print one value as a sanity check
//         if (!sol3dGlobal.empty()) {
//             std::cout << "sol(0,0,0,0) = "
//                       << sol3dGlobal[idx4(0,0,0,0)] << "\n";
//         }
// 
//         return 0;
//     }
//     catch (const std::exception& e) {
//         std::cerr << "ERROR: " << e.what() << std::endl;
//         return 1;
//     }
// }
