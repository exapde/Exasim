#include <iostream>
#include <vector>
#include <string>

#include "readbinfiles.cpp"

static void checkMatlabIndices(const std::vector<int>& i_matlab,
                               const std::vector<int>& j_matlab,
                               int p1, int ne_z)
{
    if (i_matlab.size() != j_matlab.size()) {
        throw std::runtime_error("i_matlab and j_matlab must have the same size.");
    }
    for (std::size_t k = 0; k < i_matlab.size(); ++k) {
        if (i_matlab[k] < 1 || i_matlab[k] > p1) {
            throw std::runtime_error("i_matlab[" + std::to_string(k) + "]=" +
                                     std::to_string(i_matlab[k]) +
                                     " out of range [1," + std::to_string(p1) + "]");
        }
        if (j_matlab[k] < 1 || j_matlab[k] > ne_z) {
            throw std::runtime_error("j_matlab[" + std::to_string(k) + "]=" +
                                     std::to_string(j_matlab[k]) +
                                     " out of range [1," + std::to_string(ne_z) + "]");
        }
    }
}

int main(int argc, char** argv)
{
    try {
        // ./extract2d outudg ./elempart_np ./elempartpts_np 8 9 80 2,2,2 20,40,60 1 0
        if (argc < 9) {
            std::cerr
              << "Usage:\n  " << argv[0]
              << " <sol_base> <elempart_base> <elempartpts_base> <nprocs>"
              << " <npe2d> <ne_z> <i_matlab_csv> <j_matlab_csv>"
              << " [nsteps] [stepoffsets]\n";
            return 1;
        }

        const std::string sol_base         = argv[1];
        const std::string elempart_base    = argv[2];
        const std::string elempartpts_base = argv[3];
        const int nprocs   = std::atoi(argv[4]);
        const int npe2d    = std::atoi(argv[5]);
        const int ne_z     = std::atoi(argv[6]);

        const std::vector<int> i_matlab = parseCSVInts(argv[7]);
        const std::vector<int> j_matlab = parseCSVInts(argv[8]);

        const int nsteps      = (argc >= 10) ? std::atoi(argv[9])  : 1;
        const int stepoffsets = (argc >= 11) ? std::atoi(argv[10]) : 0;

        if (nprocs <= 0)  throw std::runtime_error("nprocs must be positive");
        if (npe2d <= 0)   throw std::runtime_error("npe2d must be positive");
        if (ne_z <= 0)    throw std::runtime_error("ne_z must be positive");
        if (i_matlab.empty())
            throw std::runtime_error("i_matlab cannot be empty");
        if (i_matlab.size() != j_matlab.size())
            throw std::runtime_error("i_matlab and j_matlab must have same size");

        // ---------- Read elempartpts / elempart ----------
        std::vector<std::vector<int>> elempartpts(nprocs);
        std::vector<std::vector<int>> elempart(nprocs);

        for (int r = 0; r < nprocs; ++r) {
            elempartpts[r] = readBinInts(elempartpts_base + std::to_string(r) + ".bin");
            elempart[r]    = readBinInts(elempart_base    + std::to_string(r) + ".bin");
        }

        // ---------- Read solution ----------
        std::vector<double> sol3dGlobal;        

        for (int step = 0; step < nsteps; step++) {

            int n1 = 0, n2 = 0, ne = 0;

            // Keep sol3dGlobal's capacity across steps to reduce churn
            sol3dGlobal.clear();

            readsol_cpp(sol_base,
                        elempartpts,
                        elempart,
                        sol3dGlobal,
                        1,
                        stepoffsets + step,
                        n1, n2, ne);
    
            // ---------- Derived dimensions ----------
            const int nc  = n2;
            const int p1  = n1 / npe2d;
            const int ne2 = ne / ne_z;
    
            if (n1 % npe2d != 0)
                throw std::runtime_error("n1 not divisible by npe2d");
            if (ne % ne_z != 0)
                throw std::runtime_error("ne not divisible by ne_z");
    
            // Validate MATLAB indices against derived bounds
            checkMatlabIndices(i_matlab, j_matlab, p1, ne_z);

            // ---------- Extract 2D ----------
            std::vector<double> sol2d =
                extractSol2D(sol3dGlobal,
                             npe2d, p1, nc, ne2, ne_z,
                             i_matlab, j_matlab);
    
            const std::string out =
                "sol2d_step_" + std::to_string(stepoffsets + step) + ".bin";

            writearray2file(out, sol2d.data(), sol2d.size());

            std::cout << "Wrote " << out
                      << " (sol2d.size()=" << sol2d.size() << ")\n";
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}

