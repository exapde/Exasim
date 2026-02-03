#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstddef>
#include <sstream>
#include <iostream>
#include <iomanip>   // std::setw, std::setprecision
#include <array>     // std::array
#include <cstdint>   // std::int64_t

static inline void PrintErrorAndExit(const std::string& errmsg, const char* file, int line)
{
    int rank = 0;

#ifdef HAVE_MPI
    MPI_Comm_rank(EXASIM_COMM_WORLD, &rank);
#endif
    
    fprintf(stderr,
            "\n==============================================\n"
            "[Rank %d] ERROR: %s\n"
            "  Location: %s:%d\n"
            "==============================================\n\n",
            rank, errmsg.c_str(), file, line);
    fflush(stderr);

#ifdef HAVE_MPI
    // Abort the entire MPI job instead of trying to finalize gracefully.
    // MPI_Finalize() is unsafe after a runtime error and can hang.
    MPI_Abort(EXASIM_COMM_WORLD, EXIT_FAILURE);
#else
    exit(EXIT_FAILURE);
#endif
}

// -----------------------------------------------------------------------------
// Macro for convenient error reporting
// -----------------------------------------------------------------------------
#define error(msg)  PrintErrorAndExit((msg), __FILE__, __LINE__)

template <typename T> std::string NumberToString ( T Number )
{
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}

std::vector<int> parseCSVInts(const std::string& s)
{
    std::vector<int> v;
    std::size_t start = 0;
    while (start < s.size()) {
        std::size_t end = s.find(',', start);
        if (end == std::string::npos) end = s.size();
        v.push_back(std::stoi(s.substr(start, end - start)));
        start = end + 1;
    }
    return v;
}

void print2iarray(const int* a,
                  int m, int n,
                  const std::string& msg)
{
    std::cout << "----------------------------------------\n";
    std::cout << msg << "\n";
    std::cout << "Dimensions: (" << m << " x " << n << ")\n";
  
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << std::setw(12) << a[i + j * m] << " ";
        }
        std::cout << "\n";
    }
}

void print2darray(const double* a,
                  int m, int n,
                  const std::string& msg)
{
    std::cout << "----------------------------------------\n";
    std::cout << msg << "\n";
    std::cout << "Dimensions: (" << m << " x " << n << ")\n";

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // column-major access
            std::cout << std::setw(12)
                      << std::setprecision(6)
                      << std::scientific
                      << a[i + j * m] << " ";
        }
        std::cout << "\n";
    }
}

struct ReadoutResult {
    int npe = 0;
    int ncu = 0;   // = nc in your MATLAB code
    int Ne  = 0;   // total number of elements
    // sol stored as column-major array with layout [npe][ncu][Ne]
    // index = p + npe*(cu + ncu*e)
    std::vector<double> sol;
};

static std::vector<double> readBinDoubles(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    in.seekg(0, std::ios::end);
    std::streamoff bytes = in.tellg();
    in.seekg(0, std::ios::beg);

    if (bytes < 0 || (bytes % static_cast<std::streamoff>(sizeof(double)) != 0)) {
        throw std::runtime_error("File size is not a multiple of sizeof(double): " + filename);
    }

    std::size_t n = static_cast<std::size_t>(bytes / sizeof(double));
    std::vector<double> data(n);
    in.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(bytes));
    if (!in) {
        throw std::runtime_error("Failed reading file: " + filename);
    }

    return data;
}

static std::vector<int> readBinInts(const std::string& filename)
{    
    std::vector<double> tmp = readBinDoubles(filename);

    std::size_t n = tmp.size();
    std::vector<int> data(n);    
    for (int j = 0; j < n; j++)
      data[j] = (int) tmp[j];

    return data;
}

template <typename T> void writearray2file(std::string filename, T *a, int N)
{
    if (N>0) {
        // Open file to read
        std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }

        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );

        out.close();
    }
}

// elempartpts: size nprocs; each entry must have >=2 integers, using only first two, summed.
// elempart:    size nprocs; each entry length >= Nj, contains GLOBAL element indices.
//              NOTE: MATLAB indices are 1-based; we convert to 0-based here.
// base: prefix (e.g. "/path/outudg_np")
// npe: nodes-per-element
ReadoutResult readout_cpp(const std::string& base,
                          const std::vector<std::vector<int>>& elempartpts,
                          const std::vector<std::vector<int>>& elempart,
                          int npe)
{
    const int nprocs = static_cast<int>(elempartpts.size());
    if (static_cast<int>(elempart.size()) != nprocs) {
        throw std::runtime_error("elempart size must match elempartpts size.");
    }
    if (npe <= 0) {
        throw std::runtime_error("npe must be positive.");
    }

    int Ne = 0;
    std::vector<int> Ni(nprocs, 0);

    for (int i = 0; i < nprocs; ++i) {
        if (static_cast<int>(elempartpts[i].size()) < 2) {
            throw std::runtime_error("elempartpts[" + std::to_string(i) + "] must have at least 2 entries.");
        }
        int Nj = elempartpts[i][0] + elempartpts[i][1];  // sum(elempartpts{i}(1:2))
        Ne += Nj;
        Ni[i] = Nj;
    }

    // Read proc 0 to infer nc
    const std::string filesol0 = base + std::to_string(0) + ".bin";
    std::vector<double> udg0 = readBinDoubles(filesol0);

    const int N1 = Ni[0];
    if (N1 <= 0) {
        throw std::runtime_error("Ni(1) is non-positive; cannot infer nc.");
    }

    const std::size_t denom = static_cast<std::size_t>(npe) * static_cast<std::size_t>(N1);
    if (denom == 0 || (udg0.size() % denom) != 0) {
        throw std::runtime_error("udg0 size not divisible by (npe*N1); cannot compute nc.");
    }

    const int nc  = static_cast<int>(udg0.size() / denom);
    const int ncu = nc;

    ReadoutResult out;
    out.npe = npe;
    out.ncu = ncu;
    out.Ne  = Ne;
    out.sol.assign(static_cast<std::size_t>(npe) * ncu * Ne, 0.0);

    // Helper for sol indexing: sol(p,cu,e) with column-major [npe][ncu][Ne]
    auto solIndex = [=](int p, int cu, int e) -> std::size_t {
        return static_cast<std::size_t>(p)
             + static_cast<std::size_t>(npe) * (static_cast<std::size_t>(cu)
             + static_cast<std::size_t>(ncu) * static_cast<std::size_t>(e));
    };

    // MATLAB reshape(udg, npe, nc, Nj) is column-major:
    // tm(p,cu,j) corresponds to udg[ p + npe*(cu + nc*j) ] with 0-based p,cu,j.
    auto udgIndex = [=](int p, int cu, int j) -> std::size_t {
        return static_cast<std::size_t>(p)
             + static_cast<std::size_t>(npe) * (static_cast<std::size_t>(cu)
             + static_cast<std::size_t>(nc) * static_cast<std::size_t>(j));
    };

    for (int i = 0; i < nprocs; ++i) {
        const int Nj = Ni[i];
        const std::string filesol = base + std::to_string(i) + ".bin"; // i-1 in MATLAB, but i starts at 0 here
        std::vector<double> udg = readBinDoubles(filesol);

        const std::size_t expected = static_cast<std::size_t>(npe) * nc * static_cast<std::size_t>(Nj);
        if (udg.size() != expected) {
            std::cerr << "Size mismatch at proc " << i
                      << " npe=" << npe << " nc=" << nc << " Nj=" << Nj
                      << " udg.size=" << udg.size() << " expected=" << expected << "\n";
            continue;
        }

        if (static_cast<int>(elempart[i].size()) < Nj) {
            throw std::runtime_error("elempart[" + std::to_string(i) + "] shorter than Nj.");
        }

        for (int j = 0; j < Nj; ++j) {
            // Convert MATLAB 1-based global element index -> 0-based
            const int eGlobal0 = elempart[i][j] - 1;
            if (eGlobal0 < 0 || eGlobal0 >= Ne) {
                throw std::runtime_error("Global element index out of range in elempart at proc "
                                         + std::to_string(i) + ", j=" + std::to_string(j));
            }

            for (int cu = 0; cu < ncu; ++cu) {
                for (int p = 0; p < npe; ++p) {
                    out.sol[solIndex(p, cu, eGlobal0)] = udg[udgIndex(p, cu, j)];
                }
            }
        }
    }

    return out;
}

std::vector<std::vector<int>>
readElempart(const std::string& baseFilename, int nfiles)
{
    if (nfiles <= 0) {
        throw std::runtime_error("Number of files must be positive.");
    }

    std::vector<std::vector<int>> elempart(nfiles);

    for (int i = 0; i < nfiles; ++i) {
        std::string filename = baseFilename + std::to_string(i) + ".bin";
        elempart[i] = readBinInts(filename);
    }

    return elempart;
}

struct ReadSolMPIResult
{
    // k[rank] = {n1, n2, n3, timesteps}
    std::vector<std::array<std::int64_t, 4>> k;

    // If readSolution == true, sol is allocated and filled.
    // Stored column-major as sol(n1,n2,ne,nsteps).
    std::vector<double> sol;

    // Global dimensions for sol (meaningful if sol is filled)
    std::int64_t n1 = 0, n2 = 0, ne = 0, nsteps = 0;
};

// Read exactly `count` doubles from stream; throws on failure.
static void readDoubles(std::ifstream& in, double* dst, std::size_t count, const std::string& fname)
{
    in.read(reinterpret_cast<char*>(dst),
            static_cast<std::streamsize>(count * sizeof(double)));
    if (!in) {
        throw std::runtime_error("Failed to read doubles from file: " + fname);
    }
}

static std::int64_t fileSizeBytes(const std::string& fname)
{
    std::ifstream in(fname, std::ios::binary | std::ios::ate);
    if (!in) {
        throw std::runtime_error("Cannot open file: " + fname);
    }
    auto sz = in.tellg();
    if (sz < 0) {
        throw std::runtime_error("Failed to get file size for: " + fname);
    }
    return static_cast<std::int64_t>(sz);
}

// Main function: reads headers for all procs; optionally reads/assembles solution.
ReadSolMPIResult readsolmpi_cpp(const std::string& base,
                               int nprocs,
                               int nsteps = 1,
                               int stepoffsets = 0,
                               bool readSolution = true)
{
    if (nprocs <= 0) {
        throw std::runtime_error("nprocs must be positive.");
    }
    if (nsteps <= 0) {
        throw std::runtime_error("nsteps must be positive.");
    }
    if (stepoffsets < 0) {
        throw std::runtime_error("stepoffsets must be nonnegative.");
    }

    ReadSolMPIResult out;
    out.k.resize(static_cast<std::size_t>(nprocs));

    // ---------- Pass 1: read headers + compute timesteps per rank ----------
    for (int r = 0; r < nprocs; ++r) {
        const std::string fname = base + "_np" + std::to_string(r) + ".bin";

        std::ifstream in(fname, std::ios::binary);
        if (!in) {
            throw std::runtime_error("Cannot open file: " + fname);
        }

        double hdr[3] = {0, 0, 0};
        readDoubles(in, hdr, 3, fname);

        const std::int64_t n1 = static_cast<std::int64_t>(hdr[0]);
        const std::int64_t n2 = static_cast<std::int64_t>(hdr[1]);
        const std::int64_t n3 = static_cast<std::int64_t>(hdr[2]);

        if (n1 <= 0 || n2 <= 0 || n3 <= 0) {
            throw std::runtime_error("Invalid header (nonpositive dims) in: " + fname);
        }

        const std::int64_t N = n1 * n2 * n3;

        const std::int64_t bytes = fileSizeBytes(fname);
        if (bytes % static_cast<std::int64_t>(sizeof(double)) != 0) {
            throw std::runtime_error("File size is not multiple of 8 bytes (double): " + fname);
        }

        const std::int64_t L = bytes / static_cast<std::int64_t>(sizeof(double)); // number of doubles
        if (L < 3) {
            throw std::runtime_error("File too small to contain header: " + fname);
        }
        if ((L - 3) % N != 0) {
            throw std::runtime_error("File payload not divisible by N; corrupted or wrong format: " + fname);
        }

        const std::int64_t timesteps = (L - 3) / N;
        out.k[static_cast<std::size_t>(r)] = {n1, n2, n3, timesteps};
    }

    // If user only wants k, stop here.
    if (!readSolution) {
        return out;
    }

    // ---------- Determine global sizes ----------
    // MATLAB uses n1,n2 from header (assumed consistent); ne = sum(n3)
    const std::int64_t n1g = out.k[0][0];
    const std::int64_t n2g = out.k[0][1];

    std::int64_t ne = 0;
    for (int r = 0; r < nprocs; ++r) {
        // Optionally enforce consistency like the commented MATLAB checks:
        // if (out.k[r][0] != n1g) throw ...
        // if (out.k[r][1] != n2g) throw ...
        ne += out.k[static_cast<std::size_t>(r)][2];
    }

    out.n1 = n1g;
    out.n2 = n2g;
    out.ne = ne;
    out.nsteps = nsteps;

    // Allocate sol(n1,n2,ne,nsteps) column-major
    const std::size_t total = static_cast<std::size_t>(n1g) *
                              static_cast<std::size_t>(n2g) *
                              static_cast<std::size_t>(ne) *
                              static_cast<std::size_t>(nsteps);
    out.sol.assign(total, 0.0);

    auto solIndex = [&](std::int64_t i1, std::int64_t i2, std::int64_t ie, std::int64_t istep) -> std::size_t {
        // i1 in [0,n1), i2 in [0,n2), ie in [0,ne), istep in [0,nsteps)
        return static_cast<std::size_t>(i1)
             + static_cast<std::size_t>(n1g) * (static_cast<std::size_t>(i2)
             + static_cast<std::size_t>(n2g) * (static_cast<std::size_t>(ie)
             + static_cast<std::size_t>(ne)  * static_cast<std::size_t>(istep)));
    };

    // ---------- Pass 2: read and assemble ----------
    std::int64_t mn = 0; // MATLAB mn tracks element offset (0-based here)
    for (int r = 0; r < nprocs; ++r) {
        const std::string fname = base + "_np" + std::to_string(r) + ".bin";
        std::ifstream in(fname, std::ios::binary);
        if (!in) {
            throw std::runtime_error("Cannot open file: " + fname);
        }

        // Read header again
        double hdr[3] = {0, 0, 0};
        readDoubles(in, hdr, 3, fname);

        const std::int64_t n1 = static_cast<std::int64_t>(hdr[0]);
        const std::int64_t n2 = static_cast<std::int64_t>(hdr[1]);
        const std::int64_t n3 = static_cast<std::int64_t>(hdr[2]);
        const std::int64_t N  = n1 * n2 * n3;

        if (n1 != n1g || n2 != n2g) {
            throw std::runtime_error("n1/n2 mismatch across ranks in file: " + fname);
        }

        // MATLAB: if (stepoffsets > 0) fseek(fid, stepoffsets*N*8, 'cof');
        if (stepoffsets > 0) {
            const std::int64_t skipBytes = static_cast<std::int64_t>(stepoffsets) * N * 8;
            in.seekg(skipBytes, std::ios::cur);
            if (!in) {
                throw std::runtime_error("seekg failed (stepoffsets too large?) in file: " + fname);
            }
        }

        std::vector<double> tm(static_cast<std::size_t>(N));

        // Each step: read N doubles, place into sol(:,:,mn:mn+n3-1, step)
        for (int s = 0; s < nsteps; ++s) {
            readDoubles(in, tm.data(), static_cast<std::size_t>(N), fname);

            // tm corresponds to reshape(tm,[n1,n2,n3]) in column-major
            // tm(i1,i2,i3) index = i1 + n1*(i2 + n2*i3)
            for (std::int64_t i3 = 0; i3 < n3; ++i3) {
                const std::int64_t ieGlobal = mn + i3; // 0-based global element slice
                for (std::int64_t i2 = 0; i2 < n2; ++i2) {
                    for (std::int64_t i1 = 0; i1 < n1; ++i1) {
                        const std::size_t tmIdx =
                            static_cast<std::size_t>(i1)
                          + static_cast<std::size_t>(n1) * (static_cast<std::size_t>(i2)
                          + static_cast<std::size_t>(n2) * static_cast<std::size_t>(i3));

                        out.sol[solIndex(i1, i2, ieGlobal, s)] = tm[tmIdx];
                    }
                }
            }
        }

        mn += n3;
    }

    return out;
}

// sol3d and sol3dnew are column-major arrays of size n1*n2*Ne_total
void assembleSol3D(
    const std::vector<std::vector<int>>& elempartpts,
    const std::vector<std::vector<int>>& elempart,
    const std::vector<double>& sol3d,
    std::vector<double>& sol3dnew,
    int n1,
    int n2)
{
    const int nprocs = static_cast<int>(elempartpts.size());

    auto idx = [&](int i1, int i2, int i3) -> std::size_t {
        return static_cast<std::size_t>(i1)
             + static_cast<std::size_t>(n1) * (static_cast<std::size_t>(i2)
             + static_cast<std::size_t>(n2) * static_cast<std::size_t>(i3));
    };

    int Ne = 0;  // MATLAB Ne starts at 0 here (0-based logic)

    for (int j = 0; j < nprocs; ++j) {

        // Nj = sum(elempartpts{j}(1:2))
        const int Nj = elempartpts[j][0] + elempartpts[j][1];

        // in = (Ne+1):(Ne+Nj)  --> 0-based: Ne ... Ne+Nj-1
        for (int k = 0; k < Nj; ++k) {

            const int srcElem = Ne + k;              // sol3d(:,:,in)
            const int dstElem = elempart[j][k] - 1;  // convert 1-based → 0-based

            for (int i2 = 0; i2 < n2; ++i2) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    sol3dnew[idx(i1, i2, dstElem)] =
                        sol3d[idx(i1, i2, srcElem)];
                }
            }
        }

        Ne += Nj;
    }
}

// Original flat sol3dnew storage (same buffer, reinterpreted)
inline std::size_t idx5(
    int a, int b, int c, int d, int e,
    int npe2, int p1, int nc, int ne2)
{
    // sol3dnew(a,b,c,d,e)
    return static_cast<std::size_t>(a)
         + static_cast<std::size_t>(npe2) * (
           static_cast<std::size_t>(b)
         + static_cast<std::size_t>(p1) * (
           static_cast<std::size_t>(c)
         + static_cast<std::size_t>(nc) * (
           static_cast<std::size_t>(d)
         + static_cast<std::size_t>(ne2) * static_cast<std::size_t>(e))));
}

inline std::size_t idx3(
    int a, int b, int c,
    int npe2, int nc)
{
    // sol2dnew(a,b,c)
    return static_cast<std::size_t>(a)
         + static_cast<std::size_t>(npe2) * (
           static_cast<std::size_t>(b)
         + static_cast<std::size_t>(nc) * static_cast<std::size_t>(c));
}

std::vector<double> extractSol2D(
    const std::vector<double>& sol3dnew_flat,
    int npe2, int p1, int nc, int ne2, int ne_z,
    int i_matlab, int j_matlab)
{
    // Convert MATLAB → C++ indexing
    const int i = i_matlab - 1;
    const int j = j_matlab - 1;

    std::vector<double> sol2dnew(
        static_cast<std::size_t>(npe2) * nc * ne2);

    for (int c = 0; c < ne2; ++c) {
        for (int b = 0; b < nc; ++b) {
            for (int a = 0; a < npe2; ++a) {

                sol2dnew[idx3(a, b, c, npe2, nc)] =
                    sol3dnew_flat[idx5(
                        a, i, b, c, j,
                        npe2, p1, nc, ne2)];
            }
        }
    }
    return sol2dnew;
}

inline std::size_t idx4(int a, int b, int c, int d,
                        int n1, int n2, int n3) {
    // a + n1*(b + n2*(c + n3*d))
    return static_cast<std::size_t>(a)
         + static_cast<std::size_t>(n1) * (static_cast<std::size_t>(b)
         + static_cast<std::size_t>(n2) * (static_cast<std::size_t>(c)
         + static_cast<std::size_t>(n3) * static_cast<std::size_t>(d)));
}

// Extended version: vectorized i/j (MATLAB 1-based)
std::vector<double> extractSol2D(
    const std::vector<double>& sol3dnew_flat,
    int npe2, int p1, int nc, int ne2, int ne_z,
    const std::vector<int>& i_matlab,
    const std::vector<int>& j_matlab)
{
    if (i_matlab.size() != j_matlab.size()) {
        throw std::runtime_error("extractSol2D: i_matlab and j_matlab must have the same length.");
    }
    const int nij = static_cast<int>(i_matlab.size());

    // Basic range checks (MATLAB indexing is 1..p1 and 1..ne_z)
    for (int k = 0; k < nij; ++k) {
        if (i_matlab[k] < 1 || i_matlab[k] > p1) {
            throw std::runtime_error("extractSol2D: i_matlab out of range at entry " + std::to_string(k));
        }
        if (j_matlab[k] < 1 || j_matlab[k] > ne_z) {
            throw std::runtime_error("extractSol2D: j_matlab out of range at entry " + std::to_string(k));
        }
    }

    // Optional: sanity check sol3dnew_flat size if you know exact dims:
    // expected = npe2 * p1 * nc * ne2 * ne_z
    const std::size_t expected =
        static_cast<std::size_t>(npe2) *
        static_cast<std::size_t>(p1) *
        static_cast<std::size_t>(nc) *
        static_cast<std::size_t>(ne2) *
        static_cast<std::size_t>(ne_z);
    if (sol3dnew_flat.size() != expected) {
        throw std::runtime_error("extractSol2D: sol3dnew_flat has unexpected size.");
    }

    // Output packs nij blocks: (npe2, nc, ne2, nij)
    std::vector<double> sol2dnew(
        static_cast<std::size_t>(npe2) * nc * ne2 * nij);

    for (int k = 0; k < nij; ++k) {
        const int i = i_matlab[k] - 1;  // MATLAB -> C++
        const int j = j_matlab[k] - 1;

        for (int c = 0; c < ne2; ++c) {
            for (int b = 0; b < nc; ++b) {
                for (int a = 0; a < npe2; ++a) {
                    sol2dnew[idx4(a, b, c, k, npe2, nc, ne2)] =
                        sol3dnew_flat[idx5(a, i, b, c, j, npe2, p1, nc, ne2)];
                }
            }
        }
    }

    return sol2dnew;
}


// // ---- helpers ----
// static void readDoubles(std::ifstream& in, double* dst, std::size_t count, const std::string& fname)
// {
//     in.read(reinterpret_cast<char*>(dst),
//             static_cast<std::streamsize>(count * sizeof(double)));
//     if (!in) throw std::runtime_error("Failed to read doubles from: " + fname);
// }

// static std::int64_t fileSizeBytes(const std::string& fname)
// {
//     std::ifstream in(fname, std::ios::binary | std::ios::ate);
//     if (!in) throw std::runtime_error("Cannot open file for size: " + fname);
//     auto sz = in.tellg();
//     if (sz < 0) throw std::runtime_error("Failed to get file size: " + fname);
//     return static_cast<std::int64_t>(sz);
// }

// Column-major slice offset for A(n1,n2,ne,nsteps):
// element-slice contiguous block = n1*n2
static inline std::size_t sliceOffset4(int n1, int n2, int ne, int elem, int step)
{
    const std::size_t slice = static_cast<std::size_t>(n1) * static_cast<std::size_t>(n2);
    return slice * (static_cast<std::size_t>(elem) + static_cast<std::size_t>(ne) * static_cast<std::size_t>(step));
}

// Reads base_np{r}.bin for r=0..nprocs-1 and assembles directly into global element order using elempart.
// Output:
//   - sol3dGlobal: column-major array representing sol(n1,n2,ne,nsteps) flattened
//   - n1_out, n2_out, ne_out set from headers
void readsol_cpp(const std::string& base,
                 const std::vector<std::vector<int>>& elempartpts,
                 const std::vector<std::vector<int>>& elempart,
                 std::vector<double>& sol3dGlobal,   // OUT
                 int nsteps,
                 int stepoffsets,
                 int& n1_out,                         // OUT
                 int& n2_out,                         // OUT
                 int& ne_out)                         // OUT
{
    const int nprocs = static_cast<int>(elempartpts.size());
    if (nprocs <= 0) throw std::runtime_error("elempartpts is empty.");
    if (static_cast<int>(elempart.size()) != nprocs) throw std::runtime_error("elempart size mismatch.");
    if (nsteps <= 0) throw std::runtime_error("nsteps must be positive.");
    if (stepoffsets < 0) throw std::runtime_error("stepoffsets must be nonnegative.");

    // ---- Pass 1: read headers and compute global sizes ----
    std::vector<std::array<std::int64_t, 4>> k(static_cast<std::size_t>(nprocs)); // {n1,n2,n3,timesteps}
    std::int64_t neTotal = 0;

    for (int r = 0; r < nprocs; ++r) {
        const std::string fname = base + "_np" + std::to_string(r) + ".bin";
        std::ifstream in(fname, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open file: " + fname);

        double hdr[3] = {0,0,0};
        readDoubles(in, hdr, 3, fname);

        const std::int64_t n1 = static_cast<std::int64_t>(hdr[0]);
        const std::int64_t n2 = static_cast<std::int64_t>(hdr[1]);
        const std::int64_t n3 = static_cast<std::int64_t>(hdr[2]);
        if (n1 <= 0 || n2 <= 0 || n3 <= 0) throw std::runtime_error("Invalid header in: " + fname);

        const std::int64_t N = n1 * n2 * n3;

        const std::int64_t bytes = fileSizeBytes(fname);
        if (bytes % 8 != 0) throw std::runtime_error("File size not multiple of 8 in: " + fname);

        const std::int64_t L = bytes / 8; // # of doubles
        if (L < 3) throw std::runtime_error("File too small in: " + fname);
        if ((L - 3) % N != 0) throw std::runtime_error("Payload not divisible by N in: " + fname);

        const std::int64_t timesteps = (L - 3) / N;
        k[static_cast<std::size_t>(r)] = {n1, n2, n3, timesteps};
        neTotal += n3;

        // std::cout << "N=" << N << " L=" << L<<"\n";
        // std::cout << "n1=" << n1 << " n2=" << n2 << " n3=" << n3 << " timesteps=" << timesteps << "\n";

        // Optional consistency check with elempartpts:
        const int Nj_pts = elempartpts[r][0] + elempartpts[r][1];
        if (Nj_pts != static_cast<int>(n3)) {
            // Not fatal if your elempartpts isn't meant to match n3, but in your MATLAB it usually should.
            // If you want strict, replace with throw.
            // throw std::runtime_error("elempartpts sum != n3 in rank " + std::to_string(r));
        }
        if (static_cast<int>(elempart[r].size()) < static_cast<int>(n3)) {
            throw std::runtime_error("elempart[" + std::to_string(r) + "] shorter than n3.");
        }
        if (timesteps < static_cast<std::int64_t>(stepoffsets + nsteps)) {
            throw std::runtime_error("Requested steps exceed available timesteps in: " + fname);
        }
    }

    const int n1 = static_cast<int>(k[0][0]);
    const int n2 = static_cast<int>(k[0][1]);

    for (int r = 0; r < nprocs; ++r) {
        if (k[static_cast<std::size_t>(r)][0] != n1 || k[static_cast<std::size_t>(r)][1] != n2) {
            throw std::runtime_error("n1/n2 mismatch across rank files.");
        }
    }

    const int ne = static_cast<int>(neTotal);

    n1_out = n1;
    n2_out = n2;
    ne_out = ne;

    // Allocate sol(n1,n2,ne,nsteps) flattenedx, column-major
    const std::size_t slice = static_cast<std::size_t>(n1) * static_cast<std::size_t>(n2);
    const std::size_t total = slice * static_cast<std::size_t>(ne) * static_cast<std::size_t>(nsteps);
    sol3dGlobal.assign(total, 0.0);

    // ---- Pass 2: read each rank and scatter directly into global element indices ----
    for (int r = 0; r < nprocs; ++r) {
        const std::string fname = base + "_np" + std::to_string(r) + ".bin";
        std::ifstream in(fname, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open file: " + fname);

        double hdr[3] = {0,0,0};
        readDoubles(in, hdr, 3, fname);

        const int n3 = static_cast<int>(hdr[2]);
        const std::size_t N = slice * static_cast<std::size_t>(n3);

        // Skip stepoffsets steps (each step = N doubles)
        if (stepoffsets > 0) {
            const std::int64_t skipBytes = static_cast<std::int64_t>(stepoffsets) *
                                           static_cast<std::int64_t>(N) * 8;
            in.seekg(skipBytes, std::ios::cur);
            if (!in) throw std::runtime_error("seekg failed in: " + fname);
        }

        std::vector<double> tm(N);

        for (int s = 0; s < nsteps; ++s) {
            readDoubles(in, tm.data(), N, fname);

            // tm is reshape(tm,[n1,n2,n3]) column-major:
            // slice e is contiguous block tm[e*slice : (e+1)*slice)
            for (int e = 0; e < n3; ++e) {
                const int globalElem = elempart[r][e] - 1; // MATLAB 1-based -> 0-based
                if (globalElem < 0 || globalElem >= ne) {
                    throw std::runtime_error("elempart out of range in rank " + std::to_string(r));
                }

                const double* src = tm.data() + static_cast<std::size_t>(e) * slice;
                double* dst = sol3dGlobal.data() + sliceOffset4(n1, n2, ne, globalElem, s);

                std::copy(src, src + slice, dst);
            }
        }
    }
}