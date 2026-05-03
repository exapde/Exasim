# 03 — Generated model with the embedded solver

The generated model is the same as the prebuilt and CLI sections,
but here the user constructs an `ExasimSolver<GeneratedModel>`
explicitly inside `main.cpp` instead of letting `exasim::run`
hide it. The solver becomes one component inside a larger C++
program: optimization loops, reduced-order modeling, coupled
multi-physics, and so on.

This section uses the simplest embedded variant: `load_pdeapp`
parses `pdeapp.txt` from disk and stages mesh and configuration on
the solver. Sections 05 and 06 show embedded variants that
construct the mesh from C++ arrays without touching `pdeapp.txt`.

## Files

- `pdemodel.txt` — the PDE expressed in the SymEngine DSL.
- `pdeapp.txt` — runtime configuration parsed by `load_pdeapp`.
- `main.cpp` — embedded driver. Constructs `ExasimSolver`, calls
  `load_pdeapp(argv[1])`, calls `solve()`, finalizes.
- `CMakeLists.txt` — out-of-tree build for
  `tutorial_03_generated_embedded`.
- `grid.bin` — binary mesh for the unit square.

## Build and run

`text2code` must run before `cmake --build` so that `my_model.hpp`
and `libpdemodelserial.{so,dylib}` exist in this directory at
link time.

### CPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/03-generated-embedded \
    $EXASIM/tutorial/03-generated-embedded/pdeapp.txt
cmake --build build --target tutorial_03_generated_embedded
cd $EXASIM/tutorial/03-generated-embedded
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build/tutorial_03_generated_embedded ./pdeapp.txt
```

### GPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/03-generated-embedded \
    $EXASIM/tutorial/03-generated-embedded/pdeapp.txt
cmake --build build_gpu --target tutorial_03_generated_embedded
cd $EXASIM/tutorial/03-generated-embedded
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build_gpu/tutorial_03_generated_embedded ./pdeapp.txt
```

### MPI

Edit `pdeapp.txt` so `mpiprocs = N` matches the rank count.

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/03-generated-embedded \
    $EXASIM/tutorial/03-generated-embedded/pdeapp.txt
cmake --build build_mpi --target tutorial_03_generated_embedded
cd $EXASIM/tutorial/03-generated-embedded
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
mpirun -np 2 $EXASIM/build_mpi/tutorial_03_generated_embedded ./pdeapp.txt
```

### MPI+GPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/03-generated-embedded \
    $EXASIM/tutorial/03-generated-embedded/pdeapp.txt
cmake --build build_mpi_gpu --target tutorial_03_generated_embedded
cd $EXASIM/tutorial/03-generated-embedded
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
mpirun -np 2 $EXASIM/build_mpi_gpu/tutorial_03_generated_embedded ./pdeapp.txt
```

## File contents

### `pdemodel.txt`

```
scalars t

vectors x(2), uq(3), v(0), w(0), uhat(1), uext(1), n(2), tau(1), mu(1), eta(0)

jacobian uq, w, uhat 

hessian 

batch x, uq, v, w, uhat, n, uext

outputs Flux, Source, Tdfunc, Ubou, Fbou, Fint, Fext, FbouHdg, Initu, VisScalars, VisVectors, QoIvolume, QoIboundary

function Flux(x, uq, v, w, eta, mu, t)  
  output_size(f) = 2;
  kappa = mu[0];
  f[0] = kappa*uq[1];
  f[1] = kappa*uq[2];
end

function Source(x, uq, v, w, eta, mu, t)  
  output_size(s) = 1;
  x1 = x[0];
  x2 = x[1];
  s[0] = mul(mul(SymEngine::integer(2),mul(pi, pi)) , mul(sin(pi*x1),sin(pi*x2)));
end

function Tdfunc(x, uq, v, w, eta, mu, t)  
  output_size(m) = 1;
  ones(m);
end

function Fbou(x, uq, v, w, uhat, n, tau, eta, mu, t)
  output_size(fb) = 1;
  f = Flux(x, uq, v, w, eta, mu, t);  
  fb[0] = f[0]*n[0] + f[1]*n[1] + tau[0]*(uq[0]-uhat[0]);
end

function Ubou(x, uq, v, w, uhat, n, tau, eta, mu, t)
  output_size(ub) = 1;
  ub[0] = 0.0;
end

function FbouHdg(x, uq, v, w, uhat, n, tau, eta, mu, t)
  output_size(fb) = 1;
  fb[0] = tau[0]*(0.0 - uhat[0]);  
end

function Initu(x, eta, mu)
  output_size(ui) = 1;
  ui[0] = 0.0;
end

function VisScalars(x, uq, v, w, eta, mu, t)  
  output_size(s) = 2;
  s[0] = uq[0];
  s[1] = uq[1] + uq[2];  
end

function VisVectors(x, uq, v, w, eta, mu, t)  
  output_size(s) = 2;
  s[0] = uq[1];
  s[1] = uq[2];  
end

function QoIvolume(x, uq, v, w, eta, mu, t)  
  output_size(s) = 2;
  x1 = x[0];
  x2 = x[1];  
  t1 = Expression(pi);  
  t2 = sin(t1*x1);
  t3 = sin(t1*x2);  
  uexact = mul(t2,t3); 
  s[0] = (uq[0] - uexact)*(uq[0] - uexact);
  s[1] = uq[0];
end
  
function QoIboundary(x, uq, v, w, uhat, n, tau, eta, mu, t)
  output_size(fb) = 1;
  f = Flux(x, uq, v, w, eta, mu, t);  
  fb[0] = f[0]*n[0] + f[1]*n[1] + tau[0]*(uq[0]-uhat[0]);
end

function Fint(x, uq, v, w, uhat, n, tau, eta, mu, t)
  output_size(fb) = 2;
  fb[0] = (0-uhat[0]);
  fb[1] = (1-uhat[0]);
end

function Fext(x, uq, v, w, uhat, n, uext, tau, eta, mu, t)
  output_size(fb) = 1;
  fb[0] = (uext[0]-uhat[0]);
end
```

`Flux` returns the 2D flux `μ ∇u`. `Source` returns the
manufactured forcing `2π² sin(πx) sin(πy)`. `FbouHdg` returns the
zero-Dirichlet HDG residual `−τ ûh`. `QoIvolume` returns
`(u − u_exact)²` so the runtime can integrate the squared L²
error.

### `pdeapp.txt`

```
model = "ModelD";
modelfile = "pdemodel.txt";
meshfile = "grid.bin";
discretization = "hdg";
platform = "cpu";

builtinmodelID = 1;
gendatain = 0;
mpiprocs = 1;
debugmode = 0;
runmode = 0;
modelnumber = 0;

porder = 3;
pgauss = 6;
torder = 1;
nstage = 1;

ncu = 1;
ncw = 0;
neb = 4096;
nfb = 8192;

saveSolBouFreq = 1;
ibs = 1;

NewtonIter = 20;
NewtonTol = 1e-08;
GMRESiter = 50;
GMRESrestart = 25;
GMREStol = 1e-08;
GMRESortho = 1;
ppdegree = 1;
RBdim = 0;
matvecorder = 1;
matvectol = 0.001;
precMatrixType = 0;
preconditioner = 1;

time = 0;
tau = [1];
dt = [1];
physicsparam = [1.0];
externalparam = [0, 0];

boundaryconditions = [1];
boundaryexpressions = ["abs(y)<1e-8 || abs(y-1)<1e-8 || abs(x-1)<1e-8 || abs(x)<1e-8"];
curvedboundaries = [0];
curvedboundaryexprs = [""];
periodicboundaries1 = [];
periodicexprs1 = [];
periodicboundaries2 = [];
periodicexprs2 = [];
interfaceconditions = [];
```

`load_pdeapp(argv[1], mpirank)` parses this file and stages it on
the solver. `model = "ModelD"` selects the steady-state HDG class.
`porder = 3` and `pgauss = 6` set discretization order and
quadrature order. `physicsparam = [1.0]` sets `mu[0] = 1`.
`boundaryconditions = [1]` and `boundaryexpressions` tag every
face on the four edges of the unit square with id 1, which the
model handles as zero-Dirichlet. For multi-rank runs `mpiprocs = N`
must match the rank count.

### `main.cpp`

```cpp
#include <exasim/run.hpp>
#include <exasim/solver_facade.hpp>
#include <exasim/model.hpp>

#include "my_model.hpp"

#include <cstdio>

int main(int argc, char** argv) {
    int mpiprocs = 1, mpirank = 0, shmrank = 0;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
    MPI_Comm shmcomm;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
#endif
#ifdef HAVE_CUDA
    int nd; cudaGetDeviceCount(&nd); cudaSetDevice(shmrank % nd);
#endif
#ifdef HAVE_HIP
    int nd; hipGetDeviceCount(&nd); hipSetDevice(shmrank % nd);
#endif
    (void)shmrank;

    if (argc < 2) {
        if (mpirank == 0) std::fprintf(stderr, "Usage: %s pdeapp.txt\n", argv[0]);
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    Kokkos::initialize();
    {
        exasim::ExasimSolver<GeneratedModel> solver;
        solver.pde().saveOutputs = 1;
        solver.load_pdeapp(argv[1], mpirank);
        solver.solve(mpiprocs, mpirank);
    }
    Kokkos::finalize();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
```

Compared to section 02's three-line `main.cpp`, this version takes
explicit responsibility for MPI and Kokkos setup so the user can
weave the solve into a larger program.

The `HAVE_MPI` block initializes MPI and stores the world
communicator on the global handles `EXASIM_COMM_WORLD` and
`EXASIM_COMM_LOCAL`. The `MPI_Comm_split_type(...
MPI_COMM_TYPE_SHARED ...)` call gives each rank a node-local
communicator, and `shmrank` is its position inside that node. The
`HAVE_CUDA` and `HAVE_HIP` blocks bind one device per rank using
that local rank index, so multi-GPU nodes spread their ranks
across devices.

`exasim::ExasimSolver<GeneratedModel>` is the embedded solver
template. `solver.pde().saveOutputs = 1` tells the solver to write
`dataout/outudg_np<rank>.bin` and friends. `load_pdeapp(argv[1],
mpirank)` parses `pdeapp.txt`, reads the mesh from
`pdeapp.meshfile`, and stages everything on the solver.
`solve(mpiprocs, mpirank)` runs the Newton/GMRES loop. `Kokkos`
and `MPI` are finalized before exit.

`my_model.hpp` is generated by `text2code` from `pdemodel.txt`;
this directory does not contain a checked-in copy.

### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.16)

set(_target tutorial_03_generated_embedded)

add_executable(${_target} main.cpp)
tutorial_configure_target(${_target})

target_link_directories(${_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
```

`tutorial_configure_target` is a helper defined in
`tutorial/CMakeLists.txt` that adds the right backend defines and
libraries for the active build variant: `_CUDA` and the CUDA
runtime/cuBLAS libraries on `build_gpu` and `build_mpi_gpu`,
`_HIP` on AMD GPUs, `_MPI` on MPI-enabled builds, and
`_TEXT2CODE` everywhere. The local `target_link_directories` and
`BUILD_RPATH` overrides point the linker at this section's own
`libpdemodel*.{so,dylib}`.

### `grid.bin`

The mesh file referenced by `meshfile = "grid.bin"` in
`pdeapp.txt`. Not checked in. Regenerate with:

```bash
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
```

The script imports `SquareMesh` from
`frontends/Python/Mesh/squaremesh.py` (pure NumPy, no `gmsh`)
and writes a 16x16 Cartesian quad mesh on the unit square in the
legacy `writebin` format `[size(p), size(t), p(:), t(:)]`. See
[section gridbin](../00-overview/README.md#gridbin) of the overview
for the format details.
