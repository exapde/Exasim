# 02 — Generated model with the CLI wrapper

The same `pdemodel.txt` and `pdeapp.txt` as the prebuilt section,
but the user writes their own three-line `main.cpp` that links
directly against Exasim's header-only library. The result is a
standalone binary that the user owns: it can be renamed, installed,
or shipped, and its build system is local rather than tied to
Exasim's main install.

`text2code` reads `pdemodel.txt` and emits `my_model.hpp`
(containing the `GeneratedModel` struct) plus
`libpdemodelserial.{so,dylib}`, which `main.cpp` links against.
`main.cpp` then calls `exasim::run<GeneratedModel>(argc, argv)`,
which parses `pdeapp.txt`, runs preprocessing, solves, and writes
outputs to `dataout/`.

## Files

- `pdemodel.txt` — the PDE expressed in the SymEngine DSL.
- `pdeapp.txt` — runtime configuration consumed by `exasim::run`.
- `main.cpp` — three-line entry point.
- `CMakeLists.txt` — out-of-tree build for `tutorial_02_generated_cli`.
- `grid.bin` — binary mesh for the unit square.

## Build and run

`text2code` must run before `cmake --build` so that `my_model.hpp`
and `libpdemodelserial.{so,dylib}` exist in this directory at
link time.

### CPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/02-generated-cli \
    $EXASIM/tutorial/02-generated-cli/pdeapp.txt
cmake --build build --target tutorial_02_generated_cli
cd $EXASIM/tutorial/02-generated-cli
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build/tutorial_02_generated_cli ./pdeapp.txt
```

### GPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/02-generated-cli \
    $EXASIM/tutorial/02-generated-cli/pdeapp.txt
cmake --build build_gpu --target tutorial_02_generated_cli
cd $EXASIM/tutorial/02-generated-cli
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build_gpu/tutorial_02_generated_cli ./pdeapp.txt
```

### MPI

Edit `pdeapp.txt` so `mpiprocs = N` matches the rank count.

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/02-generated-cli \
    $EXASIM/tutorial/02-generated-cli/pdeapp.txt
cmake --build build_mpi --target tutorial_02_generated_cli
cd $EXASIM/tutorial/02-generated-cli
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
mpirun -np 2 $EXASIM/build_mpi/tutorial_02_generated_cli ./pdeapp.txt
```

### MPI+GPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/02-generated-cli \
    $EXASIM/tutorial/02-generated-cli/pdeapp.txt
cmake --build build_mpi_gpu --target tutorial_02_generated_cli
cd $EXASIM/tutorial/02-generated-cli
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
mpirun -np 2 $EXASIM/build_mpi_gpu/tutorial_02_generated_cli ./pdeapp.txt
```

## File contents

### `pdemodel.txt`

Identical in content to section 01's `pdemodel.txt`. The full
listing is reproduced here so this section is self-contained.

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

`exasim::run<GeneratedModel>(argc, argv)` reads this file from
`argv[1]`. The fields control the discretization (`porder`,
`pgauss`, `discretization = "hdg"`), the solver loop (Newton +
GMRES tolerances), the boundary tagging (every edge of the unit
square is tagged with boundary id 1, which the model handles as
zero-Dirichlet), and the physics parameters (`physicsparam = [1.0]`
sets `mu[0] = 1`). For multi-rank runs `mpiprocs = N` must match
the rank count.

### `main.cpp`

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    return exasim::run<GeneratedModel>(argc, argv);
}
```

`exasim::run<M>(argc, argv)` does the MPI / Kokkos init, parses
`argv[1]` as a `pdeapp.txt`, builds the `CSolution<M>` object,
runs the Newton/GMRES solve (or DIRK time loop for transient
problems), saves outputs, and finalizes.

`my_model.hpp` is generated by `text2code` from `pdemodel.txt`;
this directory does not contain a checked-in copy. The build
system requires `text2code` to run first.

### `CMakeLists.txt`

This is the standalone CMakeLists a real out-of-tree consumer of
Exasim would write. The in-tree tutorial build does not consume it
— it registers the same target via `tutorial/CMakeLists.txt`.

```cmake
cmake_minimum_required(VERSION 3.16)
project(tutorial_02_generated_cli CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(Exasim REQUIRED)
find_package(Kokkos REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)
target_compile_definitions(${PROJECT_NAME} PRIVATE _TEXT2CODE)
target_link_libraries(${PROJECT_NAME} PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_link_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${PROJECT_NAME} PRIVATE pdemodelserial)
set_target_properties(${PROJECT_NAME} PROPERTIES
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
```

`find_package(Exasim REQUIRED)` finds the installed Exasim from
`cmake --install $EXASIM/build --prefix /opt/exasim`; it provides
the `Exasim::headers` INTERFACE target. The user picks their own
Kokkos build, BLAS, and LAPACK — Exasim does not propagate
`find_dependency` for those.

The `target_link_directories` line points the linker at the local
`libpdemodelserial.{so,dylib}` that `text2code` emits when run with
`--out-dir .` against this directory's `pdeapp.txt`. The
`BUILD_RPATH` makes the runtime loader find that library at
launch.

To build standalone:

```bash
cmake --install $EXASIM/build --prefix /opt/exasim
$EXASIM/build/text2code --out-dir . ./pdeapp.txt
cmake -S . -B build \
      -DCMAKE_PREFIX_PATH=/opt/exasim \
      -DKokkos_DIR=/opt/exasim/external/kokkos
cmake --build build
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
./build/tutorial_02_generated_cli ./pdeapp.txt
```

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
