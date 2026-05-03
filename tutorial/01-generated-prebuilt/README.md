# 01 — Generated model with the prebuilt solver

The most automated entry point. The user writes two text files
(`pdemodel.txt` describing the math, `pdeapp.txt` configuring the
runtime), and points one of Exasim's prebuilt binaries at them.
There is no C++ to write, no `CMakeLists.txt` to author, and no
build step on the user's side.

The prebuilt binaries dispatch every PDE kernel through the
`libpdemodel.{so,dylib}` ABI that `text2code` emits from
`pdemodel.txt`. The runtime under the hood is `backend/Main/main.cpp`,
which calls `exasim::run<exasim::detail::AbiAdapter>(argc, argv)`;
the `AbiAdapter` marker type instructs the FEM templates to load
the model symbols from the shared library at startup.

The prebuilt binaries are produced by the main install when CMake
is configured with `-DWITH_TEXT2CODE=ON`:

| Backend | Binary |
|---|---|
| CPU | `cput2cEXASIM` |
| GPU | `gput2cEXASIM` |
| MPI | `cpumpit2cEXASIM` |
| MPI+GPU | `gpumpit2cEXASIM` |

## Files

- `pdemodel.txt` — the PDE expressed in the SymEngine DSL.
- `pdeapp.txt` — runtime configuration consumed by the binary.
- `grid.bin` — binary mesh for the unit square. Not checked in;
  regenerate with the tutorial generator described in
  [Files: `grid.bin`](#gridbin) below.

## Build and run

The binary's `RPATH` is hard-coded to `backend/Model/`, so this
section runs `text2code` without `--out-dir` to populate that
location with `libpdemodelserial.{so,dylib}` and friends.

### CPU

```bash
cd $EXASIM
cmake -S install -B build -DWITH_TEXT2CODE=ON
cmake --build build --target cput2cEXASIM
cd $EXASIM/tutorial/01-generated-prebuilt
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build/text2code ./pdeapp.txt
$EXASIM/build/cput2cEXASIM ./pdeapp.txt
```

### GPU

```bash
cd $EXASIM
cmake -S install -B build_gpu -DWITH_TEXT2CODE=ON -DEXASIM_NOMPI=ON -DEXASIM_CUDA=ON
cmake --build build_gpu --target gput2cEXASIM
cd $EXASIM/tutorial/01-generated-prebuilt
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build/text2code ./pdeapp.txt
$EXASIM/build_gpu/gput2cEXASIM ./pdeapp.txt
```

### MPI

Edit `pdeapp.txt` so `mpiprocs = N` matches the rank count.

```bash
cd $EXASIM
cmake -S install -B build_mpi -DWITH_TEXT2CODE=ON -DEXASIM_MPI=ON
cmake --build build_mpi --target cpumpit2cEXASIM
cd $EXASIM/tutorial/01-generated-prebuilt
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build/text2code ./pdeapp.txt
mpirun -np 2 $EXASIM/build_mpi/cpumpit2cEXASIM ./pdeapp.txt
```

### MPI+GPU

```bash
cd $EXASIM
cmake -S install -B build_mpi_gpu -DWITH_TEXT2CODE=ON -DEXASIM_MPI=ON -DEXASIM_CUDA=ON
cmake --build build_mpi_gpu --target gpumpit2cEXASIM
cd $EXASIM/tutorial/01-generated-prebuilt
mkdir -p datain dataout
python3 $EXASIM/tutorial/tools/squaregrid.py 16 ./grid.bin
$EXASIM/build/text2code ./pdeapp.txt
mpirun -np 2 $EXASIM/build_mpi_gpu/gpumpit2cEXASIM ./pdeapp.txt
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

The `vectors` and `scalars` lines declare the symbolic identifiers
text2code uses: `x` is the spatial coordinate, `uq` is the mixed
state `[u, ∂u/∂x, ∂u/∂y]`, `mu[0]` is the diffusion coefficient,
`uhat` is the HDG trace, `n` is the outward normal, `tau` is the
HDG stabilization parameter. The `outputs` line lists which kernel
functions text2code should emit. The `Flux` function returns the
2D flux `μ ∇u`, the `Source` function returns the manufactured
forcing `2π² sin(πx) sin(πy)`, and `FbouHdg` enforces the
zero-Dirichlet HDG residual `−τ ûh`. The `QoIvolume` function
returns `(u − u_exact)²` so the runtime can integrate the squared
L² error.

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

`model = "ModelD"` selects the steady-state HDG model class.
`modelfile = "pdemodel.txt"` points at the symbolic description.
`meshfile = "grid.bin"` is the binary mesh; see
[section gridbin](../00-overview/README.md#gridbin) of the
overview for where this file comes from. `porder = 3` and
`pgauss = 6` set the polynomial order and quadrature order.
`physicsparam = [1.0]` sets `mu[0] = 1`. `boundaryconditions = [1]`
and the corresponding `boundaryexpressions` tag every face on the
four edges of the unit square with boundary tag 1, which is what
`Fbou` and `FbouHdg` treat as a Dirichlet boundary. The `mpiprocs`
field must match the number of ranks at run time; for the
multi-rank variants this section's `pdeapp.txt` should be edited
to set `mpiprocs = N` before launching `mpirun -np N`.

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
