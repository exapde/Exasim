# AGENTS.md
Global instructions for AI coding agents assisting with Exasim and related
scientific computing projects.

Domain
------
This project involves:

• high-order finite element methods
• hybridizable discontinuous Galerkin (HDG) discretizations
• numerical PDE solvers
• GPU acceleration (CUDA / HIP / Kokkos)
• C++ HPC software

Agents must prioritize scientific correctness and numerical stability.

---------------------------------------------------------------------

PRIORITIES (highest first)
--------------------------

1. Correctness

Never trade correctness for performance.

Numerical methods must preserve:
• conservation
• consistency
• stability
• dimensional correctness

If uncertain about a change that may affect correctness,
the agent must explicitly state the uncertainty.

---------------------------------------------------------------------

2. Explicit Assumptions

All reasoning must clearly state:

• mathematical assumptions
• data layout assumptions
• memory ownership assumptions
• solver state assumptions

---------------------------------------------------------------------

3. GPU-aware implementation

GPU code must consider:

• memory bandwidth
• memory coalescing
• register pressure
• warp divergence
• kernel launch overhead
• host-device transfers

Prefer designs that maximize arithmetic intensity and minimize
global memory traffic.

---------------------------------------------------------------------

4. Minimal invasive changes

Avoid refactoring unrelated code.

Changes must:

• preserve existing interfaces whenever possible
• avoid altering solver structure
• minimize code footprint
• maintain backward compatibility

---------------------------------------------------------------------

5. Clear verification steps

Every proposed implementation must include:

• edge case analysis
• numerical validation plan
• regression tests
• performance sanity checks

---------------------------------------------------------------------

NUMERICAL METHOD CONTEXT
------------------------

The codebase implements HDG methods.

Key characteristics:

• element-local solves
• static condensation
• global trace system
• high-order polynomial basis

When modifying algorithms, maintain these structural properties.

The HDG solve typically involves:

1. element-local residual and Jacobian assembly
2. static condensation eliminating interior DOFs
3. global trace system solve
4. recovery of element solutions

Agents must not break this structure.

---------------------------------------------------------------------

DATA LAYOUT CONVENTIONS
-----------------------

Follow Exasim conventions:

• arrays are typically column-major
• indexing is zero-based
• flattened arrays are common
• batched operations are preferred

Avoid introducing new data layouts unless necessary.

---------------------------------------------------------------------

GPU IMPLEMENTATION GUIDELINES
------------------------------

When implementing GPU kernels:

Prefer:

• batched dense linear algebra
• BLAS operations (cuBLAS / hipBLAS)
• data-parallel loops
• kernel fusion when beneficial

Avoid:

• excessive kernel launches
• uncoalesced memory access
• host-device synchronization

Consider using:

• batched GEMM
• batched triangular solves
• shared memory for small dense blocks

---------------------------------------------------------------------

SOLVER MODIFICATION POLICY
--------------------------

Before modifying solver code, the agent must:

1. summarize the proposed design
2. list affected files
3. explain how the change interacts with the HDG solver
4. confirm assumptions about solver state

Only then propose code edits.

---------------------------------------------------------------------

CODE CHANGE FORMAT
------------------

All code proposals must include:

1. Design summary
2. Assumptions
3. Affected files
4. Implementation approach
5. Code diff or patch
6. Verification plan

---------------------------------------------------------------------

PERFORMANCE ANALYSIS
--------------------

When analyzing performance:

Check:

• arithmetic intensity
• memory bandwidth limits
• roofline position
• kernel launch counts
• GPU occupancy

Use profiler outputs (Nsight, rocprof, etc.) when available.

---------------------------------------------------------------------

DEBUGGING GUIDELINES
--------------------

When debugging numerical issues:

Check:

• boundary conditions
• flux consistency
• Jacobian correctness
• conservation violations
• NaN or negative density/energy states

Always identify likely root causes before proposing fixes.

---------------------------------------------------------------------

PAPER WRITING SUPPORT
---------------------

When assisting with manuscripts:

• maintain mathematical rigor
• clearly distinguish assumptions
• avoid unsupported claims
• use consistent notation

Prefer structured explanations.

---------------------------------------------------------------------

RESPONSE FORMAT
---------------

All responses must follow this structure:

1. Summary
2. Assumptions
3. Design reasoning
4. Proposed implementation
5. Risks or numerical concerns
6. Verification plan

---------------------------------------------------------------------

IF INFORMATION IS MISSING
-------------------------

If the agent lacks sufficient information to proceed safely:

• state the missing information
• propose diagnostic steps
• avoid speculative changes

---------------------------------------------------------------------

END OF AGENTS.md