#!/usr/bin/env python3
"""Generate grid.bin for the tutorial: an n x n Cartesian quad mesh
on the unit square, in the legacy ``writebin``-compatible binary
format the Exasim runtime expects.

Format (all values stored as float64, column-major order):
    [nd, np, nve, ne, p(:, :), t(:, :)]
where ``p`` has shape (nd, np) — node coordinates — and ``t`` has
shape (nve, ne) — 1-based element-vertex incidence.

Imports ``SquareMesh`` from ``frontends/Python/Mesh/squaremesh.py``
(pure NumPy; no gmsh dependency for the rectangular square).

Usage:
    python3 squaregrid.py <n> <output-path>
"""

import os
import sys

import numpy as np


def main(argv):
    if len(argv) != 3:
        print(f"Usage: {argv[0]} <n> <output-path>", file=sys.stderr)
        return 2

    n = int(argv[1])
    out_path = argv[2]

    this_dir = os.path.dirname(os.path.abspath(__file__))
    exasim_dir = os.path.abspath(os.path.join(this_dir, "..", ".."))
    mesh_dir = os.path.join(exasim_dir, "frontends", "Python", "Mesh")
    sys.path.insert(0, mesh_dir)
    from squaremesh import SquareMesh  # noqa: E402

    p, t = SquareMesh(n, n, 1)  # quad elements; returns p (nd, np), t (nve, ne) zero-based
    t = t + 1                   # 1-based, matching legacy MATLAB writebin convention

    header = np.array([p.shape[0], p.shape[1], t.shape[0], t.shape[1]],
                      dtype=np.float64)
    blob = np.concatenate([
        header,
        p.flatten(order="F").astype(np.float64),
        t.flatten(order="F").astype(np.float64),
    ])
    blob.tofile(out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
