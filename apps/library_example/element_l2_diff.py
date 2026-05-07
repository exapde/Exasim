#!/usr/bin/env python3
"""
HOT.7.18 — element-aware L2 diff for DG outputs.

Aggregates outudg_np<r>.bin + outelemid_np<r>.bin across ranks,
sorts by global element ID, then compares against a (serial-recorded)
baseline element-by-element with a relative L2 metric.

Usage:
    element_l2_diff.py --baseline-dir <dir> --current-dir <dir> \
        --np <N> [--rtol 1e-3] [--stem outudg]

Exits 0 on pass, non-zero on fail. Prints a one-line summary
suitable for the validate_codegen.sh harness.

Why this exists: bin-md5 baselines are sensitive to FP reduction
order (Apple Accelerate vs Intel MKL drift through GMRES+Newton
gives ~1e-3 last-bit divergence) and to MPI partition assignment
(ParMETIS shuffles local order). Element-level L2 sidesteps both:
ordering doesn't matter once we sort by global ID, and last-bit
drift sums sub-linearly into a relative metric far below the 1e-3
threshold.
"""

import argparse
import os
import struct
import sys


def read_doubles(path):
    with open(path, "rb") as f:
        data = f.read()
    if len(data) % 8 != 0:
        sys.exit(f"{path}: size {len(data)} not divisible by 8")
    return struct.unpack(f"{len(data)//8}d", data)


def read_ids(path):
    with open(path, "rb") as f:
        data = f.read()
    if len(data) % 8 != 0:
        sys.exit(f"{path}: size {len(data)} not divisible by 8")
    return struct.unpack(f"{len(data)//8}q", data)  # int64


def aggregate(udg_path_fn, id_path_fn, np_count):
    """Return dict: global_id -> tuple of doubles (one element block).

    The udg files start with a 3-double header [a0, a1, a2] written by
    `open_and_write` in solution.h (rank, npe-style metadata). We strip
    that header and split the remaining payload by ne_local equal blocks.
    """
    blocks = {}
    for r in range(np_count):
        upath = udg_path_fn(r)
        ipath = id_path_fn(r)
        if not os.path.exists(upath):
            continue
        if not os.path.exists(ipath):
            sys.exit(f"missing elemid sidecar: {ipath}")
        u_full = read_doubles(upath)
        ids = read_ids(ipath)
        if not ids:
            continue
        if len(u_full) < 3:
            sys.exit(f"{upath}: file too small ({len(u_full)} doubles, need >= 3 header)")
        u = u_full[3:]
        ne_local = len(ids)
        if len(u) % ne_local != 0:
            sys.exit(f"{upath}: {len(u)} payload doubles not divisible by ne_local={ne_local}")
        block_size = len(u) // ne_local
        for i, gid in enumerate(ids):
            if gid in blocks:
                sys.exit(f"duplicate global id {gid} (rank {r}, local {i})")
            blocks[gid] = u[i*block_size:(i+1)*block_size]
    return blocks


def element_l2(a_blocks, b_blocks):
    """Relative L2: sqrt(sum_e ||a-b||^2 / sum_e ||a||^2). Common keys only.
    Returns (rel_l2, num_compared, num_only_a, num_only_b)."""
    a_keys = set(a_blocks)
    b_keys = set(b_blocks)
    common = a_keys & b_keys
    only_a = a_keys - b_keys
    only_b = b_keys - a_keys

    err_sq = 0.0
    ref_sq = 0.0
    for k in common:
        a = a_blocks[k]
        b = b_blocks[k]
        if len(a) != len(b):
            sys.exit(f"element {k}: shape {len(a)} vs {len(b)}")
        for x, y in zip(a, b):
            err_sq += (x - y) ** 2
            ref_sq += x * x
    if ref_sq <= 0:
        return (0.0 if err_sq == 0 else float("inf"),
                len(common), len(only_a), len(only_b))
    rel = (err_sq / ref_sq) ** 0.5
    return rel, len(common), len(only_a), len(only_b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-dir", required=True)
    ap.add_argument("--current-dir", required=True)
    ap.add_argument("--np", type=int, default=1, help="rank count of CURRENT run")
    ap.add_argument("--baseline-np", type=int, default=1)
    ap.add_argument("--rtol", type=float, default=1e-3)
    ap.add_argument("--stem", default="outudg")
    args = ap.parse_args()

    base = aggregate(
        lambda r: os.path.join(args.baseline_dir, f"{args.stem}_np{r}.bin"),
        lambda r: os.path.join(args.baseline_dir, f"outelemid_np{r}.bin"),
        args.baseline_np,
    )
    cur = aggregate(
        lambda r: os.path.join(args.current_dir, f"{args.stem}_np{r}.bin"),
        lambda r: os.path.join(args.current_dir, f"outelemid_np{r}.bin"),
        args.np,
    )

    if not base:
        print(f"baseline {args.stem} aggregate empty (no outelemid sidecar?)")
        sys.exit(2)
    if not cur:
        print(f"current {args.stem} aggregate empty")
        sys.exit(2)

    rel, n, na, nb = element_l2(base, cur)
    if na or nb:
        print(f"{args.stem}: element-id mismatch — baseline-only={na}, current-only={nb}")
        sys.exit(1)
    if rel > args.rtol:
        print(f"{args.stem}: relative element-L2 = {rel:.3e} (>{args.rtol:.0e}) over {n} elements")
        sys.exit(1)
    print(f"{args.stem}: relative element-L2 = {rel:.3e} over {n} elements")
    sys.exit(0)


if __name__ == "__main__":
    main()
