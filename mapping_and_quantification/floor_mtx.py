"""
floor_mtx.py

Read a Matrix Market (.mtx) file, apply floor to all stored entries, and write a new .mtx.

Written by ChatGPT.

Usage:
  python3 floor_mtx.py input.mtx output.mtx
  python3 floor_mtx.py input.mtx output.mtx --dtype int64

Notes:
 - This script always uses floor (np.floor) and does NOT clamp negatives.
 - It will abort with an error if any negative values are found (you said we can assume none).
 - Default integer dtype is int32.
"""
import argparse
import sys
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import coo_matrix

def main():
    p = argparse.ArgumentParser(description="Apply floor to counts in a Matrix Market .mtx file")
    p.add_argument("input_mtx", help="Input .mtx file")
    p.add_argument("output_mtx", help="Output .mtx file (floored)")
    p.add_argument("--dtype", choices=["int32","int64"], default="int32",
                   help="Integer dtype for output (default: int32)")
    args = p.parse_args()

    mat = mmread(args.input_mtx)
    if not hasattr(mat, "tocoo"):
        mat = coo_matrix(mat)
    else:
        mat = mat.tocoo()

    data = mat.data.astype(np.float64)
    print(f"Input shape: {mat.shape}, nnz: {data.size}", file=sys.stderr)
    print(f"Input data min/max: {data.min():.6g} / {data.max():.6g}", file=sys.stderr)

    if np.any(data < 0):
        print("Error: input contains negative values; this script assumes no negatives. Aborting.", file=sys.stderr)
        sys.exit(2)

    floored = np.floor(data).astype(getattr(np, args.dtype))

    newmat = coo_matrix((floored, (mat.row, mat.col)), shape=mat.shape)
    mmwrite(args.output_mtx, newmat)
    print(f"Wrote floored matrix to {args.output_mtx}", file=sys.stderr)
    print(f"Output data min/max: {int(floored.min())} / {int(floored.max())}", file=sys.stderr)

if __name__ == "__main__":
    main()
