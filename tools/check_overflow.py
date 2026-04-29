#!/usr/bin/env python3
"""
Overflow checker for COLA Halo OMP based on configuration parameters.

Checks for potential integer overflow risks in:
1. Particle allocation (mem.c)
2. PM grid indexing (pm.c)
3. LPT grid indexing (lpt.c)

All checks consider the effect of MPI parallelization on local work distribution.
"""

import argparse
import sys


# Constants
INT_MAX = 2**31 - 1  # Maximum value for 32-bit signed integer
INT_MAX_SQRT = int(INT_MAX**0.5)  # ~46340


def check_particle_allocation(nc, mpi_proc, np_alloc_factor=1.25):
    """
    Check overflow risk in particle allocation (mem.c).

    The allocation formula is:
        np_alloc = np_alloc_factor * nc^2 * (local_nx + 1)

    where local_nx ≈ nc / mpi_proc (MPI decomposition in x-direction).

    Returns:
        dict with overflow status and details
    """
    import math
    local_nx = math.ceil(nc / mpi_proc)

    # np_alloc = np_alloc_factor * nc^2 * (local_nx + 1)
    np_alloc = int(np_alloc_factor * nc * nc * (local_nx + 1))

    # Check if np_alloc would overflow 32-bit int when cast
    overflow = np_alloc > INT_MAX

    # Calculate safe threshold for nc given mpi_proc and np_alloc_factor
    # np_alloc_factor * nc^2 * (nc/mpi_proc + 1) <= INT_MAX
    # Approximate: np_alloc_factor * nc^3 / mpi_proc <= INT_MAX
    # nc <= (INT_MAX * mpi_proc / np_alloc_factor)^(1/3)
    safe_nc_approx = int((INT_MAX * mpi_proc / np_alloc_factor) ** (1/3))

    return {
        "module": "Particle Allocation (mem.c)",
        "np_alloc": np_alloc,
        "local_nx": local_nx,
        "overflow": overflow,
        "safe_nc_approx": safe_nc_approx,
        "description": f"np_alloc = nc² × (local_nx+1) ≈ {np_alloc:,}",
    }


def check_pm_indexing(nc, pm_nc_factor):
    """
    Check overflow risk in PM grid indexing (pm.c).

    The PM grid uses Ngrid = nc * pm_nc_factor.
    Index calculation: (i * Ngrid + j) * (Ngrid/2 + 1) + k

    The loop variables i, j, k are 32-bit int, but the final index
    uses NgridL (size_t, 64-bit). However, intermediate calculation
    i * Ngrid + j can overflow if Ngrid >= 46341.

    Returns:
        dict with overflow status and details
    """
    ngrid = nc * pm_nc_factor

    # Check if Ngrid >= sqrt(INT_MAX) ≈ 46340
    # When j loops to Ngrid-1, i*Ngrid + j can overflow if Ngrid >= 46341
    overflow = ngrid > INT_MAX_SQRT

    # Safe threshold for nc given pm_nc_factor
    safe_nc = INT_MAX_SQRT // pm_nc_factor

    return {
        "module": "PM Grid Indexing (pm.c)",
        "ngrid": ngrid,
        "overflow": overflow,
        "safe_nc": safe_nc,
        "description": f"Ngrid = nc × pm_nc_factor = {ngrid:,}",
    }


def check_lpt_indexing(nc):
    """
    Check overflow risk in LPT grid indexing (lpt.c).

    The LPT grid uses Nmesh = nc.
    Index calculation: (i * Nmesh + j) * (Nmesh/2 + 1) + k

    Loop variables i, j, k are 32-bit int. The intermediate calculation
    i * Nmesh + j can overflow if Nmesh >= 46341.

    Note: This check is NOT affected by MPI decomposition because
    the loop over j goes from 0 to Nmesh-1 (global size).

    Returns:
        dict with overflow status and details
    """
    # Check if nc >= sqrt(INT_MAX) ≈ 46340
    overflow = nc > INT_MAX_SQRT

    return {
        "module": "LPT Grid Indexing (lpt.c)",
        "nc": nc,
        "overflow": overflow,
        "safe_nc": INT_MAX_SQRT,
        "description": f"Nmesh = nc = {nc:,}",
    }


def print_report(nc, pm_nc_factor, mpi_proc, np_alloc_factor=1.25):
    """Print formatted overflow check report."""
    print("=" * 70)
    print("COLA Halo OMP Overflow Risk Check")
    print("=" * 70)
    print(f"Configuration:")
    print(f"  nc              = {nc:,}")
    print(f"  pm_nc_factor    = {pm_nc_factor}")
    print(f"  mpi_proc        = {mpi_proc}")
    print(f"  np_alloc_factor = {np_alloc_factor}")
    print(f"  INT_MAX         = {INT_MAX:,}")
    print(f"  sqrt(INT_MAX)   ≈ {INT_MAX_SQRT:,}")
    print("-" * 70)

    # Check each module
    checks = [
        check_particle_allocation(nc, mpi_proc, np_alloc_factor),
        check_pm_indexing(nc, pm_nc_factor),
        check_lpt_indexing(nc),
    ]

    all_safe = True
    for check in checks:
        status = "✓ SAFE" if not check["overflow"] else "✗ OVERFLOW RISK"
        if check["overflow"]:
            all_safe = False

        print(f"\n{check['module']}:")
        print(f"  Status: {status}")
        print(f"  Details: {check['description']}")

        if "np_alloc" in check:
            print(f"  np_alloc    = {check['np_alloc']:,}")
            print(f"  local_nx    = {check['local_nx']:,}")
            print(f"  Safe nc     ≈ {check['safe_nc_approx']:,} (for this mpi_proc)")
        elif "ngrid" in check:
            print(f"  Ngrid       = {check['ngrid']:,}")
            print(f"  Safe nc     ≤ {check['safe_nc']:,} (for pm_nc_factor={pm_nc_factor})")
        elif "nc" in check:
            print(f"  Safe nc     ≤ {check['safe_nc']:,}")

    print("\n" + "=" * 70)
    if all_safe:
        print("Overall: ✓ All checks PASSED - No overflow risk detected")
    else:
        print("Overall: ✗ OVERFLOW RISK DETECTED - Some checks FAILED")
        print("\nRecommendations:")
        print("  - Increase MPI processes to reduce particle allocation")
        print("  - Reduce pm_nc_factor if PM indexing overflows")
        print("  - For LPT indexing overflow, code modification is needed")
        print("    (change loop variables from int to long long)")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Check for potential integer overflow risks in COLA Halo OMP"
    )
    parser.add_argument("mpi_proc", type=int, help="Number of MPI processes")
    parser.add_argument("nc", type=int, help="Coarse grid number (e.g., 1024)")
    parser.add_argument(
        "pm_nc_factor", type=int, help="PM grid factor (Ngrid = nc * pm_nc_factor)"
    )
    parser.add_argument(
        "np_alloc_factor",
        type=float,
        nargs="?",
        default=1.25,
        help="Particle allocation factor (default: 1.25)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    if args.mpi_proc <= 0 or args.nc <= 0 or args.pm_nc_factor <= 0 or args.np_alloc_factor <= 0:
        print("Error: All parameters must be positive", file=sys.stderr)
        sys.exit(1)

    print_report(args.nc, args.pm_nc_factor, args.mpi_proc, args.np_alloc_factor)


if __name__ == "__main__":
    main()
