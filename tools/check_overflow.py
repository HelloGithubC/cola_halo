#!/usr/bin/env python3
"""
Overflow checker for COLA Halo OMP based on configuration parameters.

Checks for potential integer overflow risks in:
1. Particle allocation cast to int (mem.c:69)
2. PM grid loop indices (pm.c) - int loops with size_t index
3. LPT grid loop indices (lpt.c) - int loops
4. Memory allocation size_t overflow (64-bit, safe for practical configs)

All checks consider the effect of MPI parallelization on local work distribution.
"""

import argparse
import sys


# Constants
INT_MAX = 2**31 - 1  # Maximum value for 32-bit signed integer
INT_MAX_SQRT = int(INT_MAX**0.5)  # ~46340
SIZE_T_MAX = 2**64 - 1  # Maximum value for 64-bit unsigned integer (size_t)


def check_particle_allocation(nc, mpi_proc, np_alloc_factor=1.25):
    """
    Check overflow risk in particle allocation (mem.c:69).

    The allocation formula is:
        np_alloc = (int)(np_alloc_factor * nc^2 * (local_nx + 1))

    where local_nx ≈ nc / mpi_proc (MPI decomposition in x-direction).

    The result is cast to int, so overflow occurs when exceeding INT_MAX.

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
    Check overflow risk in PM grid loop indices (pm.c).

    The PM grid uses Ngrid = nc * pm_nc_factor.
    Loop variables (Jl, iI, K) are declared as int:
        for(int Jl=0; Jl<Local_ny_td; Jl++)
          for(int iI=0; iI<Ngrid; iI++)
            for(int K=0; K<Ngrid/2+1; K++)

    However, the index calculation uses size_t (64-bit):
        size_t index= K + (NgridL/2+1)*(iI + NgridL*Jl);

    The risk is that intermediate calculations with int variables
    (e.g., iI + NgridL*Jl) could overflow if Ngrid >= 46341,
    even though the final index is size_t.

    Returns:
        dict with overflow status and details
    """
    ngrid = nc * pm_nc_factor

    # Check if Ngrid >= sqrt(INT_MAX) ≈ 46340
    # When iI and Jl loop to Ngrid-1, intermediate int calculations can overflow
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
    Check overflow risk in LPT grid loop indices (lpt.c).

    The LPT grid uses Nmesh = nc.
    Loop variables are declared as int:
        for(int i=0; i<Nmesh/2; i++)
          for(int j=0; j<i; j++)
            seedtable[i * Nmesh + j] = ...

    The index calculation i * Nmesh + j can overflow if Nmesh >= 46341,
    since i and j are int and Nmesh is int.

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


def check_memory_allocation(nc, pm_nc_factor, mpi_proc):
    """
    Check if total memory allocation would exceed size_t (64-bit).

    This is a sanity check - for practical configurations, size_t (18 EB)
    is far beyond any realistic memory需求.

    Returns:
        dict with overflow status and details
    """
    import math

    # Rough memory estimation (bytes)
    ngrid = nc * pm_nc_factor
    local_nx = math.ceil(nc / mpi_proc)

    # Particle memory (72 bytes/particle)
    np_alloc = 1.25 * nc * nc * (local_nx + 1)
    particle_mem = 72 * np_alloc

    # LPT grid memory (12 complex arrays)
    lpt_mem = 12 * (nc**3 // 2) * 8

    # PM grid memory
    pm_mem = (ngrid**3 // 2) * 8

    total_mem = particle_mem + lpt_mem + pm_mem

    # Check if total memory exceeds size_t (practically impossible)
    overflow = total_mem > SIZE_T_MAX

    return {
        "module": "Memory Allocation (size_t)",
        "total_mem_bytes": total_mem,
        "overflow": overflow,
        "description": f"Estimated total memory ≈ {total_mem / (1024**3):.1f} GB",
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
    print(f"  SIZE_T_MAX      = {SIZE_T_MAX:,} (64-bit)")
    print("-" * 70)

    # Check each module
    checks = [
        check_particle_allocation(nc, mpi_proc, np_alloc_factor),
        check_pm_indexing(nc, pm_nc_factor),
        check_lpt_indexing(nc),
        check_memory_allocation(nc, pm_nc_factor, mpi_proc),
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
        elif "total_mem_bytes" in check:
            print(f"  Safe for all practical configurations")

    print("\n" + "=" * 70)
    if all_safe:
        print("Overall: ✓ All checks PASSED - No overflow risk detected")
    else:
        print("Overall: ✗ OVERFLOW RISK DETECTED - Some checks FAILED")
        print("\nRecommendations:")
        print("  - Increase MPI processes to reduce particle allocation")
        print("  - Reduce pm_nc_factor if PM indexing overflows")
        print("  - For LPT/PM indexing overflow, code modification is needed")
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
