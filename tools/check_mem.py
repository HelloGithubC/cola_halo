#!/usr/bin/env python3
"""
Memory estimator for COLA Halo OMP based on configuration parameters.

Calculates peak memory requirements for particle data, shared memory regions
(mem1 and mem2), and provides total memory estimates with optional per-MPI
process averaging.
"""

import argparse
import sys


def format_memory(bytes_val):
    """Format memory size with appropriate unit (MB or GB)."""
    mb = bytes_val / (1024 ** 2)
    if mb < 1024:
        return f"{mb:.2f} MB"
    gb = bytes_val / (1024 ** 3)
    return f"{gb:.2f} GB"


def calculate_memory(nc, nc_pm_factor, np_alloc_factor, nproc=1):
    """
    Calculate memory requirements for COLA Halo OMP.

    Args:
        nc: Coarse grid number (e.g., 1024)
        nc_pm_factor: PM grid factor (Ngrid = nc * nc_pm_factor)
        np_alloc_factor: Particle allocation factor
        nproc: Number of MPI processes (for averaging, default=1)

    Returns:
        dict with memory breakdown in bytes
    """
    # Derived parameters
    ngrid = nc * nc_pm_factor  # PM grid resolution

    # Particle allocation: np_alloc = np_alloc_factor * nc^2 * (nc + 1)
    # This accounts for the x-direction decomposition (local_nx + 1) ≈ (nc + 1)
    np_alloc = np_alloc_factor * (nc ** 2) * (nc + 1)

    # Constants (bytes)
    sizeof_float = 4
    sizeof_complex = 8  # fftwf_complex = 2 * float
    sizeof_particle = 60  # x[3]+dx1[3]+dx2[3]+v[3]+id
    sizeof_force = 12  # float3

    # 1. Particle data (independent allocation)
    #    particles->p[np_alloc]: 60 bytes/particle
    #    particles->force[np_alloc]: 12 bytes/particle
    particle_mem = (sizeof_particle + sizeof_force) * np_alloc

    # 2. Shared memory mem1 (reuse: LPT / PM / FoF, take maximum)
    #
    # 2.1 LPT grid: 12 complex arrays, each size ≈ nc^3 / 2 (compressed)
    size_lpt_one = (nc ** 3) // 2
    mem_lpt = 12 * size_lpt_one * sizeof_complex

    # 2.2 PM grid: density grid in k-space, size ≈ Ngrid^3 / 2
    size_pm_one = (ngrid ** 3) // 2
    mem_pm = size_pm_one * sizeof_complex

    # 2.3 FoF memory (rough estimate)
    #    KD-tree nodes ≈ 2 * np_alloc, each ~40 bytes
    #    Plus additional buffers ~12 bytes/particle * np_alloc
    mem_fof = 40 * 2 * np_alloc + 12 * np_alloc

    # mem1 = max(LPT, PM, FoF)
    mem1 = max(mem_lpt, mem_pm, mem_fof)

    # 3. Shared memory mem2 (reuse: density_k / snapshot, take maximum)
    #
    # 3.1 density_k: (Ngrid/2+1) * Ngrid * Ngrid * 8 bytes
    #    (complex array with Hermitian symmetry)
    size_density_k = (ngrid // 2 + 1) * ngrid * ngrid
    mem_density_k = size_density_k * sizeof_complex

    # 3.2 snapshot: 28 bytes/particle * np_alloc
    mem_snapshot = 28 * np_alloc

    # mem2 = max(density_k, snapshot)
    mem2 = max(mem_density_k, mem_snapshot)

    # 4. PM buffer (negligible, ~30 MB for typical configs)
    #    Omitted as it's small compared to other components

    # Total peak memory (mem1 and mem2 are reused, but peak is their sum
    # since they may be simultaneously allocated during transitions)
    total_mem = particle_mem + mem1 + mem2

    # Per-MPI process memory (distributed across nproc processes)
    per_proc_mem = total_mem / nproc

    return {
        "nc": nc,
        "nc_pm_factor": nc_pm_factor,
        "ngrid": ngrid,
        "np_alloc": int(np_alloc),
        "particle_mem": particle_mem,
        "mem_lpt": mem_lpt,
        "mem_pm": mem_pm,
        "mem_fof": mem_fof,
        "mem1": mem1,
        "mem_density_k": mem_density_k,
        "mem_snapshot": mem_snapshot,
        "mem2": mem2,
        "total_mem": total_mem,
        "nproc": nproc,
        "per_proc_mem": per_proc_mem,
    }


def print_report(mem):
    """Print formatted memory report."""
    print("=" * 60)
    print("COLA Halo OMP Memory Estimation")
    print("=" * 60)
    print(f"Configuration:")
    print(f"  nc              = {mem['nc']}")
    print(f"  nc_pm_factor    = {mem['nc_pm_factor']}")
    print(f"  Ngrid (PM)      = {mem['ngrid']}")
    print(f"  np_alloc_factor = (derived from np_alloc)")
    print(f"  np_alloc        = {mem['np_alloc']:,}")
    print(f"  nproc           = {mem['nproc']}")
    print("-" * 60)

    print("Memory Breakdown:")
    print(f"  Particle data   : {format_memory(mem['particle_mem'])}")
    print(f"    - particles   : {format_memory(60 * mem['np_alloc'])}")
    print(f"    - force       : {format_memory(12 * mem['np_alloc'])}")
    print()
    print(f"  mem1 (max of below): {format_memory(mem['mem1'])}")
    print(f"    - LPT grid    : {format_memory(mem['mem_lpt'])}")
    print(f"    - PM grid     : {format_memory(mem['mem_pm'])}")
    print(f"    - FoF         : {format_memory(mem['mem_fof'])}")
    print()
    print(f"  mem2 (max of below): {format_memory(mem['mem2'])}")
    print(f"    - density_k   : {format_memory(mem['mem_density_k'])}")
    print(f"    - snapshot    : {format_memory(mem['mem_snapshot'])}")
    print("-" * 60)

    print(f"Total peak memory        : {format_memory(mem['total_mem'])}")
    if mem["nproc"] > 1:
        print(f"Per MPI process (average): {format_memory(mem['per_proc_mem'])}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Estimate memory requirements for COLA Halo OMP"
    )
    parser.add_argument("nc", type=int, help="Coarse grid number (e.g., 1024)")
    parser.add_argument(
        "nc_pm_factor", type=int, help="PM grid factor (Ngrid = nc * nc_pm_factor)"
    )
    parser.add_argument(
        "np_alloc_factor",
        type=float,
        help="Particle allocation factor (typically 1.0)",
    )
    parser.add_argument(
        "-n",
        "--nproc",
        type=int,
        default=1,
        help="Number of MPI processes (default: 1). "
        "Memory is averaged per process when > 1.",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    if args.nc <= 0 or args.nc_pm_factor <= 0 or args.np_alloc_factor <= 0:
        print("Error: All parameters must be positive", file=sys.stderr)
        sys.exit(1)

    if args.nproc <= 0:
        print("Error: nproc must be positive", file=sys.stderr)
        sys.exit(1)

    mem = calculate_memory(args.nc, args.nc_pm_factor, args.np_alloc_factor, args.nproc)
    print_report(mem)


if __name__ == "__main__":
    main()
