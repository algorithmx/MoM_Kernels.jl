# MoM_Kernels.jl – Dependency Update and Test Status (2025-11-22)

## Summary

- Performed a **moderate dependency update** by relaxing some `[compat]` bounds in `Project.toml` to allow newer minor versions within the same major series for non-MoM dependencies.
- Ran `Pkg.update()` and `Pkg.test()` with the updated compat settings.

## Changes to `[compat]` in `Project.toml`

The `[compat]` section was adjusted from narrow ranges (like `"1.0 - 1.6"`) to more permissive, major-version-based constraints, e.g.:

```toml
[compat]
FLoops = "0.2"
FastGaussQuadrature = "1"
FoldsThreads = "0.1"
GSL = "1"
IncompleteLU = "0.2"
IterativeSolvers = "0.9"
JLD2 = "0.5"
LegendrePolynomials = "0.4"
LinearMaps = "3"
MoM_Basics = "0.0.8 - 0.9"   # Kept conservative because it is local & tightly coupled
OffsetArrays = "1"
Primes = "0.5"
ProgressMeter = "1"
SpecialFunctions = "2"
StaticArrays = "1"
ThreadsX = "0.1"
UnicodePlots = "3"
julia = "1.8 - 1.20"
```

MoM-related packages are already used as local path dependencies:

- `MoM_Basics` → `path = "../MoM_Basics.jl"`
- `IterativeSolvers` → `path = "../IterativeSolvers.jl"`

No further action was required to switch them from remote to local; they are already local.

## Packages that actually updated

After changing compat, the following were upgraded via

```julia
using Pkg
Pkg.update()
```

- `ProgressMeter`: **v1.10.4 ⇒ v1.11.0**
- `SpecialFunctions`: **v2.5.1 ⇒ v2.6.1**

No other dependencies changed version.

## Test results after the update

Command:

```julia
using Pkg
Pkg.test()
```

Result:

- **55 passed, 0 failed, 10 errored, 0 broken**

All the errors share the same root cause:

- Exception: `BoundsError: attempt to access 1-element Vector{Vector{ComplexF32}} at index [2]`
- First thrown at:
  - `src/MLFMA/Precondition/SAI.jl:87`
- Involved function:
  - `sparseApproximateInversePl` (called in a threaded context; stack traces show `threadingconstructs.jl` in the call chain).
- Affected testsets (examples):
  - `Triangle, RWG` – `test/runtests.jl:14`
  - `PWC, SWG` – `test/runtests.jl:37` and a second occurrence
  - `PWC, RBF` – `test/runtests.jl:62` and a second occurrence
  - `Tetra + Hexadron, PWC` – `test/runtests.jl:80`
  - `RWG + SWG, RWG + PWC` – `test/runtests.jl:108` and `:134` (multiple errors across these groups)

These errors were present **before** the compat relaxation and remain **after** the update; they appear to be an internal logic issue in the SAI preconditioner rather than a regression caused by newer dependency versions.

## Suggestions for future bug fixing

1. **Inspect `src/MLFMA/Precondition/SAI.jl` around line 87**
   - The `BoundsError` indicates an attempt to access index `2` of a `Vector{Vector{ComplexF32}}` that currently has length `1`.
   - Likely causes:
     - Off-by-one indexing in a loop over sub-blocks or panels.
     - A mismatch between assumed block count and actual allocation (e.g. expecting 2 blocks, only 1 allocated).
     - Threaded chunking logic that uses global indices but allocates per-thread arrays with fewer entries.

2. **Add defensive checks / assertions**
   - Before indexing that vector, assert something like:
     - `@assert length(vec) ≥ expected_index`
   - If the assertion fails, log the indices and sizes to understand which mesh / testsetting triggers the mismatch.

3. **Reproduce with a minimal case**
   - Run an individual failing test (e.g. `Triangle, RWG`) and, if needed, reduce the problem size in the test so you can more easily inspect internal arrays.
   - Optionally run with a single thread (`JULIA_NUM_THREADS=1`) to see if the bug is purely threading-related or also present in serial.

4. **Review parallelization around `sparseApproximateInversePl`**
   - Check how work is divided among threads and how any per-thread arrays/vectors are allocated and indexed.
   - Make sure per-thread containers and global indices are consistent (e.g. no use of global index `i` to index a per-thread vector that only spans a subset).

5. **Re-run tests after fixes**
   - Once `SAI.jl` is corrected, re-run `Pkg.test()` to confirm that the 10 errored testsets are resolved and that the moderate dependency update remains stable.

*This note was generated to document the current dependency state, test status, and a suggested roadmap for fixing the remaining preconditioner bug.*
