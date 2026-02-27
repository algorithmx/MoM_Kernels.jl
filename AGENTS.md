# MoM_Kernels.jl - Computational Engine

Numerical computation and solvers. Consumer of MoM_Basics.jl. Depends on IterativeSolvers.jl.

## Source Tree

```
src/
├── MoM_Kernels.jl          # Main module
├── Inputs.jl               # Parameter input
├── Recorder.jl             # Performance profiling
├── DirectAlgorithm.jl      # Direct LU solver
├── FastAlgorithm.jl        # Fast multipole interface
├── Solvers.jl              # Linear solver wrappers
│
├── ZmatAndVvec/            # Impedance matrix & excitation vector
│   ├── ZmatVvec.jl         # Main assembly interface
│   ├── Singularity.jl      # Singularity handling
│   ├── Singularity/        # Singularity extraction methods
│   ├── EFIE/               # EFIE formulations (12+ variants)
│   ├── MFIE/               # MFIE formulations
│   └── CFIE/               # CFIE formulations
│
├── MLFMA/                  # Multi-Level Fast Multipole Algorithm
│   ├── MLFMAParams.jl      # MLFMA parameters
│   ├── MLFMA.jl            # Core implementation
│   ├── OctreeInfo.jl       # Octree structure
│   ├── LevelInfo.jl        # Level information
│   ├── Znear/              # Near-field matrix computation
│   ├── Precondition/       # Preconditioners (SAI, ILU)
│   └── AggOnBF/            # Basis function aggregation
│
├── PostProcess/            # Post-processing
│   ├── PostProcessing.jl
│   ├── FarField.jl         # Far-field computation
│   ├── RCS.jl              # Radar Cross Section
│   └── CurrentOnGeos.jl    # Surface current extraction
│
└── Extends/                # Extensions
    ├── ParallelParams.jl   # Parallel parameters
    ├── PartitionedVector.jl # Partitioned vectors (large memory)
    └── IOs.jl              # Extended I/O
```

## Key Components

| Module | Purpose |
|--------|---------|
| ZmatAndVvec | MoM matrix/vector assembly with singularity treatment |
| MLFMA | O(N log N) fast solver with octree partitioning |
| PostProcess | Far-field, RCS, current extraction |
| Solvers | Linear system wrappers |
