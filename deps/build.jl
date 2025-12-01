using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

Pkg.add(path = joinpath(@__DIR__, "..", "..", "IterativeSolvers.jl"))
Pkg.add(path = joinpath(@__DIR__, "..", "..", "MoM_Basics.jl"))
