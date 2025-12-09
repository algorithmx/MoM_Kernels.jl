## This file is used to build solvers
"""
Iterative solver selection
"""
function iterSolverSet(solverType::Symbol)::Function
    getproperty(IterativeSolvers, solverType)
end

"""
Save current coefficients
"""
function saveCurrent(ICurrent; str = "")
    jldsave(SimulationParams.resultDir*"ICurrent$str.jld2"; ICurrent = ICurrent)
end

"""
Load current coefficients
"""
function loadCurrent(filename)
    occursin("jld2", filename) ? load("$filename", "ICurrent") : load("$filename.jld2", "ICurrent")
end


"""
Matrix equation
Ax=b
Composite solver function
Input:
A::LinearMapType{T}, b::Vector{T}
solverT::Symbol  Solver type
"""
function solve(A::LinearMapType{T}, b::AbstractVector{T};
    solverT::Symbol=:gmres, Pl = Identity(), Pr = Identity(), rtol = 1e-3,
    maxiter = 1000, str = "", restart = 200, save_memtime = true, log = true, 
    verbose = true, args...) where{T<:Number}
    # Direct solve
    solverT  == :direct  && begin
        println("Solving matrix function with LUD.")
        try
            @clock "Direct Solver" begin
                x = A\b
            end
            save_memtime && restore_infos()
            return x, nothing
        catch
            println("Check if input is a matrix!")
            return nothing, nothing
        end
    end  

    FT = real(T)
    # Iterative solver
    solver      =   iterSolverSet(solverT)
    # Residual threshold
    resnorm0    =   FT(norm( Pl \ b))
    resnormtol  =   FT(rtol*resnorm0)

    # Iterative solve
    verbose && @info "\nSolving with $solverT, initial resnorm: $resnorm0.\n"
    @clock "Iterative Solver" begin
        x, ch       =   solver(A, b; restart = restart, abstol = resnormtol, Pl = Pl, Pr = Pr,  log = log, verbose=verbose, maxiter = maxiter, args...)
    end
    save_memtime && restore_infos()
    saveCurrent(x; str = str)
    # Relative residual results
    relresnorm  =   ch.data[:resnorm] / resnorm0

    # Command line plotting
    (verbose && SimulationParams.SHOWIMAGE)  &&  convergencePlot(relresnorm)

    # Write relative residual to file
    open(joinpath(SimulationParams.resultDir, "$(solverT)_ch$str.txt"), "w") do io
        for resi in relresnorm
            write(io, "$resi\n" )
        end
    end

    return x, ch
end

"""
Matrix equation
Ax=b
Composite solver function
Input:
A::LinearMapType{T}, b::Vector{T}
solverT::Symbol  Solver type
"""
function solve!(A::LinearMapType{T}, x::AbstractVector{T}, b::AbstractVector{T}; 
    solverT::Symbol = :gmres!, Pl = Identity(), Pr = Identity(), rtol = 1e-3, 
    maxiter = 1000, str = "", restart = 200, save_memtime = true, log = true, 
    verbose = true, args...) where{T<:Number}
    # Direct solve
    solverT  == :direct  && begin
        println("Solving matrix function with LUD.")
        try
            @clock "Direct Solver" begin
                copyto!(x, A\b)
            end
            save_memtime && restore_infos()
            return (x, nothing)
        catch
            println("Check if input is a matrix!")
            return nothing, nothing
        end
    end  

    FT = real(T)
    # Iterative solver
    solver      =   iterSolverSet(solverT)
    # Residual threshold
    resnorm0    =   FT(norm( Pl \ b))
    resnormtol  =   FT(rtol*resnorm0)

    # Iterative solve
    verbose && @info "\nSolving with $solverT, initial resnorm: $resnorm0.\n"
    @clock "Iterative Solver" begin
        x, ch       =   solver(x, A, b; restart = restart, abstol = resnormtol, Pl = Pl, Pr = Pr,  log = log, verbose = verbose, maxiter = maxiter, args...)
    end

    save_memtime && restore_infos()
    saveCurrent(x; str = str)

    # Relative residual results
    relresnorm  =   ch.data[:resnorm] / resnorm0

    # Command line plotting
    (verbose && SimulationParams.SHOWIMAGE)  &&  convergencePlot(relresnorm)

    # Write relative residual to file
    open(joinpath(SimulationParams.resultDir, "$(solverT)_ch$str.txt"), "w") do io
        for resi in relresnorm
            write(io, "$resi\n" )
        end
    end

    return x, ch

end

"""
Plot convergence curve after calculation
"""
function convergencePlot(resnorm::Vector{FT}) where{FT<:Real}

    figCnvergence = lineplot(resnorm, yscale = :log10, ylim = (minimum(resnorm), 1), xlabel = "Epoch", title = "Relative ResNorm - Epoch")

    SimulationParams.SHOWIMAGE  &&  display(figCnvergence)

    return 

end