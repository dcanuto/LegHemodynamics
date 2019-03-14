module LegHemodynamics

importall NumericalIntegration
importall MAT
importall Dierckx
importall CSV
importall Missings
importall Distributions

# conversion factors and solver parameters
include("constants.jl")
include("solverparams.jl")

# file loader
include("loadtexttree.jl")

# timer
include("buildtimer.jl")

# type definitions
include("buildbranches.jl")
include("buildall.jl")

# initialization/discretization
include("calcbranchprops.jl")
include("discretizebranches.jl")
include("assignterminals.jl")
include("discretizeperiphery.jl")

# initial conditions for each submodel
include("applybranchics.jl")
include("applyperipheryics.jl")
include("updatediscretization.jl")

# cardiac cycle counter
include("setnumbeats.jl")

# TVD RK3 + WENO3 for artery interiors
include("tvdrk3.jl")
include("invariants.jl")
include("invariantodes.jl")
include("weno3.jl")
include("smoothinds.jl")
include("weights.jl")
include("reconstruct.jl")
include("F1d.jl")
include("J1d.jl")
include("abseigs.jl")

# 0D-1D coupling and 0D updates
include("interiorbcs.jl")
include("coupleproximal.jl")
include("coupledistal.jl")
include("coupling.jl")
include("fdist.jl")
include("Jdist.jl")
include("newtondist.jl")
include("linedist.jl")

# 1D interior junction updates
include("solvesplits.jl")
include("newton.jl")
include("fsingle.jl")
include("fdouble.jl")
include("ftriple.jl")
include("Jsingle.jl")
include("Jdouble.jl")
include("Jtriple.jl")

# 1D pressure updates
include("arterialpressure.jl")

# data assimilation
include("advancetime.jl")
include("builderrors.jl")
include("paramwalk.jl")
include("applycustomics.jl")
include("rediscretizet.jl")

# main
include("main.jl")

end
