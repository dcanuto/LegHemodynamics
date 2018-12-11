module LegHemodynamics

importall NumericalIntegration
importall MAT
importall Interpolations
importall CSV
importall Missings

# conversion factors and solver parameters
include("conversions.jl")
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
# include("discretizebranches.jl")
# include("discretizeperiphery.jl")

# blood volume tracker
# include("updatevolumes.jl")

# main
include("main.jl")

# export LegSystem
# export CVTimer
# export ArterialBranches
# export SolverParams

end
