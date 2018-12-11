module LegModule

importall NumericalIntegration
importall MAT
importall Interpolations
importall CSV
importall Missings

# conversion factors and solver parameters
include("conversions.jl")
include("solverparams.jl")

# timer
include("buildtimer.jl")

# file loaders/solution struct builders
include("buildall.jl")
include("buildbranches.jl")
include("calcbranchprops.jl")
include("discretizebranches.jl")
include("discretizeperiphery.jl")

# blood volume tracker
include("updatevolumes.jl")

# main
include("main.jl")



end
