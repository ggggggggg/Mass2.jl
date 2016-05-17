# eventually these should go inside the module or their own modules and I should choose what to export. not yet
include("LJH.jl")
using LJH
using Dierckx, PyCall, HDF5, ArrayViews, Distributions, KernelDensity, PyPlot, JLD
import Graphs

include("summarize.jl")
include("mockpulses.jl")
include("runningvectors.jl")
include("histogram.jl")
include("steps.jl")
include("steps_macros.jl")
include("calibration.jl")
include("corrections.jl")
include("ljhutil.jl")
include("tes_optimal_filtering.jl")
include("glue.jl")
include("matterwatcher.jl")
module Mass2

# package code goes here

end # module
