workspace()

include("../src/ifs.jl")
include("./unit_test.jl")
include("../src/monitors.jl")

using IteratedFunctionSystem
using UnitTest
using Monitors
using PyPlot

# domain = RealSpace(Vector{Real}, 1)
# range = RealSpace(Vector{Real}, 1)
# trfm1 = Contraction(x -> 1/2 * x, domain, range, 1/2)
# trfm2 = Contraction(x -> 1/2 * x + 1/2, domain, range, 1/2)
# ifs = IFS([trfm1, trfm2])

domain = RealSpace(Vector{Real}, 2)
range = RealSpace(Vector{Real}, 2)
trfm1 = Contraction(x ->  [0.5 0.; 0. 0.5] * x + [1; 1], domain, range, 0.5)
trfm2 = Contraction(x ->  [0.5 0.; 0. 0.5] * x + [1; 50], domain, range, 0.5)
trfm3 = Contraction(x ->  [0.5 0.; 0. 0.5] * x + [25; 50], domain, range, 0.5)
ifs = IFS([trfm1, trfm2, trfm3])
res = deterministic_algorithm(ifs, monitor=update_plot, num_steps=20, num_track_points=2^12)
# plot(res[1, :], res[2, :], ".")
