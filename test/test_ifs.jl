workspace()

include("../src/ifs.jl")
include("./unit_test.jl")

using IteratedFunctionSystem
using UnitTest
using PyPlot

domain = RealSpace(Vector{Real}, 1)
range = RealSpace(Vector{Real}, 1)
trfm1 = Contraction(x -> 1/2 * x, domain, range, 1/2)
trfm2 = Contraction(x -> 1/2 * x + 1/2, domain, range, 1/2)
ifs = IFS([trfm1, trfm2])
res = deterministic_algorithm(ifs)
stem(res)
