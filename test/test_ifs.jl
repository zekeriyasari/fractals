workspace()

include("../src/ifs.jl")
include("./unit_test.jl")
include("../src/monitors.jl")

using IteratedFunctionSystem
using UnitTest
using Monitors

println("\nTest: RealSpace")
domain = RealSpace(Vector{Real}, 2)
equal(domain.dim, 2) ? pass() : fail("Expected 2, got $(domain.dim), instead")
range = RealSpace(Vector{Real}, 2)
equal(range.dim, 2) ? pass() : fail("Expected 2, got $(range.dim), instead")
println("Test: RealSpace passed.")

println("\nTest: Contraction started...")
f1(x) = [0.5 0.; 0. 0.5] * x + [1; 1]
trfm1 = Contraction(f1, domain, range, 0.5)
equal(trfm1.domain, domain) ? pass() : fail("Expected $domain, got $(trfm1.domain) instead.")
equal(trfm1.range, range) ? pass() : fail("Expected $range, got $(trfm1.range) instead.")
data = rand(2, 10)
for i = 1 : size(data)[2]
    expected = f1(data[:, i])
    calculated = trfm1.rule(data[:, i])
    equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated, instead")
end
println("Test: Contraction passed.")

println("Test: IFS started...")
# # Define contractions for Sierpinski triangle
# trfm1 = Contraction(x ->  [0.5 0.; 0. 0.5] * x + [1; 1], domain, range, 0.5)
# trfm2 = Contraction(x ->  [0.5 0.; 0. 0.5] * x + [1; 50], domain, range, 0.5)
# trfm3 = Contraction(x ->  [0.5 0.; 0. 0.5] * x + [25; 50], domain, range, 0.5)
# ifs = IFS([trfm1, trfm2, trfm3])
# res = deterministic_algorithm(ifs, monitor=update_plot, num_steps=20, num_track_points=2^10)
# Define contractions for Sierpinski triangle
trfm1 = Contraction(x ->  [0. 0.; 0. 0.16] * x + [0; 0], domain, range, 0.5)
trfm2 = Contraction(x ->  [0.85 0.04; -0.04 0.85] * x + [0; 1.6], domain, range, 0.5)
trfm3 = Contraction(x ->  [0.2 -0.26; 0.23 0.22] * x + [0; 1.6], domain, range, 0.5)
trfm4 = Contraction(x ->  [-0.15 0.28; 0.26 0.24] * x + [0; 0.44], domain, range, 0.5)
ifs = IFS([trfm1, trfm2, trfm3, trfm4])
res = deterministic_algorithm(ifs, monitor=update_plot, num_steps=8, num_track_points=2^10, multi_proc=true)
println("Test: IFS passed. ")
