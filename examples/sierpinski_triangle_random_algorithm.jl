workspace()
include("../src/ifs.jl")
include("../src/monitors.jl")
using IFS
using Base.Test
using Monitors
using PyPlot

# Sierpinski triangle ifs code
f1(x) = [0.5 0; 0 0.5] * x + [1; 1]
f2(x) = [0.5 0; 0 0.5] * x + [1; 50]
f3(x) = [0.5 0; 0 0.5] * x + [50; 50]
x0 = rand(2)
probabilities = [0.33; 0.33; 0.34]
attractor = random_algorithm([f1, f2, f3], probabilities, x0, 10000, multi_procs=false)
plot(attractor[1, :], attractor[2, :], ".", linewidth=0.1)
