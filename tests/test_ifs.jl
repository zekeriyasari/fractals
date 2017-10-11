workspace()

include("../src/ifs.jl")
using IFS
using Base.Test

@test_skip @testset "DeterministicWorker Test" begin
    # 2-dimensional fucntion test.
    A, b = [1 2; 3 4], [5; 6]
    f(x) = A * x + b
    data = [7 8; 9 10]
    @test deterministic_worker(f, data) == [A * data[:, 1] + b A * data[:, 2] + b]
    # 1-dimensional fucntion test.
    A, b = 5, 7
    f(x) = A * x + b
    data = [8]
    @test deterministic_worker(f, data) == A * data + b
    # Check MethodError is thrown
    data = 10  # Not an array.
    @test_throws MethodError deterministic_worker(f, data)
end


@test_skip @testset "DeterministicAlgrithm Test" begin
    # Sierpinski triangle ifs code
    f1(x) = [0.5 0; 0 0.5] * x + [1; 1]
    f2(x) = [0.5 0; 0 0.5] * x + [1; 50]
    f3(x) = [0.5 0; 0 0.5] * x + [50; 50]
    x0 = rand(2)
    @test deterministic_algorithm([f1, f2, f3], x0, 1, multi_procs=false) == [f1(x0) f2(x0) f3(x0)]
    @test deterministic_algorithm([f1, f2, f3], x0, 1, multi_procs=true) == [f1(x0) f2(x0) f3(x0)]
    x0 = rand(2, 1)
    @test deterministic_algorithm([f1, f2, f3], x0, 1, multi_procs=false) == [f1(x0) f2(x0) f3(x0)]
    @test deterministic_algorithm([f1, f2, f3], x0, 1, multi_procs=true) == [f1(x0) f2(x0) f3(x0)]
end


@test_skip @testset "RandomWorker Test" begin
    f1(x) = [0.5 0; 0 0.5] * x + [1; 1]
    f2(x) = [0.5 0; 0 0.5] * x + [1; 50]
    f3(x) = [0.5 0; 0 0.5] * x + [50; 50]
    x0 = rand(2)
    probs = [0.33; 0.33; 0.34]
    @test size(random_worker([f1, f2, f3], probs, x0, 10)) == (2, 10)
    @test size(random_worker([f1, f2, f3], probs, x0, 10, discard=true)) == (2, 1)
end


@testset "RandomAlgorithm Test" begin
    # Sierpinski triangle ifs code
    f1(x) = [0.5 0; 0 0.5] * x + [1; 1]
    f2(x) = [0.5 0; 0 0.5] * x + [1; 50]
    f3(x) = [0.5 0; 0 0.5] * x + [50; 50]
    x0 = rand(2)
    probabilities = [0.33; 0.33; 0.34]
    @test size(random_algorithm([f1, f2, f3], probabilities, x0, 10, multi_procs=false)) == (2, 10)
    @test size(random_algorithm([f1, f2, f3], probabilities, x0, 10, multi_procs=true)) == (2, 10)
    # Plot the results
    using PyPlot
    attractor = random_algorithm([f1, f2, f3], probabilities, x0, 10000, multi_procs=false)
    plot(attractor[1, :], attractor[2, :], ".", linewidth=0.1)
end
