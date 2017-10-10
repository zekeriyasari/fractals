"""
IteratedFunctionSystem module for analysis of fractals.
"""
module IteratedFunctionSystem

# Imports
import Base.show
import Base.Sys:CPU_CORES
using StatsBase:sample, Weights

# Exports
export Space, RealSpace
export Transformation, Contraction
export IFS
export deterministic_algorithm, random_algorithm

# Define abstract types
abstract type Space{T<:Number} end
abstract type Transformation end


"""
    RealSpace{T<:Real}(set::Type{Vector{T}}, dim::Integer)

Space consisting of real vector `set` of length `dim`.
"""
struct RealSpace{T<:Real} <: Space{T}
    set::Type{Vector{T}}
    dim::Integer
    function RealSpace{T}(set, dim) where {T<:Real}
        if dim < 1
            throw(ArgumentError("Dimension must be nonnegative integer"))
        end
        new(set, dim)
    end
end
# Constructor rule out explicit type declaration.
RealSpace(set::Type{Vector{T}}, dim::Integer) where {T<:Real} =
    RealSpace{T}(set::Type{Vector{T}}, dim::Integer)


"""
    Contraction(rule::Function, domain::Space, range::Space, contractivity::Real)

Contractivity transformation from `domain` to `range` given by the `rule` with
contractivity factor of `contractivity` which is a positive real less than 1.
"""
struct Contraction <: Transformation
    rule::Function
    domain::Space
    range::Space
    contractivity::Real
    function Contraction(rule, domain, range, contractivity)
        if 0 < contractivity < 1
            new(rule, domain, range, contractivity)
        else
            throw(ArgumentError("contractivity  must be real between 0 and 1"))
        end
    end
end


"""
    IFS(contractions::Vector{Contraction}, contractivity::Real)

Iterated function system consiting of `contractions` transformation with
`contractivity` which is maximum of contractivities of `contractions`.
"""
struct IFS
    contractions::Vector{Contraction}
    probabilities::Vector{Real}
    contractivity::Real
    function IFS(contractions, probabilities, contractivity)
        if contractivity > 1 || contractivity < 0
            error("Contractivity must be positive real between 0 and 1")
        end
        if sum(probabilities) != 1
            error("Probabilities must sum up to 1")
        end
        new(contractions, probabilities, contractivity)
    end
end
# Constructor to rule out explicit `contractivity` declarations.
function IFS(contractions)
    contractivities = [contractions[k].contractivity for k = 1 : length(contractions)]
    IFS(contractions, ones(length(contractions)) / length(contractions), maximum(contractivities))
end
# Constructor to convert row arrays to column arrays
IFS(contractions::Array{Contraction, 2}) =
    IFS(reshape(contractions, length(contractions), ))


"""
    deterministic_worker(func::Function, chunk::Array{<:Real, 2})

Deterministic worker function that applies `func` to a `chunk` of data.
"""
function deterministic_worker(func::Function, chunk::Array{<:Real, 2})
    return mapslices(func, chunk, 1)
end


"""
    deterministic_algorithm(ifs::IFS,
                            initial::Union{Void, Vector{<:Real}=nothing,
                            num_steps::Integer=15,
                            num_track_points::Integer=2^10,
                            num_update::Union{Void, Integer}=nothing,
                            monitor::Union{Void, Function}=nothing)

Deterministic algorithm to compute the attractor of the `IFS`.

Deterministic algorithm to compute the attractor of the `IFS` starting with
the initial set `initial`. During computation, `num_steps` itearations are taken.
`num_track_points` is the number of points in the set during computation of the
IFS. If the number elements in the set is increased greater than `num_track_points`,
last `num_track_points` of the set is taken and other points are discarded. If `monitor`
is provided, a remote process is launched to plot the set if the length of the is
reached to `num_update`.
"""
function deterministic_algorithm(ifs::IFS;
                                 initial::Union{Void, Vector{<:Real}}=nothing,
                                 num_steps::Int=10,
                                 num_track_points::Int=2^10,
                                 num_update::Union{Void, Int}=nothing,
                                 monitor::Union{Void, Function}=nothing,
                                 multi_proc::Bool=true)
    # Check initial set.
    domain_dim = ifs.contractions[1].domain.dim
    if initial == nothing
        initial = rand(domain_dim, 1)
    else
        if size(initial)[1] != domain_dim
            throw(DimensionMismatch("Initial does not match ifs domain dimension."))
        end
    end

    # If monitor is provided, construct a data channel
    if monitor != nothing
        channel = Channel(num_track_points)  # Construct data channel
        @schedule monitor(channel)  # Launch the remote monitor.
        # Check number of update points.
        if num_update == nothing
            num_update = Int(floor(num_steps / 5))
        end
    end

    # Compute the attractor
    num_procs = CPU_CORES
    num_contractions = length(ifs.contractions)
    set = initial
    for i = 1 : num_steps
        num_points = size(set)[2]
        next_set = zeros(size(set)[1], size(set)[2] * num_contractions)
        for j = 1 : num_contractions
            func = ifs.contractions[j].rule
            if multi_proc && num_procs > 2
                # Distribute the set to available processors.
                chunk_size = Int(floor(num_points / num_procs))
                processes = Vector{Future}(num_procs)
                for k = 1 : num_procs
                    chunk = set[:, 1 + (k - 1) * chunk_size : k * chunk_size]
                    processes[k] = @spawn deterministic_worker(func, chunk)
                end
                # Collect the results from the processors.
                temp = zeros(size(set))
                for k = 1 : num_procs
                    temp[:, 1 + (k - 1) * chunk_size : k * chunk_size] = fetch(processes[k])
                end
            else
                temp = deterministic_worker(func, set)
            end
            next_set[:, 1 + (j - 1) * num_points : j * num_points] = temp
        end  # End of for.

        # if size(next_set)[2] > num_track_points
        #     next_set = next_set[1 : Int(ceil(size(next_set)[2] / num_track_points)) : end]
        # end

        # Moitor the points
        if monitor != nothing && mod(i, num_update) == 0
            message = "iter: $i, num_points: $num_points"
            push!(channel, (next_set[1, :], next_set[2, :]))
        end  # End of if

        set = next_set

    end  # End of for

    return set

end  # End of deterministic_algorithm function


# FIXME: num_iter does not passed correctly so data size is not correct. 
function random_worker(funcs_probs::Tuple{Vector{Contraction}, Vector{Float64}},
                       initial::Array{Float64, 2},
                       num_iter::Int64;
                       discard::Bool=false)
    funcs, probs = funcs_probs
    data = zeros(size(initial)[1], num_iter)
    data[:, 1] = initial
    for i = 1 : num_iter - 1
        func = sample(funcs, Weights(probs))    # Choose a random function `func`
        data[:, i + 1] = func.rule(data[:, i])  # Apply `func` to the data.
    end
    println("Data: $data")
    println("Data Size: $(size(data))")
    if discard
        return data[:, end]
    end
    return data
end


function random_algorithm(ifs::IFS;
                          initial::Union{Void, Vector{<:Real}}=nothing,
                          num_steps::Int=2^10,
                          probabilities::Union{Vector{<:Real}, Void}=nothing,
                          num_update::Union{Void, Int}=nothing,
                          monitor::Union{Void, Function}=nothing,
                          multi_proc::Bool=true,
                          num_transients::Int=2^7)
    # Check initial set.
    domain_dim = ifs.contractions[1].domain.dim
    if initial == nothing
        initial = rand(domain_dim, 1)
    else
        if size(initial)[1] != domain_dim
            throw(DimensionMismatch("Initial does not match ifs domain dimension."))
        end
    end

    # Check probabilities
    if probabilities != nothing
        sum(probabilities) == 1 || error("Sum of probabilities must be one")
        length(probabilities) == length(ifs.contractions) || error("Length of probabilities does not match ifs.")
    else
        probabilities = ifs.probabilities
    end

    # If monitor is provided, construct a data channel
    if monitor != nothing
        channel = Channel(num_steps)  # Construct data channel
        @schedule monitor(channel)  # Launch the remote monitor.
        # Check number of update points.
        if num_update == nothing
            num_update = Int(floor(num_steps / 5))
        end
    end

    # Transient iterations
    initial = random_worker((ifs.contractions, probabilities), initial, num_transients, discard=true)
    println("Initial $initial")

    # Compute the attractor
    num_procs = CPU_CORES
    if multi_proc && num_procs > 2  # Use multiple
        iters = Int(floor(num_steps / num_procs)) * ones(Int, num_procs)
        if (num_steps % num_procs) != 0
            iters = [iters; Int(num_steps % num_procs)]
            num_procs += 1
        end
        println("Iters: $iters")
        println("num_procs: $num_procs")
        # Launch the remote processes.
        processes = Vector{Future}(num_procs)
        for i = 1 : num_procs
            processes[i] = @spawn random_worker((ifs.contractions, probabilities), initial, iters[i], discard=false)
        end
        # Collect the results from the remote processes.
        set = zeros(size(initial)[1], 1)
        for i = 1 : num_procs
            value = fetch(processes[i])
            println("Set: $set, value: $value")
            set = [set value]
        end
        set = set[:, 2:end]
    else  # Use single processing
        set = random_worker((ifs.contractions, probabilities), initial, num_steps, discard=false)
    end  # End of if-else

    return set

end  # End of deterministic_algorithm function

end # End of IteratedFunctionSystem
