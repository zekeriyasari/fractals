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
# Add method to `Base.show` to display `RealSpace` types.
# Base.show(io::IO, space::RealSpace) = print(io, "RealSpace\n
#                                                  set: $(space.set)\n
#                                                  dim: $(space.dim)")


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
# Add method to `Base.show` to display `Contraction` types.
# Base.show(io::IO, contraction::Contraction) =
#     print(io, "Contraction\n
#                rule: $(contraction.rule)\n
#                domain: $(contraction.domain)\n
#                range: $(contraction.range)\n
#                contractivity: $(contraction.contractivity)")


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
# # Add method `Base.show` to display `IFS` types.
# Base.show(io::IO, ifs::IFS) = print(io, "IFS\n
#                                          contractions: $(ifs.contractions)\n
#                                          contractivity: $(ifs.contractivity)")


"""
    deterministic_algorithm(ifs::IFS,
                            initial::Union{Void, Vector{T} where {T<:Real}}=nothing,
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
                                 initial::Union{Void, Vector{T} where {T<:Real}}=nothing,
                                 num_steps::Integer=15,
                                 num_track_points::Integer=2^10,
                                 num_update::Union{Void, Integer}=nothing,
                                 monitor::Union{Void, Function}=nothing,
                                 multi_proc::Bool=true)

    # Check number of update points.
    if num_update == nothing
        num_update = Int(floor(num_steps / 5))
    end

    # Check initial set.
    domain_dim = ifs.contractions[1].domain.dim
    if initial == nothing
        if multi_proc
            initial = rand(domain_dim, CPU_CORES)
        else
            initial = rand(domain_dim, 1)
        end
    else
        if size(initial)[1] != domain_dim
            throw(DimensionMismatch("Inital does not match ifs domain dimension."))
        end
    end

    # If monitor is provided, construct a data channel
    if monitor != nothing
        channel = Channel(num_track_points)  # Construct data channel
        @schedule monitor(channel)  # Launch the remote monitor.
    end

    # Compute the attractor
    num_procs = CPU_CORES
    num_contractions = length(ifs.contractions)
    set = initial
    for i = 1 : num_steps
        num_points = size(set)[2]
        results = Matrix{Real}(size(set)[1], size(set)[2] * num_contractions)

        if multi_proc && num_procs > 2  # If enabled and possible, use multiple cores.
            for j = 1 : num_contractions
                func = ifs.contractions[j].rule

                # Distribute the set to available processors
                chunk_size = Int(floor(num_points / num_procs))
                processes = Vector{Future}(num_procs)
                for k = 1 : num_procs
                    chunk = set[:, 1 + (k - 1) * chunk_size : k * chunk_size]
                    processes[k] = @spawn mapslices(func, chunk, 1)
                end

                # Collect the results from the processors
                result = Array{Real, 2}(size(set))
                for k = 1 : num_procs
                    value = fetch(processes[k])
                    result[:, 1 + (k - 1) * chunk_size : k * chunk_size] = value
                end
                results[:, 1 + (j - 1) * num_points : j * num_points] = result
            end  # End of for
        else  # If not possible to use multiprocessing, use single core.
            for j = 1 : num_contractions
                func = ifs.contractions[j].rule
                result = mapslices(func, set, 1)
                results[:, 1 + (j - 1) * num_points : j * num_points] = result
            end  # End of for
        end  # End of if-else.

        # # Check the number of track points
        # if size(results)[2] > num_track_points
        #     results = results[:, end - num_track_points + 1 : end]
        # end

        # Moitor the points
        if monitor != nothing && size(set)[1] == 2 && mod(i, num_update) == 0
            message = "iter: $i, num_points: $num_points"
            push!(channel, (results[1, :], results[2, :]))
        end  # End of if

        set = results

    end  # End of for

    return set

end  # End of deterministic_algorithm function

function random_worker(ifs::IFS, initial::Matrix{<:Real}, chunk_size::Integer; probabilities=nothing)
    # Check probabilities
    if probabilities == nothing
        probabilities = ifs.probabilities
    end

    data = Matrix{Real}(ifs.contractions[1].domain.dim, chunk_size + 1)
    data[:, 1] = initial
    for i = 1 : chunk_size
        func = sample(ifs.contractions, Weights(probabilities))
        data[:, i + 1] = func.rule(data[:, i])
    end
    return data
end  # End of random_worker

function random_algorithm(ifs::IFS;
                          initial::Union{Void, Vector{T} where {T<:Real}}=nothing,
                          probabilities::Union{Vector{<:Real}, Void}=nothing,
                          num_steps::Integer=15,
                          num_update::Union{Void, Integer}=nothing,
                          monitor::Union{Void, Function}=nothing,
                          multi_proc::Bool=false)
    # Check the probabilities
    if probabilities != nothing
        sum(probabilities) != 1 && error("Sum of probabilities must be 1")
        length(probabilities) != length(ifs.contractions) && error("Length of probabilities does not match to ifs.contractions.")
    else
        probabilities = ifs.probabilities
    end  # End of if-else.

    # Check initial set.
    domain_dim = ifs.contractions[1].domain.dim
    if initial == nothing
        if multi_proc
            initial = rand(domain_dim, CPU_CORES)
        else
            initial = rand(domain_dim, 1)
        end
    else
        if size(initial)[1] != domain_dim
            throw(DimensionMismatch("Inital does not match ifs domain dimension."))
        end
    end

    # If monitor is provided, construct a data channel
    if monitor != nothing
        channel = Channel(num_track_points)  # Construct data channel
        @schedule monitor(channel)  # Launch the remote monitor.
    end

    # Compute the attractor
    num_procs = CPU_CORES
    attractor = Matrix{Real}(size(initial)[1], num_steps + 1)
    attractor[:, 1] = initial
    if multi_proc && num_procs > 2
        chunk_size = Int(floor(num_steps / num_procs))
    else
        attractor = random_worker(ifs, initial, num_steps, probabilities=probabilities)
    end  # End of if-else.

    return attractor

end  # End of random_algorithm function.

end # End of IteratedFunctionSystem
