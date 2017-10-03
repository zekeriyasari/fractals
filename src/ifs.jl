
module IteratedFunctionSystem

# TODO: Define Hausforff distance between sets
# TODO: Terminate

# Imports
import Base.show

# Exports
export Space, RealSpace
export Transformation, Contraction
export IFS
export deterministic_algorithm

# Define abstract types
abstract type Space{T<:Number} end
abstract type Transformation end


# Define RealSpace
struct RealSpace{T<:Real} <: Space{T}
    set::Type{Vector{T}}
    dim::Integer
    function RealSpace{T}(set, dim) where {T<:Real}
        if dim < 1
            throw(ArgumentError("Dimension must be nonnegative integer"))
        end
        new(set, dim)
    end  # End of function RealSpace
end  # End of type RealSpace

RealSpace(set::Type{Vector{T}}, dim::Integer) where {T<:Real} =
    RealSpace{T}(set::Type{Vector{T}}, dim::Integer)

Base.show(io::IO, space::RealSpace) = print(io, "RealSpace\n
                                                 \tset: $(space.set)\n
                                                 \tdim: $(space.dim)")


# Define contraction transformations
struct Contraction <: Transformation
    rule::Function
    domain::Space
    range::Space
    contractivity::Number
    function Contraction(rule, domain, range, contractivity)
        if 0 < contractivity < 1
            new(rule, domain, range, contractivity)
        else
            throw(ArgumentError("contractivity  must be real between 0 and 1"))
        end
    end
end

Base.show(io::IO, contraction::Contraction) =
    print(io, "Contraction\n
               \trule: $(contraction.rule)\n
               \tdomain: $(contraction.domain)\n
               \trange: $(contraction.range)\n
               \tcontractivity: $(contraction.contractivity)")


# Define IteratedFunctionSystem
struct IFS
    contractions::Vector{Contraction}
    contractivity::Real
    function IFS(contractions, contractivity)
        if  0 < contractivity < 1
            new(contractions, contractivity)
        else
            throw("Contractivity must be real between 0 and 1")
        end
    end
end

function IFS(contractions)
    contractivities = [contractions[k].contractivity for k = 1 : length(contractions)]
    IFS(contractions, maximum(contractivities))
end

IFS(contractions::Array{Contraction, 2}) =
    IFS(reshape(contractions, length(contractions), ))

Base.show(io::IO, ifs::IFS) = print(io, "IFS\n
                                         \tcontractions: $(ifs.contractions)\n
                                         \tcontractivity: $(ifs.contractivity)")


# Define deterministic algorithm for iterated function system.
function deterministic_algorithm(ifs::IFS;
                                 initial::Union{Void, Vector{T} where {T<:Real}}=nothing,
                                 num_steps::Integer=10,
                                 num_track_points::Integer=2^10,
                                 num_update::Union{Void, Integer}=nothing,
                                 monitor::Union{Void, Function}=nothing)
    # Check number of update points.
    if num_update == nothing
        num_update = floor(num_steps / 5)
    end

    # Check initial set.
    domain_dim = ifs.contractions[1].domain.dim
    if initial == nothing
        initial = randn(domain_dim, 1)
    else
        if size(initial)[1] != domain_dim
            throw(DimensionMismatch("Inital does not match ifs domain dimension."))
        end
    end

    # If monitor is provided, construct a data channel
    if monitor != nothing
        channel = Channel(num_track_points)
    end

    # Compute the attractor
    num_procs = nprocs()
    num_contractions = length(ifs.contractions)
    set = initial
    for i = 1 : num_steps
        num_points = length(set)
        results = Matrix(size(set)[1], size(set)[2] * num_contractions)
        for j = 1 : num_contractions
            func = ifs.contractions[j].rule

            # Distribute the set to available processors
            chunk_size = floor(num_points / num_procs)
            processes = Vector{Future}(num_procs)
            for k = 1 : num_procs
                chunk = set[:, 1 + (k - 1) * chunk_size : k * chunk_size]
                processes[k] = @spawn mapslices(func, chunk, 1)
            end

            # Collect the results from the processors
            result = Array{Real, 2}(size(set))
            for k = 1 : num_procs
                result[:, 1 + (k - 1) * chunk_size : k * chunk_size] = fetch(processes[k])
            end
            results[:, 1 + (j - 1) * num_points : j * num_points] = result
            println("i: $i results: results")
            # Check the number of track points
            if size(results)[2] > num_track_points
                set = results[:, end - num_track_points : end]
            end

            # Moitor the points
            if monitor != nothing && size(set)[1] == 2 && mod(i, num_update)
                # message = "Iteration: $i Points: $(length(set))"
                # monitor(set, message)
                push!(channel, (set[1, :], set[2, :]))
                push!(channel, nothing)
                @schedule monitor(channel)
            end  # End of if
        end  # End of for
    end  # End of for

    return set

end  # End of deterministic_algorithm function

end # End of IteratedFunctionSystem
