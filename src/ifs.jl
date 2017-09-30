
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
                                 max_iter::Integer=5)
    # Check initial set.
    domain_dim = ifs.contractions[1].domain.dim  # Contraction domains assumed to be same.
    if initial == nothing
        initial = randn(domain_dim, 1)
    else
        if size(initial)[1] != domain_dim
            throw(DimensionMismatch("Inital does not match ifs domain dimension."))
        end
    end

    # Compute the attractor
    # TODO: Use Hausforff distance as a stopping criteria.
    num_contractions = length(ifs.contractions)
    x = initial
    for i = 1 : max_iter
        processes = [@spawn mapslices(ifs.contractions[k].rule, x, 1) for k = 1 : num_contractions]
        x = reduce(hcat, fetch(processes[k]) for k = 1 : num_contractions)
    end

    return x

end  # End of deterministic_algorithm function

end # End of IteratedFunctionSystem
