
module IFS

# Imports
import Base.Sys.CPU_CORES
import StatsBase.sample
import StatsBase.Weights


# Exports
export deterministic_worker, deterministic_algorithm
export random_worker, random_algorithm


function slice(data::Array{<:Real, 2}, num_chunks::Int)
    chunk_size = round(Int, size(data)[2] / num_chunks)
    chunks = Vector{Array{<:Real, 2}}(num_chunks)
    for k = 1 : num_chunks - 1
        chunks[k] = data[:, 1 : chunk_size]
        data = data[:, chunk_size + 1 : end]
    end
    chunks[num_chunks] = data
end


function deterministic_worker(func::Function, chunk::AbstractArray{<:Real})
    return mapslices(func, chunk, 1)
end


function deterministic_algorithm(funcs::AbstractArray{Function},
                                 initial::AbstractArray{<:Real},
                                 num_iter::Integer;
                                 multi_procs::Bool=true)
    # Check the initial
    if ndims(initial) == 1
        initial = reshape(initial, length(initial), 1)
    end
    dim = size(initial)[1]

    # Compute the attractor
    num_funcs = length(funcs)
    num_procs = CPU_CORES
    set = initial
    for i = 1 : num_iter
        num_points = size(set)[2]
        next_set = zeros(dim, size(set)[2] * num_funcs)
        for j = 1 : num_funcs
            if multi_procs && num_procs > 2 && num_points >= num_procs
                # Distribute the processes
                chunks = slice(set, num_procs)
                processes = Vector{Future}(num_procs)
                for k = 1 : num_procs
                    processes[k] = @spawn deterministic_worker(funcs[j], chunks[j])
                end
                # Collect the results from the processes
                temp = zeros(size(set))
                for k = 1 : num_procs
                    temp[:, 1 + (k - 1) * chunk_size : k * chunk_size] = fetch(processes[k])
                end
            else
                temp = deterministic_worker(funcs[j], set)
            end
            next_set[:, 1 + (j - 1) * num_points : j * num_points] = temp
        end
        set = next_set
    end
    return set
end  # deterministic_algorithm


function random_worker(funcs::AbstractArray{Function},
                       probabilities::AbstractArray{<:Real},
                       initial::AbstractArray{<:Real},
                       num_iter::Integer;
                       discard::Bool=false)
    if ndims(initial) == 1
       initial = reshape(initial, length(initial), 1)
    end

    data = zeros(size(initial)[1], num_iter)
    data[:, 1] = initial
    for i = 1 : num_iter - 1
        func = sample(funcs, Weights(probabilities))    # Choose a random function `func`
        data[:, i + 1] = func(data[:, i])               # Apply `func` to the data.
    end
    if discard
        data =  data[:, end : end]
    end
    return data
end  # random_worker


function random_algorithm(funcs::AbstractArray{Function},
                          probabilities::AbstractArray{<:Real},
                          initial::AbstractArray{<:Real},
                          num_iter::Integer;
                          num_transients::Integer=100,
                          multi_procs::Bool=true)
    # Check the initial
    if ndims(initial) == 1
     initial = reshape(initial, length(initial), 1)
    end

    # Transient iterations
    initial = random_worker(funcs, probabilities, initial, num_transients, discard=true)

    # Compute the attractor
    num_procs = CPU_CORES
    if multi_procs && num_procs > 2 && num_iter >= num_procs  # Use multiple
        iters = slice(ones(Int, 1, num_iter), num_procs)
        # Launch the remote processes.
        processes = Vector{Future}(num_procs)
        for i = 1 : num_procs
            processes[i] = @spawn random_worker(funcs, probabilities, initial, iters[i], discard=false)
        end
        # Collect the results from the remote processes.
        set = zeros(size(initial)[1], 1)
        for i = 1 : num_procs
            value = fetch(processes[i])
            set = [set value]
        end
        set = set[:, 2:end]
    else
        set = random_worker(funcs, probabilities, initial, num_iter, discard=false)
    end  # End of if-else

    return set

end  # random_algorithm

end  # module
