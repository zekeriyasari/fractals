# Test script to test remote plot

include("../src/monitors.jl")

using Monitors

function test_update_plot(initial::Vector{<:Real};
           num_track_points::Integer=2^10,
           monitor::Union{Void, Function}=nothing)
    # Construct a data channel
    channel = Channel(num_track_points)
    @schedule monitor(channel)  # Launch the remote processes

    # Construct data and push into the channel
    x = linspace(0, 6 * pi, 100)
    for phase in linspace(0, 10 * pi, 100)
        x_data = x
        y_data = sin.(x_data + phase)
        push!(channel, (x_data, y_data))  # Push data to channel
    end
    push!(channel, nothing)  # Poisson pill for the channel
end

test_update_plot([i for i = 0 : 0.1 : 6 * pi], monitor=update_plot)
