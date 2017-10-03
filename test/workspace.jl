using PyPlot

function update_plot(channel)
    fig = figure()
    ax = fig[:add_subplot](111)
    x = rand()
    y = rand()
    curve, = ax[:plot]([], [], ".")
    while true
        data = take!(channel)
        if data == nothing
            break
        end
        x_data, y_data = data
        curve[:set_xdata](x_data)
        curve[:set_ydata](y_data)
        fig[:canvas][:draw]()
        fig[:canvas][:flush_events]()
        ax[:autoscale_view]()
        ax[:relim]()
        # sleep(0.1)
    end
end

function f()
    channel = Channel(500)
    x = linspace(0, 6 * pi, 100)
    y = sin.(x)
    for phase in linspace(0, 10 * pi, 100)
        x_data = x
        y_data = sin.(x_data + phase)
        push!(channel, (x_data, y_data))
    end
    push!(channel, nothing)
    @schedule update_plot(channel)
end

f()
