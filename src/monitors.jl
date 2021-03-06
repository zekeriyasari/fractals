"""
Monitors module for data visualization
"""
module Monitors

# TODO: Three dimensional monitoring plots are to be included to the module.

using PyPlot

export update_plot


"""
    update_plot(channel::Channel{Any})

Updates the data of the plot.

`update_plot` is basically a scope to visualise two dimensional data
Plotting is performed in a remote process.
"""
function update_plot(channel::Channel{Any};
                     timeout::Real=0.1,
                     sleep_time::Real=0.)
    # Construct the plot window
    fig = figure()
    ax = fig[:add_subplot](111)
    # ax[:grid]()
    ax[:set_autoscaley_on](true)

    # Dummy plot for initialization
    curve, = ax[:plot]([], [], ".")

    # Update the plot curve by listening the channel
    while true
        # Take the data in the channel if available
        if isready(channel)
            data = take!(channel)
        else
            try
                wait(timeout)
            catch
                error("No data available in the channel.")
            end
        end

        # Check the poission pill to terminate listening the channel.
        if data == nothing
            break
        end
        # data = take!(channel)
        # Update the curve of the plot with the data taken from the channel.
        if length(data) == 2
            x_data, y_data = data
        elseif length(data) == 3
            x_data, y_data, message = data
        end

        curve[:set_xdata](x_data)
        curve[:set_ydata](y_data)
        fig[:canvas][:draw]()
        fig[:canvas][:flush_events]()
        ax[:set_xlim](minimum(x_data), maximum(x_data))
        ax[:set_ylim](minimum(y_data), maximum(y_data))
        ax[:autoscale_view]()
        ax[:relim]()
        if length(data) == 3
            ax[:set_title](message)
        end

        # If provided, freeze the plot
        if sleep_time != 0
            sleep(sleep_time)
        end
    end  # End of while
end  # End of update_plot function

end  # End of the module
