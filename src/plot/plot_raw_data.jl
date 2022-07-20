using Plots
using DelimitedFiles

mutable struct Plot_parameters
    border::String


    function Plot_parameters()
        this = new()

        return this
    end
end


function plot_czechia(plt::Plots.Plot, border_data::String)
    border::String = border_data

    data = readdlm(border)
    # new plot
    plot!(plt, title = "Czechia", xlabel = "geodetic longitude", ylabel = "geodetic latitude")
    plot!(plt, xlims = (11.80, 19.20), xticks = 12:0.5:19.5)
    plot!(plt, ylims = (48.30, 51.30), yticks = 48.5:0.25:51.25)
    plot!(plt, data[:,3], data[:,2], label="Czechia boundaries")

end

function plot_grid(plt::Plots.Plot)
    plot!(plt)
end

function plot_points(plt::Plots.Plot, x::Vector{T}, y::Vector{T}) where (T<:AbstractFloat)
    plot!(plt, x, y)
end
