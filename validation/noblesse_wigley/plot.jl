using CairoMakie
using CSV

"""
    plot_validation(output_directory)

Create paper-style wave-cut and hydrodynamic-coefficient figures from validation CSV files.
"""
function plot_validation(output_directory)
    coefficients = CSV.File(joinpath(output_directory, "coefficients.csv")) |> collect

    wave_figure = Figure(size=(1100, 900))
    for (index, row) in enumerate(coefficients)
        axis = Axis(
            wave_figure[(index - 1) ÷ 2 + 1, (index - 1) % 2 + 1],
            title="F = $(row.Fn)",
            xlabel=index > 4 ? "x/L" : "",
            ylabel=isodd(index) ? "η = e/F²" : "",
            limits=((-0.6, 0.6), (-0.12, 0.25)),
        )
        filename = "wavecut_Fn_$(replace(string(row.Fn), '.' => 'p')).csv"
        wave_data = CSV.File(joinpath(output_directory, filename)) |> collect
        lines!(axis, getproperty.(wave_data, :x), getproperty.(wave_data, :eta); linewidth=2)
        vlines!(axis, [-0.5, 0.5]; color=(:gray, 0.45), linestyle=:dot)
    end
    save(joinpath(output_directory, "noblesse_wigley_wavecuts.png"), wave_figure)

    coefficient_figure = Figure(size=(1000, 750))
    Fn = getproperty.(coefficients, :Fn)
    quantities = (
        (:c_wave_near, "Cwave,near"),
        (:sinkage, "h/L"),
        (:trim, "trim [rad]"),
    )
    for (index, (field, label)) in enumerate(quantities)
        axis = Axis(coefficient_figure[index, 1]; xlabel="F", ylabel=label)
        scatterlines!(axis, Fn, getproperty.(coefficients, field); markersize=9)
    end
    save(joinpath(output_directory, "noblesse_wigley_coefficients.png"), coefficient_figure)
    nothing
end
