using Pkg

const VALIDATION_DIRECTORY = @__DIR__
const REPOSITORY_ROOT = normpath(joinpath(VALIDATION_DIRECTORY, "..", ".."))

Pkg.develop(path=REPOSITORY_ROOT)
Pkg.instantiate()

include(joinpath(VALIDATION_DIRECTORY, "NoblesseWigleyValidation.jl"))
using .NoblesseWigleyValidation

parse_bool(value) = lowercase(value) in ("1", "true", "yes", "on")

nx = parse(Int, get(ENV, "NK_NX", "40"))
nz = parse(Int, get(ENV, "NK_NZ", "8"))
element = Symbol(get(ENV, "NK_ELEMENT", "quadrilateral"))
filter_waves = parse_bool(get(ENV, "NK_FILTER", "true"))
contour = parse_bool(get(ENV, "NK_CONTOUR", "true"))
solver = Symbol(get(ENV, "NK_SOLVER", "direct"))
output_directory = get(ENV, "NK_OUTPUT", joinpath(VALIDATION_DIRECTORY, "results"))

froude_numbers = if haskey(ENV, "NK_FROUDE_NUMBERS")
    parse.(Float64, split(ENV["NK_FROUDE_NUMBERS"], ','))
else
    NOBLESSE_FROUDE_NUMBERS
end

run_validation(
    output_directory;
    froude_numbers,
    nx,
    nz,
    element,
    filter=filter_waves,
    contour,
    solver,
)

if parse_bool(get(ENV, "NK_PLOT", "true"))
    include(joinpath(VALIDATION_DIRECTORY, "plot.jl"))
    plot_validation(output_directory)
end
