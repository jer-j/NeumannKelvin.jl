module NoblesseWigleyValidation

using CSV
using GeometryBasics
using NeumannKelvin
using StaticArrays

export NOBLESSE_FROUDE_NUMBERS,
    linear_hydrodynamic_coefficients,
    run_validation,
    solve_wigley_case,
    waterline_halfbreadth,
    wave_cut,
    wigley_panels,
    wigley_point

const NOBLESSE_FROUDE_NUMBERS = (0.250, 0.267, 0.289, 0.316, 0.354, 0.408)

"""
    wigley_point(x, z; beam=0.1, draft=1 / 16)

Return a point on the negative-`y` half of the Wigley hull used by Noblesse et al.
The length is one, the bow is at `x=0.5`, and the undisturbed free surface is `z=0`.
"""
function wigley_point(x, z; beam=0.1, draft=1 / 16)
    shape_x = max(zero(x), one(x) - 4x^2)
    shape_z = max(zero(z), one(z) - (z / draft)^2)
    SA[x, -beam * shape_x * shape_z / 2, z]
end

"""
    wigley_panels(nx=40, nz=8; element=:quadrilateral, beam=0.1, draft=1 / 16)

Create a uniform panelization of the negative-`y` half of the Noblesse Wigley hull.
`element=:quadrilateral` uses one quadrature panel per `(x,z)` cell.
`element=:triangle` splits every cell into two flat triangles.

The paper uses `nx=200`, `nz=20`, `dx=0.005`, and `dz=0.003125` with flat
triangles. Lower resolutions are intended for convergence studies and CI smoke tests.
"""
function wigley_panels(
    nx::Integer=40,
    nz::Integer=8;
    element::Symbol=:quadrilateral,
    beam=0.1,
    draft=1 / 16,
)
    nx > 0 || throw(ArgumentError("nx must be positive"))
    nz > 0 || throw(ArgumentError("nz must be positive"))
    element in (:quadrilateral, :triangle) ||
        throw(ArgumentError("element must be :quadrilateral or :triangle"))

    T = promote_type(typeof(float(beam)), typeof(float(draft)))
    beam_T, draft_T = T(beam), T(draft)
    dx, dz = one(T) / nx, draft_T / nz
    surface(x, z) = wigley_point(x, z; beam=beam_T, draft=draft_T)

    if element === :quadrilateral
        xc = range(-T(0.5) + dx / 2; step=dx, length=nx)
        zc = range(-draft_T + dz / 2; step=dz, length=nz)
        return Table(measure.(surface, xc, zc', dx, dz))
    end

    first_panel = measure(
        surface(-T(0.5), -draft_T),
        surface(-T(0.5) + dx, -draft_T),
        surface(-T(0.5) + dx, -draft_T + dz),
    )
    panels = typeof(first_panel)[]
    sizehint!(panels, 2nx * nz)
    for iz in 1:nz, ix in 1:nx
        x0, x1 = -T(0.5) + (ix - 1) * dx, -T(0.5) + ix * dx
        z0, z1 = -draft_T + (iz - 1) * dz, -draft_T + iz * dz
        p00, p10 = surface(x0, z0), surface(x1, z0)
        p01, p11 = surface(x0, z1), surface(x1, z1)
        push!(panels, measure(p00, p10, p11))
        push!(panels, measure(p00, p11, p01))
    end
    Table(panels)
end

"""
    waterline_halfbreadth(x; beam=0.1)

Return the Wigley half-breadth at `z=0`. Points beyond the bow and stern return zero.
"""
waterline_halfbreadth(x; beam=0.1) = beam * max(zero(x), one(x) - 4x^2) / 2

"""
    solve_wigley_case(Fn; kwargs...)

Solve one steady Neumann-Kelvin Wigley-hull case and return the solved system.
The half hull is reflected across `y=0` by `sym_axes=2`.
"""
function solve_wigley_case(
    Fn;
    nx::Integer=40,
    nz::Integer=8,
    element::Symbol=:quadrilateral,
    filter::Bool=true,
    contour::Bool=true,
    solver::Symbol=:direct,
    atol=1e-7,
    verbose::Bool=false,
)
    panels = wigley_panels(nx, nz; element)
    system = NKPanelSystem(panels; ℓ=Fn^2, sym_axes=2, filter, contour)
    if solver === :direct
        return directsolve!(system; verbose)
    elseif solver === :gmres
        return gmressolve!(system; atol, verbose)
    end
    throw(ArgumentError("solver must be :direct or :gmres"))
end

"""
    linear_hydrodynamic_coefficients(system; beam=0.1)

Calculate Noblesse's linear near-field coefficients from `phi_x` on the half hull.
The returned sinkage and trim are the equilibrium values obtained from linear hydrostatic
restoring coefficients for the fixed Wigley waterplane.
"""
function linear_hydrodynamic_coefficients(system; beam=0.1)
    panels = system.body
    velocity = u(system)
    phi_x = getindex.(velocity, 1) .- system.U[1]

    c_wave_near = 2 * sum(panels.n[i][1] * phi_x[i] * panels.dA[i] for i in eachindex(phi_x))
    c_lift = 2 * sum(panels.n[i][3] * phi_x[i] * panels.dA[i] for i in eachindex(phi_x))
    c_pitch = 2 * sum(
        (panels.n[i][1] * panels.x[i][3] - panels.n[i][3] * panels.x[i][1]) *
        phi_x[i] * panels.dA[i] for i in eachindex(phi_x)
    )

    waterplane_area = 2 * beam / 3
    waterplane_second_moment = beam / 30
    Fn2 = system.args.ℓ
    sinkage = Fn2 * c_lift / waterplane_area
    trim = Fn2 * c_pitch / waterplane_second_moment

    (; c_wave_near, c_lift, c_pitch, sinkage, trim)
end

"""
    wave_cut(system; x=range(-0.6, 0.6, 241), beam=0.1, offset=nothing)

Evaluate `eta=e/F^2=phi_x` just outside the hull waterline. A small transverse offset keeps
the field points outside the source-panel boundary and makes the numerical convention explicit.
"""
function wave_cut(
    system;
    x=range(-0.6, 0.6; length=241),
    beam=0.1,
    offset=nothing,
)
    panel_scale = sqrt(sum(system.body.dA) / length(system.body))
    transverse_offset = isnothing(offset) ? panel_scale / 20 : offset
    eta = map(x) do x_coordinate
        y = -waterline_halfbreadth(x_coordinate; beam) - transverse_offset
        ζ(SA[x_coordinate, y, zero(x_coordinate)], system)
    end
    (; x=collect(x), eta, offset=transverse_offset)
end

"""
    run_validation(output_directory; kwargs...)

Run the six Froude-number cases reported by Noblesse et al. and write coefficient and
wave-cut CSV files. Results are overwritten within `output_directory`.
"""
function run_validation(
    output_directory;
    froude_numbers=NOBLESSE_FROUDE_NUMBERS,
    nx::Integer=40,
    nz::Integer=8,
    element::Symbol=:quadrilateral,
    filter::Bool=true,
    contour::Bool=true,
    solver::Symbol=:direct,
    atol=1e-7,
    verbose::Bool=true,
)
    mkpath(output_directory)
    coefficients = NamedTuple[]

    for Fn in froude_numbers
        verbose && println("Solving Noblesse Wigley case Fn=$(Fn), nx=$(nx), nz=$(nz)")
        elapsed = @elapsed system = solve_wigley_case(
            Fn; nx, nz, element, filter, contour, solver, atol, verbose=false
        )
        loads = linear_hydrodynamic_coefficients(system)
        cut = wave_cut(system)

        push!(coefficients, (
            Fn,
            nx,
            nz,
            element=String(element),
            filter,
            contour,
            solver=String(solver),
            panels=length(system.body),
            elapsed_seconds=elapsed,
            loads...,
        ))
        CSV.write(
            joinpath(output_directory, "wavecut_Fn_$(replace(string(Fn), '.' => 'p')).csv"),
            (; x=cut.x, eta=cut.eta, offset=fill(cut.offset, length(cut.x))),
        )
    end

    coefficient_path = joinpath(output_directory, "coefficients.csv")
    CSV.write(coefficient_path, coefficients)
    coefficient_path
end

end
