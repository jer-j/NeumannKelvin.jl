using Pkg
using Test

const VALIDATION_DIRECTORY = @__DIR__
const REPOSITORY_ROOT = normpath(joinpath(VALIDATION_DIRECTORY, "..", ".."))

Pkg.develop(path=REPOSITORY_ROOT)
Pkg.instantiate()

include(joinpath(VALIDATION_DIRECTORY, "NoblesseWigleyValidation.jl"))
using .NoblesseWigleyValidation

@testset "Noblesse Wigley geometry" begin
    quadrilaterals = NoblesseWigleyValidation.wigley_panels(8, 4; element=:quadrilateral)
    triangles = NoblesseWigleyValidation.wigley_panels(8, 4; element=:triangle)

    @test length(quadrilaterals) == 32
    @test length(triangles) == 64
    @test all(panel -> panel.x[3] < 0, quadrilaterals)
    @test all(panel -> panel.n[2] < 0, quadrilaterals)
    @test all(panel -> panel.n[2] < 0, triangles)
    @test NoblesseWigleyValidation.waterline_halfbreadth(0.0) == 0.05
    @test NoblesseWigleyValidation.waterline_halfbreadth(0.5) == 0.0
    @test NoblesseWigleyValidation.waterline_halfbreadth(0.6) == 0.0
end

@testset "Coarse NK smoke case" begin
    system = NoblesseWigleyValidation.solve_wigley_case(0.316; nx=6, nz=2, verbose=false)
    coefficients = NoblesseWigleyValidation.linear_hydrodynamic_coefficients(system)
    cut = NoblesseWigleyValidation.wave_cut(system; x=range(-0.6, 0.6; length=9))

    @test all(isfinite, values(coefficients))
    @test all(isfinite, cut.eta)
    @test cut.offset > 0
end
