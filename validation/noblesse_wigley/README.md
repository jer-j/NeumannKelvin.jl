# Noblesse Wigley-hull NK validation

This validation applies the existing `NKPanelSystem` to the Wigley hull and Froude-number
cases used in Section 4 of Noblesse et al. (2013), *The NM theory of ship waves*. It is an
NK baseline for later comparisons with NM and NN implementations. It does not claim to
reproduce the NM method.

Reference: F. Noblesse, F. Huang, and C. Yang, "The Neumann-Michell theory of ship waves,"
*Journal of Engineering Mathematics*, vol. 79, pp. 51-71, 2013.

The nondimensional hull is

```math
y = \pm \frac{b}{2}(1-4x^2)\left(1-\frac{z^2}{d^2}\right),
\qquad b=0.1,\quad d=\frac{1}{16}.
```

The validation cases are

```math
F\in\{0.250,0.267,0.289,0.316,0.354,0.408\}.
```

The output quantities follow the paper's normalization:

```math
\eta=\frac{e}{F^2}=\phi_x,
```

```math
C_{\mathrm{wave}}^{\mathrm{near}}
=\int_{\Sigma_H}n_x\phi_x\,dS,
```

```math
h=F^2\frac{C_{\mathrm{lift}}}{a_0^W},
\qquad
\tau=F^2\frac{C_{\mathrm{pitch}}}{a_2^W}.
```

For this Wigley waterplane,

```math
a_0^W=\frac{2b}{3},
\qquad
a_2^W=\frac{b}{30}.
```

## Run

From the repository root:

```bash
julia --project=validation/noblesse_wigley validation/noblesse_wigley/run.jl
```

The default run uses a moderate `40 x 8` quadrilateral discretization of one half hull.
The following environment variables control the calculation:

| Variable | Default | Meaning |
|---|---:|---|
| `NK_NX` | `40` | Longitudinal cells |
| `NK_NZ` | `8` | Vertical cells |
| `NK_ELEMENT` | `quadrilateral` | `quadrilateral` or `triangle` |
| `NK_FILTER` | `true` | Apply the package's panel-scale wave filter |
| `NK_CONTOUR` | `true` | Apply the Baar waterline contour correction |
| `NK_SOLVER` | `direct` | `direct` or `gmres` |
| `NK_PLOT` | `true` | Create PNG figures with CairoMakie |
| `NK_FROUDE_NUMBERS` | all six cases | Comma-separated subset |
| `NK_OUTPUT` | `validation/noblesse_wigley/results` | Output directory |

For example:

```bash
NK_NX=80 NK_NZ=12 NK_FILTER=false NK_CONTOUR=true \
    julia --project=validation/noblesse_wigley \
    validation/noblesse_wigley/run.jl
```

The nominal paper spacing is obtained with `NK_NX=200` and `NK_NZ=20`. This is an
expensive Kelvin-panel calculation because the present implementation has no fast summation
for the Kelvin Green function. A mesh sequence should be completed before interpreting the
paper-resolution result.

## Interpretation

The paper publishes the experimental and NM results graphically, not as numerical tables.
The scripts therefore produce reproducible NK curves in the same coordinates and
normalization, but they do not embed estimated values digitized from the figures. Quantitative
error norms should only be added when an authoritative numerical data set or documented
digitization is available.

The calculated wave cut is evaluated slightly outside the hull waterline. The offset is stored
in every wave-cut CSV file. This avoids evaluating the source-panel representation precisely
on a waterline edge and makes the numerical convention reproducible.

`coefficients.csv` also reports two numerical diagnostics. `c_wave_bernoulli` uses the complete
perturbation-pressure expression instead of the linear pressure used for `c_wave_near`.
`bc_max_abs` and `bc_rms` measure the no-penetration residual at the collocation points.

## Reproducibility checks

The checked reference data contain a three-level mesh sequence and a complete two-by-two
sensitivity study for `filter` and `contour` at `80 x 12`. They are stored in
[`reference/mesh_convergence.csv`](reference/mesh_convergence.csv) and
[`reference/model_sensitivity_80x12.csv`](reference/model_sensitivity_80x12.csv).

For the `80 x 12`, filtered, contour-corrected run, the main coefficients are:

| `F` | `C_wave^near` | `h` | `tau` |
|---:|---:|---:|---:|
| 0.250 | 0.00008287 | 0.0011766 | 0.0003736 |
| 0.267 | 0.00007105 | 0.0013720 | 0.0006879 |
| 0.289 | 0.00012157 | 0.0016647 | 0.0013992 |
| 0.316 | 0.00010735 | 0.0019334 | 0.0013763 |
| 0.354 | 0.00011918 | 0.0026010 | 0.0032487 |
| 0.408 | 0.00015724 | 0.0037684 | 0.0097756 |

The direct solves satisfy the discrete body boundary condition to about `7e-17` RMS. Between
the `40 x 8` and `80 x 12` meshes, sinkage changes by 2.7-4.0 percent, while the drag and trim
remain more oscillatory and can change by roughly 12 percent at individual Froude numbers.
The wave filter has a modest effect on this mesh. The Baar contour correction has a much larger
effect on drag and trim, especially at the higher Froude numbers.

These calculations reproduce the Noblesse Wigley geometry, Froude-number cases, and output
normalizations, but they do not reproduce the paper's NM curves. In particular, the present
NK near-field drag lies visibly below the NM and experimental curves in Fig. 3. This is a
formulation-level result, not a linear-solver residual: the package uses an indirect
Kelvin-source equation and an optional panel correction, whereas the paper solves the NM
equation. A quantitative experimental error norm still requires a documented digitization or
the original tabulated measurements.

The existing `NKPanelSystem` is an indirect Kelvin-source formulation. The `contour=true`
option is the package's Baar-type source-panel correction and is not the complete direct NK
waterline integral. Both `NK_CONTOUR=true/false` and `NK_FILTER=true/false` should therefore
be included in the final sensitivity study.
