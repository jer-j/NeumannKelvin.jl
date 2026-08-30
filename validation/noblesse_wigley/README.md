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

The existing `NKPanelSystem` is an indirect Kelvin-source formulation. The `contour=true`
option is the package's Baar-type source-panel correction and is not the complete direct NK
waterline integral. Both `NK_CONTOUR=true/false` and `NK_FILTER=true/false` should therefore
be included in the final sensitivity study.
