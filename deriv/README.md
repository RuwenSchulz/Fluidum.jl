# deriv/: symbolic derivations of the generated matrices

Some matrix functions in `src/Matrix/` are written by a Wolfram Language script in this
directory. They sit in the matrix file of their system, between the lines
`# symbolic derivation in: deriv/<script>` and `# end of symbolic derivation`; the rest of the file is
hand-written. The script derives the quasi-linear system covariantly from its tensor equations,
runs its own consistency checks, and rewrites the block only if every check passed (exit status 1
otherwise). Do not edit a generated block by hand: change the script and regenerate.

| script | writes | system |
|---|---|---|
| `derive_2p1d_viscous.wls` | `two_d_viscous_matrix_derived` in `2d_viscous.jl` | Israel–Stewart viscous hydro plus the DNMR term δ_ππ, 2+1D Cartesian Milne (and 1+1D cylindrical as the control) |
| `derive_hq_2p1d.wls` | `two_d_HQ_BG_matrix`, `two_d_HQ_live_rows` in `2d_viscous_HQ_BG.jl` | charm first moment (α, ν^x, ν^y) on a prescribed 2+1D background, Fick closure |
| `derive_two_way.wls` (+ `_gates`, `_export`) | `one_d_HQ_BG_consistent_matrix`, `one_d_HQ_BG_m2_matrix`, `two_d_HQ_BG_consistent_matrix`, `two_d_HQ_BG_m2_matrix`, `one_d_two_way_matrix(_conformal)`, `two_d_two_way_matrix(_conformal)` | charm first moment (consistent closure) and second moment on a prescribed background, and the two-way coupled system (charm in the total T^{μν}); 1+1D cylindrical and 2+1D Cartesian; EOS kept symbolic |

Helper: `julia_source_reader.wl` reads a matrix out of a `.jl` file in `src/Matrix/`, so the
scripts can compare their result with the code itself rather than with a hand transcription, and
writes the generated block back.

## How each derivation is checked

Each script compares its own result against an independent known answer before it writes anything:

- the cylindrical (1+1D) reduction must reproduce the existing hand-written 1-D matrix entry by
  entry (`one_d_viscous_matrix`, `one_d_viscous_matrix_fugacity_BG_only`), read from the `.jl`
  file;
- metric compatibility, transversality and tracelessness of every constructed tensor hold
  identically;
- the system is quasi-linear (no second derivatives survive);
- the Bjorken and slab limits reduce to the known 0+1D and decoupled equations, and the 2+1D
  systems reduce to the 1+1D ones on axisymmetric flow;
- for the two-way system, the fluid rows equal ∇_μ(T_m + T_Q)^{μν} computed by an independent
  route, the second moment is passive, and there is no momentum exchange at rest.

The Julia tests in `test/` check the generated functions against closed forms.

## Running

Requires Mathematica / `wolframscript` (no extra packages). From the package root:

```sh
wolframscript -file deriv/derive_2p1d_viscous.wls          # seconds
wolframscript -file deriv/derive_hq_2p1d.wls
wolframscript -file deriv/derive_two_way.wls               # ~5 min; --charts=bj,cyl,car
```

A regeneration with unchanged scripts is byte-identical to the committed files.
`FLUIDUM_DERIV_RENDER=1` also writes the rendered derivation into `rendered/` (see below).
`derive_two_way.wls` uses a scratch directory, `TWOWAY_OUT` (default
`$TemporaryDirectory/derive_two_way_out`), and runs Julia for two of its checks in the project
`FLUIDUM_DERIV_PROJECT` (default the package root).

## Rendered output

`rendered/` holds the derivations in readable form, written by the scripts themselves from the
matrices, sources and check results they computed:

- `derivations.pdf`: one section per script (setup, every nonzero matrix and source entry, the
  table of checks). `derivations.tex` is the master; `viscous.tex`, `hq_2p1d.tex` and
  `two_way.tex` are the generated fragments.
- `<name>.nb`: the same content as notebooks, with the script source at the end; written locally,
  not tracked.

For the two-way script, the 1+1D charm systems and the 1+1D conformal two-way system are written
out in full; the generic-EOS 1+1D two-way system and all 2+1D systems are shown as field layout
and sparsity pattern (their entries are in the generated blocks).

The rendering code is `render.wl`. It is loaded only when `FLUIDUM_DERIV_RENDER=1`; without it
the scripts behave as before. Every displayed entry is compared with the derived expression at
random exact points before it is written. To regenerate everything:

```sh
deriv/render_all.sh
```

or one script at a time, then the PDF:

```sh
FLUIDUM_DERIV_RENDER=1 wolframscript -file deriv/<script>.wls
cd deriv/rendered && latexmk -pdf derivations.tex
```
