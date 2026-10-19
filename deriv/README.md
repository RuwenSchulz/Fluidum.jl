# deriv/: symbolic derivations of the generated matrices

Each `*_generated.jl` file in `src/Matrix/` is written by one Wolfram Language script in this
directory. The script derives the quasi-linear system covariantly from its tensor equations,
runs its own consistency checks, writes the Julia file, and exits with status 1 if any check failed.
`derive_two_way.wls` is stricter: it writes into `src/Matrix/` only when every check passed, and
otherwise into its scratch directory. Do not edit a generated file by hand: change the script and
regenerate.

| script | writes | system |
|---|---|---|
| `derive_2p1d_viscous.wls` | `viscous_generated.jl` | Israel–Stewart viscous hydro plus the DNMR term δ_ππ, 1+1D cylindrical and 2+1D Cartesian Milne |
| `derive_hq_2p1d.wls` | `HQ_2p1d_BG_generated.jl` | charm first moment (α, ν^x, ν^y) on a prescribed 2+1D background, Fick closure |
| `derive_hq_consistent_2p1d.wls` | `HQ_2p1d_BG_consistent_generated.jl` | the same with the consistent (full ∇P) first moment |
| `derive_two_way.wls` (+ `_gates`, `_export`) | `two_way_generated.jl` | ideal and viscous medium, charm first and second moment, and the two-way coupled system (charm in the total T^{μν}); 0+1D, 1+1D cylindrical, 2+1D Cartesian; EOS kept symbolic |

Helpers: `julia_source_reader.wl` reads a matrix out of a `.jl` file in `src/Matrix/`, so the
scripts can compare their result with the code itself rather than with a hand transcription.
`refpoints_csv.wl` writes reference states (the derivation evaluated at random rational points) as CSV.

## How each derivation is checked

Each script compares its own result against an independent known answer before it writes anything:

- the cylindrical (1+1D) reduction must reproduce the existing hand-written 1-D matrix entry by
  entry (`one_d_viscous_matrix`, `HQ_const_BG.jl`, the consistent 1-D rows), read from the `.jl`
  file;
- metric compatibility, transversality and tracelessness of every constructed tensor hold
  identically;
- the system is quasi-linear (no second derivatives survive);
- the Bjorken and slab limits reduce to the known 0+1D and decoupled equations;
- for the two-way system, the fluid rows equal ∇_μ(T_m + T_Q)^{μν} computed by an independent
  route, and there is no momentum exchange at rest.

The Julia tests in `test/` and `bench/` evaluate the generated functions against these
derivations numerically.

## Running

Requires Mathematica / `wolframscript` (no extra packages). From the package root:

```sh
wolframscript -file deriv/derive_2p1d_viscous.wls          # seconds
wolframscript -file deriv/derive_hq_2p1d.wls
wolframscript -file deriv/derive_hq_consistent_2p1d.wls
wolframscript -file deriv/derive_two_way.wls               # ~5 min; --charts=bj,cyl,car
```

A regeneration with unchanged scripts is byte-identical to the committed files.

Environment variables (all optional):

| variable | default | meaning |
|---|---|---|
| `FLUIDUM_REFPOINTS_DIR` | `$TemporaryDirectory/fluidum_refpoints` | where the reference-state CSVs go |
| `FLUIDUM_DERIV_PROJECT` | the package root | Julia project for the Julia-side checks of `derive_two_way.wls` (needs StaticArrays and SpecialFunctions) |
| `FLUIDUM_EXTERNAL_REFPOINTS` | unset | directory with the 1+1D reference CSVs of the first- and second-moment derivations; when unset, the three checks that use them are skipped |
| `TWOWAY_OUT`, `TWOWAY_JL`, `TWOWAY_REFDIR`, `FLUIDUM_MATRIX_DIR` | see the header of `derive_two_way.wls` | scratch and output locations of the two-way derivation |
