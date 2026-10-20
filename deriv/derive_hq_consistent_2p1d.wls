#!/usr/bin/env wolframscript
(* ============================================================================
   derive_hq_consistent_2p1d.wls -- the consistent charm first moment in 2+1D.

   Companion of derive_hq_2p1d.wls, which derives the grad-alpha-only row.  This
   file derives the consistent (full grad-P) row in Cartesian Milne; the 1-D
   cylindrical consistent row is `hqc_matrices` in HQ_const_BG_consistent.jl.
   The grad-alpha-only row has no temperature-gradient channel; on a realistic
   Pb+Pb cell that channel is larger than the alpha drive.

   Physics.  Charm is a passive tracer with T_Q^{mu nu} the charm stress
   tensor.  The consistent first moment is the Delta-projected divergence

       (D_s/T) Delta^mu_lambda Nabla_nu T_Q^{nu lambda} + nu^mu = 0,
       T_Q^{ab} = (n h - P0) u^a u^b + P0 Delta^{ab} + h (u^a nu^b + nu^a u^b),
       P0 = n T,   tau_n = D_s h / T,   kappa = D_s dn/dalpha.

   Its Navier-Stokes limit at rest is

       nu^i = -kappa d_i alpha - (D_s/T)(n + T dn/dT) d_i T

   The first term is Fick's law; the second is the full-grad-P channel, with
   `n + T dn/dT = n h / T`.

   Checks (printed with the prefix HC-)
     G0  metric compatibility, both charts
     G1  divergence check: Nabla.V for V=(X(tau),0,0,0) must be X' + X/tau
     G2  u.u == -1 and u.nu == 0 identically, both charts
     G3  the consistent system is quasi-linear
     C1  cylindrical reduction: the static, at-rest, nu = 0 source is (Ds/T) d_r P0,
         the source of the 1-D consistent row `hqc_matrices`.  This is the control
         for the Cartesian chart.
     C2  the at-rest NS limit is  -kappa d_i alpha - (Ds/T)(n + T dn_dT) d_i T
     C3  consistent - grad-alpha-only == the named extra terms, identically (x-row):
         (Ds/T)(n + T dn/dT) grad_perp T,  tau_n n a^mu,  tau_n Delta (nu.Nabla) u,
         nu^mu (tau_n theta + (Ds/T) D h),  with tau_n = Ds h/T, kappa = Ds dn/dalpha
     C4  the slab limit decouples nu^y
     C5  setting the extra terms to zero recovers the grad-alpha-only row

   Run:    wolframscript -file deriv/derive_hq_consistent_2p1d.wls
   Writes: src/Matrix/HQ_2p1d_BG_consistent_generated.jl, and
           refpoints_hq2d_consistent.csv (the derivation at 12 random states) to
           $FLUIDUM_REFPOINTS_DIR (default a temporary directory)
   ============================================================================ *)

$HistoryLength = 0;
pass = 0; fail = 0; gateLog = {};
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  AppendTo[gateLog, {"HC-" <> id, TrueQ[ok], msg}];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], "HC-" <> id, "  ", msg]);
PKG = ParentDirectory[DirectoryName[$InputFileName]];   (* the Fluidum.jl package root *)
DERIV = DirectoryName[$InputFileName];
REFDIR = With[{v = Environment["FLUIDUM_REFPOINTS_DIR"]},
  If[StringQ[v] && v =!= "", v, FileNameJoin[{$TemporaryDirectory, "fluidum_refpoints"}]]];
Get[FileNameJoin[{DERIV, "refpoints_csv.wl"}]];

(* ------------------------------------------------------------------ geometry *)
mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam|>];

divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] +
  Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];

DcovV[geo_, u_, W_] := Table[Sum[u[[l]] D[W[[m]], geo["x"][[l]]], {l, 4}] +
  Sum[u[[l]] geo["Gam"][[m, l, s]] W[[s]], {l, 4}, {s, 4}], {m, 4}];

(* Nabla_nu T^{nu lambda} for a symmetric contravariant rank-2 tensor, with both
   Christoffel contractions written out.  A basis-indexed covariant derivative drops
   the Christoffel/sqrt(g) part and loses the geometric dilution
   tau_n nu (u^tau/tau + u^r/r); G1 and C1 would catch that. *)
divT2[geo_, T_] := Table[
  Sum[D[T[[a, b]], geo["x"][[a]]], {a, 4}] +
  Sum[geo["Gam"][[a, a, s]] T[[s, b]], {a, 4}, {s, 4}] +
  Sum[geo["Gam"][[b, a, s]] T[[a, s]], {a, 4}, {s, 4}], {b, 4}];

metricCompat[geo_] := Simplify @ Table[
  D[geo["g"][[m, n2]], geo["x"][[l]]] - Sum[geo["Gam"][[s, l, m]] geo["g"][[s, n2]], {s, 4}]
  - Sum[geo["Gam"][[s, l, n2]] geo["g"][[m, s]], {s, 4}], {l, 4}, {m, 4}, {n2, 4}];

CYL = mkGeom[{tau, r, phi, eta}, {-1, 1, r^2, tau^2}];
CAR = mkGeom[{tau, xx, yy, eta}, {-1, 1, 1, tau^2}];

Print["== G0/G1  geometry =="];
gate["G0cyl", metricCompat[CYL] === ConstantArray[0, {4, 4, 4}], "cylindrical Milne: Nabla g == 0"];
gate["G0car", metricCompat[CAR] === ConstantArray[0, {4, 4, 4}], "Cartesian   Milne: Nabla g == 0"];
gate["G1cyl", Simplify[divV[CYL, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "cylindrical: Nabla.V = X' + X/tau"];
gate["G1car", Simplify[divV[CAR, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "Cartesian:   Nabla.V = X' + X/tau"];

(* --------------------------------------------------------------- equations *)
(* Delta^{ab} = g^{ab} + u^a u^b ; Delta^a_b = delta^a_b + u^a u_b *)
DeltaUp[geo_, u_] := geo["ginv"] + Outer[Times, u, u];
DeltaMixed[geo_, u_] := IdentityMatrix[4] + Outer[Times, u, geo["g"] . u];

(* The grad-alpha-only row, as in derive_hq_2p1d.wls, for the C3/C5 comparison. *)
shippedRow[geo_, uv_, nuv_, alf_, kap_, tn_] := Module[{Dnu, acc, nda},
  Dnu = DcovV[geo, uv, nuv];
  acc = DcovV[geo, uv, uv];
  nda = Sum[(geo["g"] . acc)[[m]] nuv[[m]], {m, 4}];
  Table[tn (Dnu[[m]] - uv[[m]] nda) + nuv[[m]] +
    kap Sum[DeltaUp[geo, uv][[m, b]] D[alf, geo["x"][[b]]], {b, 4}], {m, 4}]];

(* The consistent row: (Ds/T) Delta^mu_lambda Nabla_nu T_Q^{nu lambda} + nu^mu. *)
consistentRow[geo_, uv_, nuv_, nfun_, Tfun_, hfun_, Dsv_] := Module[{P0, epsQ, TQ, dv},
  P0 = nfun Tfun;
  epsQ = nfun hfun - P0;
  TQ = Table[epsQ uv[[a]] uv[[b]] + P0 DeltaUp[geo, uv][[a, b]] +
             hfun (uv[[a]] nuv[[b]] + nuv[[a]] uv[[b]]), {a, 4}, {b, 4}];
  dv = divT2[geo, TQ];
  Table[(Dsv/Tfun) Sum[DeltaMixed[geo, uv][[m, l]] dv[[l]], {l, 4}] + nuv[[m]], {m, 4}]];

continuity[geo_, uv_, nuv_, nfun_] := divV[geo, nfun uv + nuv];

coefs[eqs_, vars_, cs_] := Module[{As, src},
  As = Table[Simplify @ Table[D[eqs[[i]], D[vars[[j]], c]], {i, Length[eqs]}, {j, Length[vars]}], {c, cs}];
  src = Simplify[eqs - Sum[As[[k]] . Table[D[v, cs[[k]]], {v, vars}], {k, Length[cs]}]];
  {As, src}];

(* ================================================== CHART 1: cylindrical control *)
Print["== G2/G3/C1  cylindrical control =="];
ucyl  = {Sqrt[1 + urb[tau, r]^2], urb[tau, r], 0, 0};
nucyl = {urb[tau, r] nur[tau, r]/Sqrt[1 + urb[tau, r]^2], nur[tau, r], 0, 0};
gate["G2cyl-uu", Simplify[Sum[(CYL["g"] . ucyl)[[m]] ucyl[[m]], {m, 4}] + 1] === 0, "u.u == -1"];
gate["G2cyl-un", Simplify[Sum[(CYL["g"] . ucyl)[[m]] nucyl[[m]], {m, 4}]] === 0, "u.nu == 0"];

nCyl = nn[Tb[tau, r], al[tau, r]];
hCyl = hq[Tb[tau, r]];
e2cC = consistentRow[CYL, ucyl, nucyl, nCyl, Tb[tau, r], hCyl, Ds];
e2cS = shippedRow[CYL, ucyl, nucyl, al[tau, r], kappa, taun];
e1c  = continuity[CYL, ucyl, nucyl, nCyl];

eqsC = {e1c, e2cC[[2]]}; varsC = {al[tau, r], nur[tau, r]};
sdC = Flatten @ Table[D[eqsC, D[v, c1], D[v2, c2]], {v, varsC}, {v2, varsC}, {c1, {tau, r}}, {c2, {tau, r}}];
gate["G3cyl", Simplify[sdC] === ConstantArray[0, Length[sdC]], "cylindrical consistent system is quasi-linear"];

(* C1: the at-rest, static, nu=0 source must be (Ds/T) d_r P0, the source of the
   1-D consistent row hqc_matrices.  That is the control. *)
statC = {urb[tau, r] -> 0, Derivative[1, 0][urb][tau, r] -> 0, Derivative[0, 1][urb][tau, r] -> 0,
         nur[tau, r] -> 0, Derivative[1, 0][nur][tau, r] -> 0, Derivative[0, 1][nur][tau, r] -> 0,
         Derivative[1, 0][Tb][tau, r] -> 0, Derivative[1, 0][al][tau, r] -> 0};
srcCyl = Simplify[e2cC[[2]] /. statC];
tgtCyl = Simplify[(Ds/Tb[tau, r]) (
    Tb[tau, r] Derivative[0, 1][nn][Tb[tau, r], al[tau, r]] Derivative[0, 1][al][tau, r] +
    (nn[Tb[tau, r], al[tau, r]] + Tb[tau, r] Derivative[1, 0][nn][Tb[tau, r], al[tau, r]]) *
      Derivative[0, 1][Tb][tau, r])];
gate["C1", Simplify[srcCyl - tgtCyl] === 0,
  "cylindrical static nu=0 source == (Ds/T) d_r P0  (== the 1-D hqc_matrices source)"];
If[Simplify[srcCyl - tgtCyl] =!= 0, Print["      diff = ", Simplify[srcCyl - tgtCyl]]];

(* ================================================ CHART 2: Cartesian 2+1D *)
Print["== G2/G3/C2-C5  Cartesian 2+1D =="];
ut2   = Sqrt[1 + ux2[tau, xx, yy]^2 + uy2[tau, xx, yy]^2];
ucar  = {ut2, ux2[tau, xx, yy], uy2[tau, xx, yy], 0};
nucar = {(ux2[tau, xx, yy] nux[tau, xx, yy] + uy2[tau, xx, yy] nuy[tau, xx, yy])/ut2,
         nux[tau, xx, yy], nuy[tau, xx, yy], 0};
gate["G2car-uu", Simplify[Sum[(CAR["g"] . ucar)[[m]] ucar[[m]], {m, 4}] + 1] === 0, "u.u == -1"];
gate["G2car-un", Simplify[Sum[(CAR["g"] . ucar)[[m]] nucar[[m]], {m, 4}]] === 0, "u.nu == 0"];

nCar = nn[T2[tau, xx, yy], al2[tau, xx, yy]];
hCar = hq[T2[tau, xx, yy]];
e2kC = consistentRow[CAR, ucar, nucar, nCar, T2[tau, xx, yy], hCar, Ds];
e2kS = shippedRow[CAR, ucar, nucar, al2[tau, xx, yy], kappa, taun];
e1k  = continuity[CAR, ucar, nucar, nCar];

eqsK  = {e1k, e2kC[[2]], e2kC[[3]]};
varsK = {al2[tau, xx, yy], nux[tau, xx, yy], nuy[tau, xx, yy]};
sdK = Flatten @ Table[D[eqsK, D[v, c1], D[v2, c2]], {v, varsK}, {v2, varsK},
                      {c1, {tau, xx, yy}}, {c2, {tau, xx, yy}}];
gate["G3car", Simplify[sdK] === ConstantArray[0, Length[sdK]],
  "Cartesian consistent system is quasi-linear"];

(* C2: at rest, static, nu = 0, the x-row source must be
   -( -kappa d_x alpha - (Ds/T)(n + T dn_dT) d_x T ) in the source-on-LHS
   convention, i.e. the closed-form Navier-Stokes limit. *)
statK = {ux2[tau, xx, yy] -> 0, uy2[tau, xx, yy] -> 0,
         Derivative[a_, b_, c_][ux2][tau, xx, yy] :> 0,
         Derivative[a_, b_, c_][uy2][tau, xx, yy] :> 0,
         nux[tau, xx, yy] -> 0, nuy[tau, xx, yy] -> 0,
         Derivative[a_, b_, c_][nux][tau, xx, yy] :> 0,
         Derivative[a_, b_, c_][nuy][tau, xx, yy] :> 0,
         Derivative[1, 0, 0][T2][tau, xx, yy] -> 0,
         Derivative[0, 0, 1][T2][tau, xx, yy] -> 0,
         Derivative[1, 0, 0][al2][tau, xx, yy] -> 0,
         Derivative[0, 0, 1][al2][tau, xx, yy] -> 0};
srcCar = Simplify[e2kC[[2]] /. statK];
tgtCar = Simplify[(Ds/T2[tau, xx, yy]) (
    T2[tau, xx, yy] Derivative[0, 1][nn][T2[tau, xx, yy], al2[tau, xx, yy]] *
      Derivative[0, 1, 0][al2][tau, xx, yy] +
    (nn[T2[tau, xx, yy], al2[tau, xx, yy]] +
     T2[tau, xx, yy] Derivative[1, 0][nn][T2[tau, xx, yy], al2[tau, xx, yy]]) *
      Derivative[0, 1, 0][T2][tau, xx, yy])];
gate["C2", Simplify[srcCar - tgtCar] === 0,
  "Cartesian at-rest NS source == kappa d_x alpha + (Ds/T)(n + T dn_dT) d_x T"];
If[Simplify[srcCar - tgtCar] =!= 0, Print["      diff = ", Simplify[srcCar - tgtCar]]];

(* C3: consistent - grad-alpha-only == the named terms, with tau_n = Ds h/T, kappa = Ds dn_da. *)
subTie = {taun -> Ds hq[T2[tau, xx, yy]]/T2[tau, xx, yy],
          kappa -> Ds Derivative[0, 1][nn][T2[tau, xx, yy], al2[tau, xx, yy]]};
DeltaU  = DeltaUp[CAR, ucar];
DeltaM  = DeltaMixed[CAR, ucar];
accCar  = DcovV[CAR, ucar, ucar];
thetaCar = divV[CAR, ucar];
DTcar   = Sum[ucar[[l]] D[T2[tau, xx, yy], CAR["x"][[l]]], {l, 4}];
gradTp  = Table[Sum[DeltaU[[m, b]] D[T2[tau, xx, yy], CAR["x"][[b]]], {b, 4}], {m, 4}];
nugradu = Table[Sum[DeltaM[[m, l]] Sum[nucar[[s]] (D[ucar[[l]], CAR["x"][[s]]] +
            Sum[CAR["Gam"][[l, s, q]] ucar[[q]], {q, 4}]), {s, 4}], {l, 4}], {m, 4}];
DlnHcar = Sum[ucar[[l]] D[hq[T2[tau, xx, yy]], CAR["x"][[l]]], {l, 4}];
taunS   = Ds hq[T2[tau, xx, yy]]/T2[tau, xx, yy];
nD      = nn[T2[tau, xx, yy], al2[tau, xx, yy]];
n10     = Derivative[1, 0][nn][T2[tau, xx, yy], al2[tau, xx, yy]];
extrasK = Table[
   (Ds/T2[tau, xx, yy]) (nD + T2[tau, xx, yy] n10) gradTp[[m]] +
   taunS nD accCar[[m]] +
   taunS nugradu[[m]] +
   nucar[[m]] (taunS thetaCar + (Ds/T2[tau, xx, yy]) DlnHcar), {m, 4}];
c3 = Simplify[(e2kC[[2]] - (e2kS[[2]] /. subTie)) - extrasK[[2]]];
gate["C3", c3 === 0, "consistent - shipped == the four named terms, exactly (x-row)"];
If[c3 =!= 0, Print["      residual = ", c3]];

(* C4: slab limit decouples nu^y *)
{AsK, srcK} = coefs[eqsK, varsK, {tau, xx, yy}];
{AtK, AxK, AyK} = AsK;
zeroY = {uy2[tau, xx, yy] -> 0, nuy[tau, xx, yy] -> 0,
         Derivative[a_, b_, c_][uy2][tau, xx, yy] :> 0, Derivative[a_, b_, c_][nuy][tau, xx, yy] :> 0,
         Derivative[a_, b_, 1][ux2][tau, xx, yy] :> 0, Derivative[a_, b_, 1][T2][tau, xx, yy] :> 0,
         Derivative[a_, b_, 1][al2][tau, xx, yy] :> 0, Derivative[a_, b_, 1][nux][tau, xx, yy] :> 0};
od[M_] := Simplify[{M[[1, 3]], M[[2, 3]], M[[3, 1]], M[[3, 2]]} /. zeroY];
gate["C4", Union @ Flatten[{od[AtK], od[AxK]}] === {0},
  "slab limit (u^y -> 0, d_y -> 0) decouples nu^y"];

(* C5: dropping the extras recovers the grad-alpha-only row identically *)
c5 = Simplify[(e2kC[[2]] - extrasK[[2]]) - (e2kS[[2]] /. subTie)];
gate["C5", c5 === 0, "consistent row minus the named terms == the shipped row"];

(* ------------------------------------------------------------------ export *)
jrepl = {
  ux2[tau, xx, yy] -> ux, uy2[tau, xx, yy] -> uy, T2[tau, xx, yy] -> T, al2[tau, xx, yy] -> al,
  nux[tau, xx, yy] -> nuxF, nuy[tau, xx, yy] -> nuyF,
  Derivative[1, 0, 0][ux2][tau, xx, yy] -> dtux, Derivative[0, 1, 0][ux2][tau, xx, yy] -> dxux,
  Derivative[0, 0, 1][ux2][tau, xx, yy] -> dyux,
  Derivative[1, 0, 0][uy2][tau, xx, yy] -> dtuy, Derivative[0, 1, 0][uy2][tau, xx, yy] -> dxuy,
  Derivative[0, 0, 1][uy2][tau, xx, yy] -> dyuy,
  Derivative[1, 0, 0][T2][tau, xx, yy] -> dtT, Derivative[0, 1, 0][T2][tau, xx, yy] -> dxT,
  Derivative[0, 0, 1][T2][tau, xx, yy] -> dyT,
  Derivative[1][hq][__] -> hp, hq[__] -> h,
  Derivative[0, 1][nn][__] -> dn_da, Derivative[1, 0][nn][__] -> dn_dT, nn[__] -> n};
jstr[e_] := StringReplace[ToString[InputForm[Simplify[Simplify[e] //. jrepl]]],
   {"Sqrt[" -> "sqrt(", "]" -> ")", "[" -> "("}];
mat2jl[M_, name_] := "    " <> name <> " = SMatrix{3,3}((\n" <>
   StringRiffle[Flatten @ Table["        " <> jstr[M[[i, j]]], {j, 3}, {i, 3}], ",\n"] <> "))\n";
vec2jl[V_, name_] := "    " <> name <> " = SVector{3}((\n" <>
   StringRiffle[Table["        " <> jstr[V[[i]]], {i, 3}], ",\n"] <> "))\n";

hdr = "# ============================================================================
# HQ_2p1d_BG_consistent_generated.jl -- GENERATED by
# deriv/derive_hq_consistent_2p1d.wls.  DO NOT EDIT.
#
# The consistent charm first moment phi = (alpha, nu^x, nu^y) on a prescribed
# 2+1D background, Milne (tau, x, y, eta_s), metric diag(-1, 1, 1, tau^2):
#
#     Nabla_mu (n u^mu + nu^mu) = 0
#     (D_s/T) Delta^mu_lambda Nabla_nu T_Q^{nu lambda} + nu^mu = 0
#
# with T_Q^{ab} = (n h - P0) u^a u^b + P0 Delta^{ab} + h (u^a nu^b + nu^a u^b),
# P0 = n T.  This is the full-grad-P row: its at-rest Navier-Stokes limit is
#
#     nu^i = -kappa d_i alpha - (D_s/T)(n + T dn/dT) d_i T
#
# The second term is absent from `hq2d_matrices` (the grad-alpha-only row).
# On a realistic Pb+Pb cell it is larger than the alpha drive.
#
# Convention: source-on-the-LHS, A_t d_tau phi + A_x d_x phi + A_y d_y phi + src = 0,
# as every matrix in this directory.
#
# Checks of the generator (prefixed HC-): C1 the cylindrical reduction reproduces the
# 1-D consistent source; C2 the at-rest NS limit is the closed form above; C3/C5
# consistent - grad-alpha-only is the named extra terms, identically.
#
# Known limitation: kappa enters these rows as Ds*dn_da, and the caller
# (HQ_2p1d_BG_m2.jl) passes dn_da after adding the regulator dmn_eps, so the
# regulator also enters kappa, as in the 1-D consistent rows (HQ_const_BG_consistent.jl).
# The effect is largest where n -> 0.
# ============================================================================

\"\"\"
    hq2d_matrices_consistent(al, nuxF, nuyF, tau, T, ux, uy, dtT, dxT, dyT,
                             dtux, dxux, dyux, dtuy, dxuy, dyuy,
                             n, dn_dT, dn_da, taun, kappa, h, hp, Ds)

(A_t, A_x, A_y, src) for the consistent 2+1D charm first moment.  Same signature
as `hq2d_matrices` plus the enthalpy per particle `h`, its derivative `hp`, and
`Ds`, the three quantities the full-grad-P channel needs.
\"\"\"
function hq2d_matrices_consistent(al, nuxF, nuyF, tau, T, ux, uy, dtT, dxT, dyT,
                                  dtux, dxux, dyux, dtuy, dxuy, dyuy,
                                  n, dn_dT, dn_da, taun, kappa, h, hp, Ds)
";
body = mat2jl[AtK, "At"] <> "\n" <> mat2jl[AxK, "Ax"] <> "\n" <> mat2jl[AyK, "Ay"] <> "\n" <>
       vec2jl[srcK, "src"];
ftr = "    return (At, Ax, Ay, src)\nend\n";

outfile = FileNameJoin[{PKG, "src", "Matrix", "HQ_2p1d_BG_consistent_generated.jl"}];
Export[outfile, hdr <> body <> ftr, "Text"];
Print["\n  wrote ", outfile];

(* reference states: the generated routine, evaluated from the derivation *)
SeedRandom[20261013];
rq[a_, b_, d_] := RandomInteger[{a, b}]/d;
rp = {};
Do[Module[{inC, sC},
    inC = <|al -> -rq[1, 80, 10], nuxF -> rq[-9, 9, 10^4], nuyF -> rq[-9, 9, 10^4], tau -> rq[5, 50, 10],
      T -> rq[15, 50, 100], ux -> rq[-6, 6, 10], uy -> rq[-6, 6, 10], dtT -> -rq[1, 99, 1000],
      dxT -> rq[-99, 99, 1000], dyT -> rq[-99, 99, 1000], dtux -> rq[-99, 99, 1000], dxux -> rq[-99, 99, 1000],
      dyux -> rq[-99, 99, 1000], dtuy -> rq[-99, 99, 1000], dxuy -> rq[-99, 99, 1000], dyuy -> rq[-99, 99, 1000],
      n -> rq[1, 99, 1000], taun -> rq[10, 60, 10], kappa -> rq[1, 99, 10^4],
      h -> rq[15, 30, 10], hp -> rq[1, 30, 10], Ds -> rq[1, 30, 100]|>;
    inC = Join[KeyValueMap[{SymbolName[#1], #1, #2} &, inC],
      {{"dn_dT", Verbatim[dn_dT] /. Verbatim -> Identity, rq[1, 99, 1000]},
       {"dn_da", Verbatim[dn_da] /. Verbatim -> Identity, rq[1, 99, 1000]}}];
    sC[e_] := Simplify[e //. jrepl];
    rp = Join[rp, refRows["hq2d_matrices_consistent", p, inC,
      Join[matRows["At", sC[AtK]], matRows["Ax", sC[AxK]], matRows["Ay", sC[AyK]], vecRows["src", sC[srcK]]]]]],
  {p, 12}];
If[fail == 0, writeRefpoints[FileNameJoin[{REFDIR, "refpoints_hq2d_consistent.csv"}], rp], Print["  NOT writing the reference states: a gate failed"]];
Print["\n  ", pass, "/", pass + fail, " gates passed"];
If[Environment["FLUIDUM_DERIV_RENDER"] === "1",   (* deriv/rendered/: PDF fragment + notebook *)
  Get[FileNameJoin[{DERIV, "render.wl"}]]; renderHQC[]];
If[fail > 0, Exit[1]];
