#!/usr/bin/env wolframscript
(* ============================================================================
   derive_hq_2p1d.wls -- the charm sector in 2+1 dimensions, symbolically.

   Derives the quasi-linear system for phi = (alpha, nu^x, nu^y) on a prescribed
   background in Milne coordinates (tau, x, y, eta), boost invariant, and exports
   it as Julia.

   Physics.  Charm is a passive tracer: N^mu = n u^mu + nu^mu with
       (1)  Nabla_mu N^mu = 0                                    [charge]
       (2)  tau_n Delta^mu_b D nu^b + nu^mu + kappa Delta^{mu b} d_b alpha = 0
   D = u.Nabla, Delta^{mu nu} = g^{mu nu} + u^mu u^nu, kappa = D_s dn/dalpha.
   Equation (2) is the grad-alpha-only first moment of HQ_const_BG.jl; its
   Navier-Stokes limit is Fick's law nu = -kappa grad_perp alpha.

   Two charts.  The same routine is run in cylindrical Milne (tau,r,phi,eta) with
   azimuthal symmetry, where it must reproduce the hand-written
   `one_d_viscous_matrix_fugacity_BG_only`, and then in Cartesian Milne (tau,x,y,eta),
   which is the result.  The 1-D case is derived here only as that control.

   Checks (printed with the prefix H-)
     G0  metric compatibility Nabla_lambda g_{mu nu} == 0, both charts
     G1  divergence check: Nabla.V for V=(X(tau),0,0,0) must be X' + X/tau
     G2  transversality u.nu == 0 and normalisation u.u == -1, identically
     G3  the system is quasi-linear: no second derivatives of the fields survive
     G4  cylindrical reduction == the hand-written 1-D matrix, entry by entry, up to
         the overall Jacobian u^tau^3 r tau (against a transcription of it, below)
     G4src  the same, against HQ_const_BG.jl read from the file (julia_source_reader.wl),
         so that the comparison follows any edit of the code
     G5  the Cartesian slab limit (u^y -> 0, d_y -> 0) decouples nu^y
     G6a/G6b  the live-background (10-field) rows are the same system with the
         medium derivatives moved from the source into the matrices

   Run:    wolframscript -file deriv/derive_hq_2p1d.wls
   Writes: src/Matrix/HQ_2p1d_BG_generated.jl, and refpoints_hq2d.csv (the derivation
           evaluated at 12 random states) to $FLUIDUM_REFPOINTS_DIR (default a
           temporary directory)
   ============================================================================ *)

$HistoryLength = 0;
pass = 0; fail = 0;
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], "H-" <> id, "  ", msg]);
PKG = ParentDirectory[DirectoryName[$InputFileName]];   (* the Fluidum.jl package root *)
DERIV = DirectoryName[$InputFileName];
REFDIR = With[{v = Environment["FLUIDUM_REFPOINTS_DIR"]},
  If[StringQ[v] && v =!= "", v, FileNameJoin[{$TemporaryDirectory, "fluidum_refpoints"}]]];
Get[FileNameJoin[{DERIV, "julia_source_reader.wl"}]];
Get[FileNameJoin[{DERIV, "refpoints_csv.wl"}]];

mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam|>];

divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] +
  Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];

DcovV[geo_, u_, W_] := Table[Sum[u[[l]] D[W[[m]], geo["x"][[l]]], {l, 4}] +
  Sum[u[[l]] geo["Gam"][[m, l, s]] W[[s]], {l, 4}, {s, 4}], {m, 4}];

metricCompat[geo_] := Simplify @ Table[
  D[geo["g"][[m, nn2]], geo["x"][[l]]] - Sum[geo["Gam"][[s, l, m]] geo["g"][[s, nn2]], {s, 4}]
  - Sum[geo["Gam"][[s, l, nn2]] geo["g"][[m, s]], {s, 4}], {l, 4}, {m, 4}, {nn2, 4}];

CYL = mkGeom[{tau, r, phi, eta}, {-1, 1, r^2, tau^2}];
CAR = mkGeom[{tau, xx, yy, eta}, {-1, 1, 1, tau^2}];

Print["== G0/G1  geometry =="];
gate["G0cyl", metricCompat[CYL] === ConstantArray[0, {4, 4, 4}], "cylindrical Milne: Nabla g == 0"];
gate["G0car", metricCompat[CAR] === ConstantArray[0, {4, 4, 4}], "Cartesian   Milne: Nabla g == 0"];
gate["G1cyl", Simplify[divV[CYL, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "cylindrical: Nabla.V = X' + X/tau for V=(X(tau),0,0,0)"];
gate["G1car", Simplify[divV[CAR, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "Cartesian:   Nabla.V = X' + X/tau"];

(* the two equations, for any chart *)
buildEqs[geo_, uvec_, nuvec_, alf_, nfun_, kap_, tn_] := Module[{Dnu, acc, nda, e1, e2},
  Dnu = DcovV[geo, uvec, nuvec];
  acc = DcovV[geo, uvec, uvec];
  nda = Sum[(geo["g"] . acc)[[m]] nuvec[[m]], {m, 4}];
  e1 = divV[geo, nfun uvec + nuvec];
  e2 = Table[tn (Dnu[[m]] - uvec[[m]] nda) + nuvec[[m]] +
       kap Sum[(geo["ginv"] + Outer[Times, uvec, uvec])[[m, b]] D[alf, geo["x"][[b]]], {b, 4}], {m, 4}];
  {e1, e2}];

coefs[eqs_, vars_, cs_] := Module[{As, src},
  As = Table[Simplify @ Table[D[eqs[[i]], D[vars[[j]], c]], {i, Length[eqs]}, {j, Length[vars]}], {c, cs}];
  src = Simplify[eqs - Sum[As[[k]] . Table[D[v, cs[[k]]], {v, vars}], {k, Length[cs]}]];
  {As, src}];

(* ---------------------------------------- CHART 1: cylindrical (the control) *)
Print["== G2/G3/G4  cylindrical control =="];
ucyl  = {Sqrt[1 + urb[tau, r]^2], urb[tau, r], 0, 0};
nucyl = {urb[tau, r] nur[tau, r]/Sqrt[1 + urb[tau, r]^2], nur[tau, r], 0, 0};
gate["G2cyl-uu", Simplify[Sum[(CYL["g"] . ucyl)[[m]] ucyl[[m]], {m, 4}] + 1] === 0, "u.u == -1"];
gate["G2cyl-un", Simplify[Sum[(CYL["g"] . ucyl)[[m]] nucyl[[m]], {m, 4}]] === 0, "u.nu == 0"];
{e1c, e2c} = buildEqs[CYL, ucyl, nucyl, al[tau, r], nn[Tb[tau, r], al[tau, r]], kappa, taun];
eqsC = {e1c, e2c[[2]]}; varsC = {al[tau, r], nur[tau, r]};
sdC = Flatten @ Table[D[eqsC, D[v, c1], D[v2, c2]], {v, varsC}, {v2, varsC}, {c1, {tau, r}}, {c2, {tau, r}}];
gate["G3cyl", Simplify[sdC] === ConstantArray[0, Length[sdC]], "cylindrical system is quasi-linear"];
{AsC, srcC} = coefs[eqsC, varsC, {tau, r}];
{AtC, AxC} = AsC;

shippedAt = {{(1 + ur^2)^2 r tau dnmu,        ur (1 + ur^2) r tau},
             {kappa ur (1 + ur^2)^2 r tau,    (1 + ur^2)^2 taun r tau}};
shippedAx = {{ur (1 + ur^2)^(3/2) r tau dnmu, (1 + ur^2)^(3/2) r tau},
             {kappa (1 + ur^2)^(5/2) r tau,   ur (1 + ur^2)^(3/2) taun r tau}};
shippedSrc = {
  nuF (Sqrt[1 + ur^2] tau + ur^2 Sqrt[1 + ur^2] tau + r (ur + ur^3 + dtur tau)) +
    (1 + ur^2) nF (ur Sqrt[1 + ur^2] tau + r (1 + ur^2 + dtur ur tau + drur Sqrt[1 + ur^2] tau)) +
    (1 + ur^2) (dtT + dtT ur^2 + drT ur Sqrt[1 + ur^2]) r tau dnT,
  (Sqrt[1 + ur^2] - ur (dtur taun + dtur ur^2 taun + ur Sqrt[1 + ur^2] (-1 + drur taun))) r nuF tau};

(* Two passes on purpose: a single rule list would try to rewrite Derivative[0,1][nn][Tb[tau,r],
   al[tau,r]] before its arguments have been renamed, so the nn rules would never fire and the
   comparison would fail on notation.  Pass 1 renames the background, pass 2 the EoS calls,
   matched on any arguments. *)
subBG = {Derivative[1, 0][urb][tau, r] -> dtur, Derivative[0, 1][urb][tau, r] -> drur, urb[tau, r] -> ur,
         Derivative[1, 0][Tb][tau, r] -> dtT, Derivative[0, 1][Tb][tau, r] -> drT, Tb[tau, r] -> TT,
         nur[tau, r] -> nuF, al[tau, r] -> alF};
subEOS = {Derivative[0, 1][nn][__] -> dnmu, Derivative[1, 0][nn][__] -> dnT, nn[__] -> nF};
red[e_] := Simplify[(e /. subBG) /. subEOS];
jac = (1 + ur^2)^(3/2) r tau;
dAt = Simplify[red[AtC] jac - shippedAt];
dAx = Simplify[red[AxC] jac - shippedAx];
dSr = Simplify[red[srcC] jac - shippedSrc];
gate["G4-At",  dAt === {{0, 0}, {0, 0}}, "cylindrical A_t * u^tau^3 r tau == Fluidum's shipped A_t"];
gate["G4-Ax",  dAx === {{0, 0}, {0, 0}}, "cylindrical A_x * u^tau^3 r tau == Fluidum's shipped A_x"];
gate["G4-src", dSr === {0, 0},           "cylindrical src  * u^tau^3 r tau == Fluidum's shipped source"];
If[dAt =!= {{0,0},{0,0}}, Print["      dAt  = ", dAt]];
If[dAx =!= {{0,0},{0,0}}, Print["      dAx  = ", dAx]];
If[dSr =!= {0,0},         Print["      dsrc = ", dSr]];
(* H-G4src: the same comparison against the .jl file *)
shF = juliaMatrices[FileNameJoin[{PKG, "src", "Matrix", "HQ_const_BG.jl"}],
  "one_d_viscous_matrix_fugacity_BG_only"] /. {jl["X"][1] -> alF, jl["X"][2] -> nuF, jl["tau"] -> tau,
  jl["r"] -> r, jl["ur"] -> ur, jl["dtT"] -> dtT, jl["drT"] -> drT, jl["drur"] -> drur, jl["dtur"] -> dtur,
  jl["n"] -> nF, jl["dnUdmu"] -> dnmu, jl["dnUdT"] -> dnT, jl["tauDiffusion"] -> taun, jl["kappaG"] -> kappa};
gate["G4src", Simplify[{red[AtC] jac - shF["At"], red[AxC] jac - shF["Ax"], red[srcC] jac - shF["source"]}] ===
  {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {0, 0}},
  "cylindrical A_t, A_x, src * u^tau^3 r tau == HQ_const_BG.jl READ FROM THE FILE"];

(* ---------------------------------------- CHART 2: Cartesian 2+1D (the result) *)
Print["== G2/G3/G5  Cartesian 2+1D =="];
ut2   = Sqrt[1 + ux2[tau, xx, yy]^2 + uy2[tau, xx, yy]^2];
ucar  = {ut2, ux2[tau, xx, yy], uy2[tau, xx, yy], 0};
nucar = {(ux2[tau, xx, yy] nux[tau, xx, yy] + uy2[tau, xx, yy] nuy[tau, xx, yy])/ut2,
         nux[tau, xx, yy], nuy[tau, xx, yy], 0};
gate["G2car-uu", Simplify[Sum[(CAR["g"] . ucar)[[m]] ucar[[m]], {m, 4}] + 1] === 0, "u.u == -1"];
gate["G2car-un", Simplify[Sum[(CAR["g"] . ucar)[[m]] nucar[[m]], {m, 4}]] === 0, "u.nu == 0"];
{e1k, e2k} = buildEqs[CAR, ucar, nucar, al2[tau, xx, yy], nn[T2[tau, xx, yy], al2[tau, xx, yy]], kappa, taun];
eqsK = {e1k, e2k[[2]], e2k[[3]]};
varsK = {al2[tau, xx, yy], nux[tau, xx, yy], nuy[tau, xx, yy]};
sdK = Flatten @ Table[D[eqsK, D[v, c1], D[v2, c2]], {v, varsK}, {v2, varsK}, {c1, {tau, xx, yy}}, {c2, {tau, xx, yy}}];
gate["G3car", Simplify[sdK] === ConstantArray[0, Length[sdK]], "Cartesian system is quasi-linear"];
{AsK, srcK} = coefs[eqsK, varsK, {tau, xx, yy}];
{AtK, AxK, AyK} = AsK;

slab = {Derivative[a_, b_, c_] /; c > 0 -> (0 &)};
zeroY = {uy2[tau, xx, yy] -> 0, nuy[tau, xx, yy] -> 0,
         Derivative[a_, b_, c_][uy2][tau, xx, yy] :> 0, Derivative[a_, b_, c_][nuy][tau, xx, yy] :> 0,
         Derivative[a_, b_, 1][ux2][tau, xx, yy] :> 0, Derivative[a_, b_, 1][T2][tau, xx, yy] :> 0,
         Derivative[a_, b_, 1][al2][tau, xx, yy] :> 0, Derivative[a_, b_, 1][nux][tau, xx, yy] :> 0};
od[M_] := Simplify[{M[[1, 3]], M[[2, 3]], M[[3, 1]], M[[3, 2]]} /. zeroY];
gate["G5", Union @ Flatten[{od[AtK], od[AxK]}] === {0},
  "slab limit (u^y -> 0, d_y -> 0) decouples nu^y from (alpha, nu^x)"];

(* ----------------------- CHART 2b: the same equations with T, u^i as fields -- *)
(* For the live-background (10-field) solver the medium is evolved, so d_tau T, d_x u^x, ...
   are unknowns: they move from the source into the matrices, as extra columns.  The system
   is the same, re-sorted; G6 checks this. *)
Print["== G6  the live-background (10-field) re-sorting =="];
varsL = {T2[tau, xx, yy], ux2[tau, xx, yy], uy2[tau, xx, yy],
         al2[tau, xx, yy], nux[tau, xx, yy], nuy[tau, xx, yy]};
{AsL, srcL} = coefs[eqsK, varsL, {tau, xx, yy}];
{AtL, AxL, AyL} = AsL;
(* G6a: the last three columns are exactly the background-mode 3x3 matrices *)
gate["G6a", Simplify[AtL[[All, 4 ;; 6]] - AtK] === ConstantArray[0, {3, 3}] &&
            Simplify[AxL[[All, 4 ;; 6]] - AxK] === ConstantArray[0, {3, 3}] &&
            Simplify[AyL[[All, 4 ;; 6]] - AyK] === ConstantArray[0, {3, 3}],
  "the (alpha, nu^x, nu^y) columns of the 10-field rows == the background-mode matrices"];
(* G6b: putting the medium columns back onto the background derivatives recovers the BG source *)
bgput = AtL[[All, 1 ;; 3]] . {D[T2[tau,xx,yy],tau], D[ux2[tau,xx,yy],tau], D[uy2[tau,xx,yy],tau]} +
        AxL[[All, 1 ;; 3]] . {D[T2[tau,xx,yy],xx],  D[ux2[tau,xx,yy],xx],  D[uy2[tau,xx,yy],xx]} +
        AyL[[All, 1 ;; 3]] . {D[T2[tau,xx,yy],yy],  D[ux2[tau,xx,yy],yy],  D[uy2[tau,xx,yy],yy]};
gate["G6b", Simplify[srcL + bgput - srcK] === {0, 0, 0},
  "src_live + (medium columns).(background derivatives) == src_background: same system, re-sorted"];

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
  Derivative[0, 1][nn][__] -> dn_da, Derivative[1, 0][nn][__] -> dn_dT, nn[__] -> n};
jstr[e_] := StringReplace[ToString[InputForm[Simplify[Simplify[e] //. jrepl]]],
   {"Sqrt[" -> "sqrt(", "]" -> ")", "[" -> "("}];
mat2jl[M_, name_] := "    " <> name <> " = SMatrix{3,3}((\n" <>
   StringRiffle[Flatten @ Table["        " <> jstr[M[[i, j]]], {j, 3}, {i, 3}], ",\n"] <> "))\n";
vec2jl[V_, name_] := "    " <> name <> " = SVector{3}((\n" <>
   StringRiffle[Table["        " <> jstr[V[[i]]], {i, 3}], ",\n"] <> "))\n";

hdr = "# ============================================================================
# HQ_2p1d_BG_generated.jl -- GENERATED by deriv/derive_hq_2p1d.wls.  DO NOT EDIT.
#
# The charm first moment phi = (alpha, nu^x, nu^y) on a prescribed 2+1D background,
# Milne (tau, x, y, eta_s), boost invariant, metric diag(-1, 1, 1, tau^2):
#     Nabla_mu (n u^mu + nu^mu) = 0
#     tau_n Delta^mu_b D nu^b + nu^mu + kappa Delta^{mu b} d_b alpha = 0
# in the source-on-the-LHS convention  A_t d_tau phi + A_x d_x phi + A_y d_y phi + src = 0,
# the convention of every matrix in this directory.
#
# The same generator reproduces the hand-written 1-D cylindrical matrix of HQ_const_BG.jl
# entry by entry (check H-G4, and H-G4src against the file itself) before it is applied
# to the Cartesian chart.
# ============================================================================

\"\"\"
    hq2d_matrices(al, nuxF, nuyF, tau, T, ux, uy, dtT, dxT, dyT,
                  dtux, dxux, dyux, dtuy, dxuy, dyuy, n, dn_dT, dn_da, taun, kappa)

(A_t, A_x, A_y, src) for the 2+1D charm first moment on a prescribed background.
\"\"\"
function hq2d_matrices(al, nuxF, nuyF, tau, T, ux, uy, dtT, dxT, dyT,
                       dtux, dxux, dyux, dtuy, dxuy, dyuy,
                       n, dn_dT, dn_da, taun, kappa)
";
body = mat2jl[AtK, "At"] <> "\n" <> mat2jl[AxK, "Ax"] <> "\n" <> mat2jl[AyK, "Ay"] <> "\n" <> vec2jl[srcK, "src"];

mat36jl[M_, name_] := "    " <> name <> " = SMatrix{3,6}((\n" <>
   StringRiffle[Flatten @ Table["        " <> jstr[M[[i, j]]], {j, 6}, {i, 3}], ",\n"] <> "))\n";
hdr2 = "\n\n\"\"\"
    hq2d_live_rows(al, nuxF, nuyF, tau, T, ux, uy, n, dn_dT, dn_da, taun, kappa)

The same three charm equations with the medium promoted from background to fields: the returned
3x6 blocks are the coefficients of d(T, u^x, u^y, alpha, nu^x, nu^y), in that column order, and
`src` is what is left over.  Columns 4:6 are the background-mode 3x3 matrices, and putting
columns 1:3 back onto the background derivatives recovers the background-mode source (checks
H-G6a/H-G6b of deriv/derive_hq_2p1d.wls).  The two forms are one system, sorted two ways.

Used by `matrix2d_visc_HQ_BG!` to fill rows 8:10 of the 10-field live solver, whose columns are
(T, ux, uy, piyy, pizz, pixy, piB, alpha, nu^x, nu^y): the six columns land at 1,2,3,8,9,10.
\"\"\"
function hq2d_live_rows(al, nuxF, nuyF, tau, T, ux, uy, n, dn_dT, dn_da, taun, kappa)
";
body2 = mat36jl[AtL, "AtL"] <> "\n" <> mat36jl[AxL, "AxL"] <> "\n" <> mat36jl[AyL, "AyL"] <> "\n" <> vec2jl[srcL, "srcL"];
out = hdr <> body <> "\n    return (At, Ax, Ay, src)\nend\n" <>
      hdr2 <> body2 <> "\n    return (AtL, AxL, AyL, srcL)\nend\n";
path = FileNameJoin[{PKG, "src", "Matrix", "HQ_2p1d_BG_generated.jl"}];
Export[path, out, "Text"];
Print["\n  wrote ", path];

(* reference states: both generated routines, evaluated from the derivation *)
SeedRandom[20261013];
rq[a_, b_, d_] := RandomInteger[{a, b}]/d;
rp = {};
Do[Module[{inB, inL, sB},
    inB = <|al -> -rq[1, 80, 10], nuxF -> rq[-9, 9, 10^4], nuyF -> rq[-9, 9, 10^4], tau -> rq[5, 50, 10],
      T -> rq[15, 50, 100], ux -> rq[-6, 6, 10], uy -> rq[-6, 6, 10], dtT -> -rq[1, 99, 1000],
      dxT -> rq[-99, 99, 1000], dyT -> rq[-99, 99, 1000], dtux -> rq[-99, 99, 1000], dxux -> rq[-99, 99, 1000],
      dyux -> rq[-99, 99, 1000], dtuy -> rq[-99, 99, 1000], dxuy -> rq[-99, 99, 1000], dyuy -> rq[-99, 99, 1000],
      n -> rq[1, 99, 1000], taun -> rq[10, 60, 10], kappa -> rq[1, 99, 10^4]|>;
    (* dn_dT / dn_da are written as patterns by jrepl (that is how their names reach the .jl) *)
    inB = Join[KeyValueMap[{SymbolName[#1], #1, #2} &, inB],
      {{"dn_dT", Verbatim[dn_dT] /. Verbatim -> Identity, rq[1, 99, 1000]},
       {"dn_da", Verbatim[dn_da] /. Verbatim -> Identity, rq[1, 99, 1000]}}];
    sB[e_] := e //. jrepl;
    rp = Join[rp, refRows["hq2d_matrices", p, inB,
      Join[matRows["At", sB[AtK]], matRows["Ax", sB[AxK]], matRows["Ay", sB[AyK]], vecRows["src", sB[srcK]]]]];
    inL = Select[inB, MemberQ[{"al", "nuxF", "nuyF", "tau", "T", "ux", "uy", "n", "dn_dT", "dn_da", "taun", "kappa"}, #[[1]]] &];
    rp = Join[rp, refRows["hq2d_live_rows", p, inL,
      Join[matRows["AtL", sB[AtL]], matRows["AxL", sB[AxL]], matRows["AyL", sB[AyL]], vecRows["srcL", sB[srcL]]]]]],
  {p, 12}];
If[fail == 0, writeRefpoints[FileNameJoin[{REFDIR, "refpoints_hq2d.csv"}], rp], Print["  NOT writing the reference states: a gate failed"]];
Print["\n", pass, " passed, ", fail, " failed  (", pass + fail, " gates evaluated)"];
Exit[If[fail == 0, 0, 1]];
