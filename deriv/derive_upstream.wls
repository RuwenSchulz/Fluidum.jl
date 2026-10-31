#!/usr/bin/env wolframscript
(* ============================================================================
   derive_upstream.wls -- the matrices no solver path of this branch relies on, against
   the same tensor equations.  two_d_viscous_matrix stays wrong on purpose: `legacy = true` is its reason to exist.

     two_d_viscous_matrix                   (2d_viscous.jl)                   medium, opt-in legacy = true
     two_d_viscous_HQ_matrix                (2d_viscous_fugacity.jl)          medium + charm, 10 fields
     one_d_ideal_matrix_ruwen7              (1d_ideal_HQ_cilindrical_p_br_fugacity.jl)  ideal, charm in T^{mu nu}
     one_d_ideal_matrix_ruwen6              (1d_ideal_HQ_cilindrical_p_br_density.jl)   the same in (T, u^r, n)

   Physics of the p_br files: P_tot = P(T) + P_h(T, alpha), e + P_tot = T d_T P_tot at fixed alpha
   (which is Ts + mu n, since mu = alpha T), T^{mu nu} = (e + P_tot) uu + P_tot g,
   Nabla_mu (n u^mu) = 0.  ruwen6 takes n = d_alpha P_h / T as the field and eliminates alpha.

   Checks (printed with the prefix U-)
     U-L*   legacy against the derived matrix read from its generated block: the entries that differ
     U-H*   the 10-field medium block is the derived matrix at dpp = 0; its charm rows are Fick with kappa = Ds n
     U-P7*  ruwen7 raw, entry by entry
     U-P6*  ruwen6 raw

   Run:    wolframscript -file deriv/derive_upstream.wls
   Writes: nothing into src/
   ============================================================================ *)

$HistoryLength = 0;
pass = 0; fail = 0; gateLog = {};
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  AppendTo[gateLog, {"U-" <> id, TrueQ[ok], msg}];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], "U-" <> id, "  ", msg]);
PKG = ParentDirectory[DirectoryName[$InputFileName]];
DERIV = DirectoryName[$InputFileName];
MAT = FileNameJoin[{PKG, "src", "Matrix"}];
Get[FileNameJoin[{DERIV, "julia_source_reader.wl"}]];

mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam|>];
divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] + Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];
DcovV[geo_, u_, W_] := Table[Sum[u[[l]] D[W[[m]], geo["x"][[l]]], {l, 4}] +
  Sum[u[[l]] geo["Gam"][[m, l, s]] W[[s]], {l, 4}, {s, 4}], {m, 4}];
divT[geo_, P_] := Table[Sum[D[P[[m, n2]], geo["x"][[m]]], {m, 4}]
   + Sum[geo["Gam"][[m, m, s]] P[[s, n2]], {m, 4}, {s, 4}]
   + Sum[geo["Gam"][[n2, m, s]] P[[m, s]], {m, 4}, {s, 4}], {n2, 4}];
coefs[eqs_, vars_, cs_] := Module[{As, src},
  As = Table[Simplify @ Table[D[eqs[[i]], D[vars[[j]], c]], {i, Length[eqs]}, {j, Length[vars]}], {c, cs}];
  src = Simplify[eqs - Sum[As[[k]] . Table[D[v, cs[[k]]], {v, vars}], {k, Length[cs]}]];
  {As, src}];
CYL = mkGeom[{tau, rr, ph, et}, {-1, 1, rr^2, tau^2}];
CAR = mkGeom[{tau, xx, yy, et}, {-1, 1, 1, tau^2}];
nz[d_] := Position[d, x_ /; x =!= 0, {2}, Heads -> False];

(* ------------------------------------------------------- the legacy medium *)
Print["== legacy two_d_viscous_matrix =="];
fL = {T, ux, uy, pyy, pee, pxy, PB};
rlCommon = {jl["tau"] -> tau, jl["dtp"] -> p1, jl["dtdtp"] -> p2, jl["zeta"] -> zeta, jl["visc"] -> visc,
  jl["tauS"] -> tauS, jl["tauB"] -> tauB};
rlLeg = Join[Table[jl["u"][k] -> fL[[k]], {k, 7}], rlCommon, {jl["p"] -> p0}];
rlDer = Join[Thread[jl /@ {"T", "ux", "uy", "piyy", "pizz", "pixy", "piB"} -> fL], rlCommon, {jl["pres"] -> p0}];
leg = juliaMatrices[FileNameJoin[{MAT, "2d_viscous.jl"}], "two_d_viscous_matrix"] /. rlLeg;
der = juliaMatrices[FileNameJoin[{MAT, "2d_viscous.jl"}], "two_d_viscous_matrix_derived"] /. rlDer;
dLeg = AssociationMap[Simplify[leg[#] - der[#]] &, {"At", "Ax", "Ay"}];
wrong = Flatten[Table[{k, #} & /@ nz[dLeg[k]], {k, {"At", "Ax", "Ay"}}], 1];
gate["L-count", Length[wrong] == 10, "exactly ten of the 147 entries of A_t, A_x, A_y differ from the derived matrix (dpp = 0)"];
gate["L-place", wrong === {{"At", {4, 3}}, {"At", {6, 2}}, {"At", {6, 3}}, {"Ax", {4, 3}}, {"Ax", {6, 2}}, {"Ax", {6, 3}},
    {"Ay", {4, 2}}, {"Ay", {4, 3}}, {"Ay", {6, 2}}, {"Ay", {6, 3}}},
  "they are in the velocity columns of the pi^{yy} and pi^{xy} rows"];
gate["L-src", Simplify[leg["source"] - der["source"]] === ConstantArray[0, 7], "the source is right"];

(* --------------------------------------------------- the 10-field upstream *)
Print["== two_d_viscous_HQ_matrix =="];
fH = Join[fL, {al, nx, ny}];
hq = juliaMatrices[FileNameJoin[{MAT, "2d_viscous_fugacity.jl"}], "two_d_viscous_HQ_matrix"] /.
  Join[Table[jl`getindex[jl["u"], k] -> fH[[k]], {k, 10}], Table[jl["u"][k] -> fH[[k]], {k, 10}], rlCommon,
    Thread[jl /@ {"T", "ux", "uy", "piyy", "pizz", "pixy", "piB"} -> fL],
    {jl["p"] -> p0, jl["n"] -> n, jl["dtn"] -> nT, jl["dmn"] -> na, jl["tauDiff"] -> taun, jl["Ds"] -> Ds}];
gate["H-med", Simplify[{hq["At"][[1 ;; 7, 1 ;; 7]] - der["At"], hq["Ax"][[1 ;; 7, 1 ;; 7]] - der["Ax"], hq["Ay"][[1 ;; 7, 1 ;; 7]] - der["Ay"]}] ===
    ConstantArray[0, {3, 7, 7}] && Simplify[hq["source"][[1 ;; 7]] - der["source"]] === ConstantArray[0, 7],
  "its medium block == the derived matrix at dpp = 0, raw (fixed 2026-10-30: it carried the ten legacy entries)"];
gate["H-pass", hq["At"][[1 ;; 7, 8 ;; 10]] === ConstantArray[0, {7, 3}] && hq["Ax"][[1 ;; 7, 8 ;; 10]] === ConstantArray[0, {7, 3}] &&
   hq["Ay"][[1 ;; 7, 8 ;; 10]] === ConstantArray[0, {7, 3}], "the medium does not see charm"];
utK = Sqrt[1 + ux2[tau, xx, yy]^2 + uy2[tau, xx, yy]^2];
uK = {utK, ux2[tau, xx, yy], uy2[tau, xx, yy], 0};
nuK = {(ux2[tau, xx, yy] nux[tau, xx, yy] + uy2[tau, xx, yy] nuy[tau, xx, yy])/utK, nux[tau, xx, yy], nuy[tau, xx, yy], 0};
accK = DcovV[CAR, uK, uK]; DnuK = DcovV[CAR, uK, nuK];
ndaK = Sum[(CAR["g"] . accK)[[m]] nuK[[m]], {m, 4}];
e1K = divV[CAR, nn[T2[tau, xx, yy], al2[tau, xx, yy]] uK + nuK];
e2K = Table[taun (DnuK[[m]] - uK[[m]] ndaK) + nuK[[m]] +
   Ds nn[T2[tau, xx, yy], al2[tau, xx, yy]] Sum[(CAR["ginv"] + Outer[Times, uK, uK])[[m, b]] D[al2[tau, xx, yy], CAR["x"][[b]]], {b, 4}], {m, 4}];
eqsH = {e1K, e2K[[2]], e2K[[3]]};
varsH = {T2[tau, xx, yy], ux2[tau, xx, yy], uy2[tau, xx, yy], al2[tau, xx, yy], nux[tau, xx, yy], nuy[tau, xx, yy]};
{AsH, srcH} = coefs[eqsH, varsH, {tau, xx, yy}];
subH = {T2[tau, xx, yy] -> T, ux2[tau, xx, yy] -> ux, uy2[tau, xx, yy] -> uy, al2[tau, xx, yy] -> al, nux[tau, xx, yy] -> nx, nuy[tau, xx, yy] -> ny,
  Derivative[1, 0][nn][__] -> nT, Derivative[0, 1][nn][__] -> na, nn[__] -> n};
colsH = {1, 2, 3, 8, 9, 10};
gate["H-charm", Simplify[Table[(AsH[[k]] /. subH) - hq[{"At", "Ax", "Ay"}[[k]]][[8 ;; 10, colsH]], {k, 3}]] === ConstantArray[0, {3, 3, 6}],
  "rows 8-10 == N and Fick, tau_n Delta D nu + nu + Ds n Delta grad(alpha) = 0, raw"];
gate["H-charmsrc", Simplify[(srcH /. subH) - hq["source"][[8 ;; 10]]] === {0, 0, 0}, "and their source, raw"];
gate["H-kappa", !FreeQ[hq["At"], Ds n], "the argument Ds enters as kappa = Ds n, which is D_s dn/dalpha for one unit-charge Boltzmann species"];

(* ------------------------------------------------------ the p_br kernels *)
Print["== one_d_ideal_matrix_ruwen7, ruwen6 =="];
utC = Sqrt[1 + ur[tau, rr]^2];
uC = {utC, ur[tau, rr], 0, 0};
tmnP[Tf_, af_] := Module[{Ptot = pp[Tf] + ph[Tf, af], W = Tf (pp'[Tf] + Derivative[1, 0][ph][Tf, af])},
  Table[W uC[[m]] uC[[b]] + Ptot CYL["ginv"][[m, b]], {m, 4}, {b, 4}]];
subP = {TT[tau, rr] -> T, ur[tau, rr] -> u, al[tau, rr] -> al0, YY[tau, rr] -> Y,
  Derivative[2][pp][_] -> p2, Derivative[1][pp][_] -> p1, pp[_] -> p0,
  Derivative[2, 0][ph][__] -> pTT, Derivative[1, 1][ph][__] -> pTa, Derivative[0, 2][ph][__] -> paa,
  Derivative[1, 0][ph][__] -> pT, Derivative[0, 1][ph][__] -> pa, ph[__] -> ph0,
  Derivative[1, 0][nn][__] -> nT, Derivative[0, 1][nn][__] -> na, nn[__] -> n, rr -> R};
rlP = {jl["X"][1] -> T, jl["X"][2] -> u, jl["tau"] -> tau, jl["r"] -> R, jl["dPUdT"] -> p1, jl["dPUdTdT"] -> p2, jl["dPhqUdT"] -> pT,
  jl["dPhqUdalpha"] -> pa, jl["ddPhqUdTdT"] -> pTT, jl["ddPhqUdTdalpha"] -> pTa, jl["ddPhqUdalphadalpha"] -> paa, jl["n"] -> n,
  jl["dnUdmu"] -> na, jl["dnUdT"] -> nT};
eqs7 = Join[divT[CYL, tmnP[TT[tau, rr], al[tau, rr]]][[1 ;; 2]], {divV[CYL, nn[TT[tau, rr], al[tau, rr]] uC]}];
vars7 = {TT[tau, rr], ur[tau, rr], al[tau, rr]};
{As7, src7} = coefs[eqs7, vars7, {tau, rr}];
{At7, Ar7, S7} = {As7[[1]], As7[[2]], src7} //. subP;
r7 = juliaMatrices[FileNameJoin[{MAT, "1d_ideal_HQ_cilindrical_p_br_fugacity.jl"}], "one_d_ideal_matrix_ruwen7"] /.
  Join[rlP, {jl["X"][3] -> al0}];
gate["P7", Simplify[{At7 - r7["At"], Ar7 - r7["Ax"]}] === ConstantArray[0, {2, 3, 3}] && Simplify[S7 - r7["source"]] === {0, 0, 0},
  "ruwen7 == the derivation with e + P_tot = T d_T P_tot, raw, entry by entry"];
(* ruwen6: alpha(T, Y) from Y = d_alpha P_h / T.  Its gradients come from differentiating that definition *)
dadT = (pa/T - pTa)/paa; dadY = T/paa;
elim[e_] := e /. {Derivative[1, 0][al][tau, rr] -> dadTs D[TT[tau, rr], tau] + dadYs D[YY[tau, rr], tau],
   Derivative[0, 1][al][tau, rr] -> dadTs D[TT[tau, rr], rr] + dadYs D[YY[tau, rr], rr]};
gate["P6-chain", Simplify[{D[Derivative[0, 1][ph][tt, aa]/tt, tt], D[Derivative[0, 1][ph][tt, aa]/tt, aa]} /.
     {Derivative[1, 1][ph][tt, aa] -> pTa, Derivative[0, 2][ph][tt, aa] -> paa, Derivative[0, 1][ph][tt, aa] -> pa, tt -> T}] ===
    Simplify[{pTa/T - pa/T^2, paa/T}] &&
    Simplify[-(pTa/T - pa/T^2)/(paa/T) - dadT] === 0 && Simplify[1/(paa/T) - dadY] === 0,
  "implicit function: d_T alpha at fixed n and d_n alpha at fixed T from n = d_alpha P_h / T"];
eqs6 = Join[elim /@ divT[CYL, tmnP[TT[tau, rr], al[tau, rr]]][[1 ;; 2]], {divV[CYL, YY[tau, rr] uC]}];
vars6 = {TT[tau, rr], ur[tau, rr], YY[tau, rr]};
{As6, src6} = coefs[eqs6, vars6, {tau, rr}];
{At6, Ar6, S6} = {As6[[1]], As6[[2]], src6} //. subP /. {dadTs -> dadT, dadYs -> dadY};
r6 = juliaMatrices[FileNameJoin[{MAT, "1d_ideal_HQ_cilindrical_p_br_density.jl"}], "one_d_ideal_matrix_ruwen6"] /.
  Join[rlP, {jl["X"][3] -> Y}];
gate["P6-cols", Simplify[{At6[[All, 1 ;; 2]] - r6["At"][[All, 1 ;; 2]], Ar6[[All, 1 ;; 2]] - r6["Ax"][[All, 1 ;; 2]]}] === ConstantArray[0, {2, 3, 2}],
  "ruwen6, columns 1-2: == the derivation, raw"];
gate["P6-row3", Simplify[{At6[[3]] - r6["At"][[3]], Ar6[[3]] - r6["Ax"][[3]]}] === ConstantArray[0, {2, 3}] && Simplify[S6 - r6["source"]] === {0, 0, 0},
  "ruwen6, row 3 and the source: == the derivation, raw"];
gate["P6-c3", Simplify[{At6[[1 ;; 2, 3]] - r6["At"][[1 ;; 2, 3]], Ar6[[1 ;; 2, 3]] - r6["Ax"][[1 ;; 2, 3]]}] === ConstantArray[0, {2, 2}],
  "ruwen6, column 3 of rows 1-2: == the derivation, raw (fixed 2026-10-30: A_t and A_r were exchanged there)"];
gate["P6-ctl", Simplify[At6[[1 ;; 2, 3]] - Ar6[[1 ;; 2, 3]]] =!= {0, 0}, "control: A_t and A_r differ in that column, so the exchange would be seen"];

Print["\n", pass, " passed, ", fail, " failed  (", pass + fail, " checks evaluated)"];
If[Environment["FLUIDUM_DERIV_RENDER"] === "1",
  Get[FileNameJoin[{DERIV, "render.wl"}]]; renderUpstream[]];
Exit[If[fail == 0, 0, 1]];
