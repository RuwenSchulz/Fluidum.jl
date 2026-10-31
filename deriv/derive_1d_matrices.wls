#!/usr/bin/env wolframscript
(* ============================================================================
   derive_1d_matrices.wls -- every remaining 1+1D cylindrical matrix, against a covariant
   derivation, read out of the Julia source.

   Derives from the tensor equations, in one chart with explicit Christoffels, and compares with the
   .jl files themselves (julia_source_reader.wl).  It writes nothing into src/: these matrices are
   hand-written.

   Cylindrical Milne (tau, r, phi, eta), metric diag(-1, 1, r^2, tau^2), radial flow
   u = (sqrt(1+u_r^2), u_r, 0, 0), source on the left: A_t d_tau phi + A_r d_r phi + src = 0.

   Equations (D = u.Nabla, Delta = g + uu, sigma, theta as usual)
     medium   Nabla_m T^{mn} = 0,  T = (e+p+Pi) uu + (p+Pi) g + pi,  e = T p' - p
              tauS Delta^{mn}_{ab} D pi^{ab} + pi^{mn} + 2 eta sigma^{mn} + dpp pi^{mn} theta = 0
              tauB D Pi + Pi + zeta theta = 0
     charm    Nabla_m (n u^m + nu^m) = 0
              tau_n Delta^m_b D nu^b + nu^m + kappa Delta^{mb} d_b alpha = 0

   Checks (printed with the prefix M-)
     M-V5    5-field one_d_viscous_matrix (1d_viscous_cilindrical.jl): raw entries, symbolically, for
             symbolic dpp.  M-V5-0 / M-V5-1 state the two switch values; M-V5=V7m: at dpp = 4/3 tauS
             it equals the 7-field :MIS medium block, raw, entry by entry.
     M-BG    2-field one_d_viscous_matrix_fugacity_BG_only (1d_viscous_HQ_cilindrical_BG.jl): raw
             entries up to the Jacobian u^tau^3 r tau.
     M-V7m   the medium block (rows/cols 1-5) of the 7-field one_d_viscous_matrix_fugacity
             (1d_viscous_HQ_cilindrical_fugacity.jl): raw entries, with dpp solved for.
     M-V7c   the whole 7-field system including the charm rows 6-7, in the solved form
             A_t^-1 A_r, A_t^-1 src, exactly at random rational states.
     M-V7e   the charm rows of the file are the derived ones with d_tau T and d_tau u^r eliminated
             with the medium rows, raw: row_c -> row_c - A_t[c,med] A_t[med,med]^-1 row_med.
     M-G7/M-G7e  the gamma variant one_d_viscous_matrix_gamma (1d_viscous_HQ_cilindrical_gamma.jl):
             phi_6 = gamma, n = gamma n0(T), kappa0 multiplying d gamma.  Same two checks.
     M-I3n/M-I3a  one_d_ideal_matrix_ruwen4 (ideal medium + passive density n) and
             one_d_ideal_matrix_ruwen5 (ideal medium + passive n(T, alpha)), solved form;
             M-I3ne/M-I3ae the same elimination, raw.
   Each check also requires the derived system to be quasi-linear (no second derivatives).

   Run:    wolframscript -file deriv/derive_1d_matrices.wls      (~1 min)
   Writes: matrix_refpoints_1d.csv into FLUIDUM_REFPOINTS_DIR, when that is set: long format
           (system, point, key, value), the kernel's inputs and the derived solved form at 24 random
           states per system ("visc5" at dpp = 0, "visc5d" at dpp = 4/3 tauS), for Julia tests.
   ============================================================================ *)

$HistoryLength = 0;
pass = 0; fail = 0; gateLog = {};
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  AppendTo[gateLog, {id, TrueQ[ok], msg}];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], id, "  ", msg]);
PKG = ParentDirectory[DirectoryName[$InputFileName]];   (* the Fluidum.jl package root *)
DERIV = DirectoryName[$InputFileName];
MAT = FileNameJoin[{PKG, "src", "Matrix"}];
Get[FileNameJoin[{DERIV, "julia_source_reader.wl"}]];
SeedRandom[20261013];

(* ---------------------------------------------------------------- geometry *)
mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam|>];
divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] + Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];
gradV[geo_, V_] := Table[D[V[[m]], geo["x"][[l]]] + Sum[geo["Gam"][[m, l, s]] V[[s]], {s, 4}], {l, 4}, {m, 4}];
DcovV[geo_, u_, W_] := Table[Sum[u[[l]] gradV[geo, W][[l, m]], {l, 4}], {m, 4}];
DcovT[geo_, u_, P_] := Table[Sum[u[[l]] (D[P[[m, n2]], geo["x"][[l]]]
   + Sum[geo["Gam"][[m, l, s]] P[[s, n2]], {s, 4}] + Sum[geo["Gam"][[n2, l, s]] P[[m, s]], {s, 4}]), {l, 4}],
   {m, 4}, {n2, 4}];
divT[geo_, P_] := Table[Sum[D[P[[m, n2]], geo["x"][[m]]], {m, 4}]
   + Sum[geo["Gam"][[m, m, s]] P[[s, n2]], {m, 4}, {s, 4}]
   + Sum[geo["Gam"][[n2, m, s]] P[[m, s]], {m, 4}, {s, 4}], {n2, 4}];
metricCompat[geo_] := Simplify @ Table[
  D[geo["g"][[m, n2]], geo["x"][[l]]] - Sum[geo["Gam"][[s, l, m]] geo["g"][[s, n2]], {s, 4}]
  - Sum[geo["Gam"][[s, l, n2]] geo["g"][[m, s]], {s, 4}], {l, 4}, {m, 4}, {n2, 4}];

Print["== geometry =="];
CYL = mkGeom[{tau, rr, ph, et}, {-1, 1, rr^2, tau^2}];
gate["M-G0", metricCompat[CYL] === ConstantArray[0, {4, 4, 4}], "cylindrical Milne: Nabla g == 0"];
gate["M-G1", Simplify[divV[CYL, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "Nabla.V = X' + X/tau for V = (X(tau),0,0,0): the Christoffels are in"];

(* ------------------------------------------------------------- the systems *)
utC = Sqrt[1 + ur[tau, rr]^2];
uC  = {utC, ur[tau, rr], 0, 0};
prrC = -utC^2 (pff[tau, rr] + pee[tau, rr]); ptrC = (ur[tau, rr]/utC) prrC; pttC = (ur[tau, rr]/utC) ptrC;
piC = {{pttC, ptrC, 0, 0}, {ptrC, prrC, 0, 0}, {0, 0, pff[tau, rr]/rr^2, 0}, {0, 0, 0, pee[tau, rr]/tau^2}};
gate["M-G2", Simplify[{uC . CYL["g"] . uC + 1, Total[CYL["g"] piC, 2]}] === {0, 0} &&
             Simplify[(CYL["g"] . uC) . piC] === {0, 0, 0, 0}, "u.u = -1, pi traceless and transverse"];

(* medium rows: (energy, r-momentum, pi^phi_phi, pi^eta_eta, Pi), Fluidum's normalisation *)
medium[pimat_, PiS_, dppC_] := Module[{g = CYL["g"], gi = CYL["ginv"], Dl, Dud, th, gu, Dttr, sig, ee, Tmn, cons, Dpi, Esh},
  Dl = gi + Outer[Times, uC, uC]; Dud = Dl . g;
  th = divV[CYL, uC]; gu = gi . gradV[CYL, uC];
  Dttr[a_, b_, c_, d_] := (1/2) (Dud[[a, c]] Dud[[b, d]] + Dud[[a, d]] Dud[[b, c]]) - (1/3) Dl[[a, b]] (g . Dl . g)[[c, d]];
  sig = Table[Sum[Dttr[a, b, c, d] gu[[c, d]], {c, 4}, {d, 4}], {a, 4}, {b, 4}];
  ee = TT[tau, rr] pp'[TT[tau, rr]] - pp[TT[tau, rr]];
  Tmn = Table[(ee + pp[TT[tau, rr]] + PiS) uC[[m]] uC[[n2]] + (pp[TT[tau, rr]] + PiS) gi[[m, n2]] + pimat[[m, n2]], {m, 4}, {n2, 4}];
  cons = divT[CYL, Tmn];
  Dpi = DcovT[CYL, uC, pimat];
  Esh = (Table[tauS Sum[Dttr[a, b, c, d] Dpi[[c, d]], {c, 4}, {d, 4}] + pimat[[a, b]] + 2 visc sig[[a, b]]
          + dppC pimat[[a, b]] th, {a, 4}, {b, 4}]) . g;
  {cons[[1]], cons[[2]], Esh[[3, 3]]/rr^2, Esh[[4, 4]]/tau^2,
   tauB (D[PiS, tau] uC[[1]] + D[PiS, rr] uC[[2]]) + PiS + zeta th}];
(* charm rows: number conservation + the shipped (grad-alpha) first moment, r component *)
charm[ndens_, alf_] := Module[{nuv, Dnu, acc, nda, e2},
  nuv = {ur[tau, rr] nur[tau, rr]/utC, nur[tau, rr], 0, 0};
  Dnu = DcovV[CYL, uC, nuv]; acc = DcovV[CYL, uC, uC];
  nda = Sum[(CYL["g"] . acc)[[m]] nuv[[m]], {m, 4}];
  e2 = Table[tauDiff (Dnu[[m]] - uC[[m]] nda) + nuv[[m]] +
       kappa Sum[(CYL["ginv"] + Outer[Times, uC, uC])[[m, b]] D[alf, CYL["x"][[b]]], {b, 4}], {m, 4}];
  {divV[CYL, ndens uC + nuv], e2[[2]]}];

coefs[eqs_, vars_] := Module[{As, src},
  As = Table[Table[D[eqs[[i]], D[vars[[j]], c]], {i, Length[eqs]}, {j, Length[vars]}], {c, {tau, rr}}];
  src = eqs - Sum[As[[k]] . Table[D[v, {tau, rr}[[k]]], {v, vars}], {k, 2}];
  {Simplify[As[[1]]], Simplify[As[[2]]], Simplify[src]}];
quasiLinear[eqs_, vars_] := Simplify[Flatten @ Table[D[eqs, D[v, c1], D[v2, c2]], {v, vars}, {v2, vars},
  {c1, {tau, rr}}, {c2, {tau, rr}}]] === ConstantArray[0, Length[eqs] Length[vars]^2 4];

(* derivation symbols -> the flat names both sides are compared in *)
subD = {TT[tau, rr] -> uT, ur[tau, rr] -> uR, pff[tau, rr] -> pFF, pee[tau, rr] -> pEE, PiB[tau, rr] -> pBB,
  al[tau, rr] -> aL, nur[tau, rr] -> nuR, nf[tau, rr] -> nF, gm[tau, rr] -> gM,
  Derivative[2][pp][_] -> dtdtp, Derivative[1][pp][_] -> dtp, pp[_] -> pres,
  Derivative[1, 0][nn][__] -> dtn, Derivative[0, 1][nn][__] -> dmn, nn[__] -> n,
  Derivative[1][n0][_] -> dtn, n0[_] -> n, rr -> R};
flat[e_] := e //. subD;

(* exact solved form at a rational state: the comparison that survives row operations *)
solved[{At_, Ax_, src_}] := With[{Ai = Inverse[At]}, {Ai . Ax, Ai . src}];
(* u_r chosen so that sqrt(1 + u_r^2) is rational: u_r = (k^2-1)/(2k) *)
rstate[] := Module[{k = RandomInteger[{5, 20}]/RandomInteger[{6, 12}]},
  <|uT -> RandomInteger[{15, 50}]/100, uR -> (k^2 - 1)/(2 k), pFF -> RandomInteger[{-9, 9}]/1000,
    pEE -> RandomInteger[{-9, 9}]/1000, pBB -> RandomInteger[{-9, 9}]/2000, aL -> -RandomInteger[{1, 80}]/10,
    nuR -> RandomInteger[{-9, 9}]/10^4, nF -> RandomInteger[{1, 99}]/1000, gM -> RandomInteger[{1, 99}]/100,
    tau -> RandomInteger[{5, 50}]/10, R -> RandomInteger[{5, 80}]/10,
    pres -> RandomInteger[{5, 20}]/10, dtp -> RandomInteger[{20, 60}]/10, dtdtp -> RandomInteger[{5, 30}]/10,
    zeta -> RandomInteger[{0, 9}]/100, visc -> RandomInteger[{1, 99}]/100, tauS -> RandomInteger[{2, 12}]/10,
    tauB -> RandomInteger[{2, 12}]/10, n -> RandomInteger[{1, 99}]/1000, dtn -> RandomInteger[{1, 99}]/1000,
    dmn -> RandomInteger[{1, 99}]/1000, tauDiff -> RandomInteger[{10, 60}]/10, kappa -> RandomInteger[{1, 99}]/10^4,
    dtT -> -RandomInteger[{1, 99}]/1000, drT -> -RandomInteger[{1, 99}]/1000,
    drur -> RandomInteger[{1, 99}]/1000, dtur -> RandomInteger[{1, 99}]/1000|>];
exactSame[sysA_, sysB_, npts_] := Module[{worst = 0, st, a, b},
  Do[st = rstate[]; a = solved[sysA /. st]; b = solved[sysB /. st];
    worst = Max[worst, Max[Abs[Flatten[a - b]]]], {npts}];
  worst];

(* the charm rows of the .jl files are the derived rows with d_tau of the medium fields eliminated:
   row_c -> row_c - A_t[c,med] A_t[med,med]^-1 row_med.  Compared raw, exactly, at random states. *)
elimL[{At_, Ax_, src_}, nm_] := With[{c = nm + 1 ;; Length[At], m = 1 ;; nm}, At[[c, m]] . Inverse[At[[m, m]]]];
elimRows[sys : {At_, Ax_, src_}, nm_] := With[{c = nm + 1 ;; Length[At], m = 1 ;; nm, L = elimL[sys, nm]},
  {At[[c]] - L . At[[m]], Ax[[c]] - L . Ax[[m]], src[[c]] - L . src[[m]]}];
elimWorst[der_, sh_, nm_, npts_] := BlockRandom[Module[{worst = 0, st, a, b},
  Do[st = rstate[]; a = elimRows[der /. st, nm];
    b = {sh[[1, nm + 1 ;;]], sh[[2, nm + 1 ;;]], sh[[3, nm + 1 ;;]]} /. st;
    worst = Max[worst, Max[Abs[Flatten[a - b]]]], {npts}];
  worst]];

(* ---- export collector: long format (system, point, key, value) ---- *)
fmt[x_] := ToString[NumberForm[N[x, 25], 18, NumberFormat -> (If[#3 == "", #1, #1 <> "e" <> #3] &),
  ExponentFunction -> (If[-3 <= # <= 3, Null, #] &)]];
rows = {};
exportSys[name_, sys_, inputs_, npts_] := Do[Module[{st = rstate[], sv},
    sv = solved[sys /. st];
    Do[AppendTo[rows, {name, p, "in:" <> ToString[k], fmt[st[k]]}], {k, inputs}];
    Do[AppendTo[rows, {name, p, "M:" <> ToString[i] <> ":" <> ToString[j], fmt[sv[[1, i, j]]]}],
      {i, Length[sv[[2]]]}, {j, Length[sv[[2]]]}];
    Do[AppendTo[rows, {name, p, "S:" <> ToString[i], fmt[sv[[2, i]]]}], {i, Length[sv[[2]]]}]],
  {p, npts}];

(* ========================================================= M-V5 (control) *)
Print["== M-V5  5-field one_d_viscous_matrix (control) =="];
vars5 = {TT[tau, rr], ur[tau, rr], pff[tau, rr], pee[tau, rr], PiB[tau, rr]};
eqs5 = medium[piC, PiB[tau, rr], 0];
gate["M-V5q", quasiLinear[eqs5, vars5], "quasi-linear"];
der5 = flat /@ coefs[eqs5, vars5];
jl5 = {jl["u"][1] -> uT, jl["u"][2] -> uR, jl["u"][3] -> pFF, jl["u"][4] -> pEE, jl["u"][5] -> pBB,
  jl["tau"] -> tau, jl["R"] -> R, jl["p"] -> pres, jl["dtp"] -> dtp, jl["dtdtp"] -> dtdtp, jl["zeta"] -> zeta,
  jl["visc"] -> visc, jl["tauS"] -> tauS, jl["tauB"] -> tauB};
sh5a = juliaMatrices[FileNameJoin[{MAT, "1d_viscous_cilindrical.jl"}], "one_d_viscous_matrix"] /. jl5;
(* the shipped kernel as a function of dpp: IS part + dpp * (dAt, dAx, dsource), dpp kept a symbol *)
sh5d[d_] := {sh5a["At"] + d sh5a["dAt"], sh5a["Ax"] + d sh5a["dAx"], sh5a["source"] + d sh5a["dsource"]};
sh5 = sh5d[0];
eqs5d = medium[piC, PiB[tau, rr], dpp];
der5d = flat /@ coefs[eqs5d, vars5];
zero5 = {ConstantArray[0, {5, 5}], ConstantArray[0, {5, 5}], ConstantArray[0, 5]};
gate["M-V5", Simplify[der5d - sh5d[dpp]] === zero5,
  "derived IS + dpp pi theta == the shipped 5-field A_t, A_r, src, RAW entry by entry, read from the .jl, SYMBOLIC dpp"];
gate["M-V5-0", Simplify[der5 - sh5d[0]] === zero5, "FLUIDUM_DELTA_PIPI=0: dpp = 0, plain IS == the shipped IS part"];
gate["M-V5-1", Simplify[(der5d /. dpp -> 4 tauS/3) - sh5d[4 tauS/3]] === zero5,
  "FLUIDUM_DELTA_PIPI=1 (default): dpp = 4/3 tauS"];
exportSys["visc5", der5, {uT, uR, pFF, pEE, pBB, tau, R, pres, dtp, dtdtp, zeta, visc, tauS, tauB}, 24];

(* ======================================================= M-BG (control) *)
Print["== M-BG  2-field one_d_viscous_matrix_fugacity_BG_only (control) =="];
bgR = {Derivative[1, 0][ur][tau, rr] -> dtur, Derivative[0, 1][ur][tau, rr] -> drur,
       Derivative[1, 0][TT][tau, rr] -> dtT, Derivative[0, 1][TT][tau, rr] -> drT};
varsBG = {al[tau, rr], nur[tau, rr]};
eqsBG = charm[nn[TT[tau, rr], al[tau, rr]], al[tau, rr]];
gate["M-BGq", quasiLinear[eqsBG, varsBG], "quasi-linear"];
derBG = flat /@ ((coefs[eqsBG, varsBG]) /. bgR);
shBG = Values[juliaMatrices[FileNameJoin[{MAT, "1d_viscous_HQ_cilindrical_BG.jl"}], "one_d_viscous_matrix_fugacity_BG_only"]] /.
  {jl["u"][1] -> aL, jl["u"][2] -> nuR, jl["tau"] -> tau, jl["r"] -> R, jl["ur"] -> uR, jl["dtT"] -> dtT,
   jl["drT"] -> drT, jl["drur"] -> drur, jl["dtur"] -> dtur, jl["n"] -> n, jl["dnUdmu"] -> dmn,
   jl["dnUdT"] -> dtn, jl["tauDiffusion"] -> tauDiff, jl["kappaG"] -> kappa};
jacBG = (1 + uR^2)^(3/2) R tau;
dBG = Simplify[jacBG derBG - shBG];
gate["M-BG", dBG === {ConstantArray[0, {2, 2}], ConstantArray[0, {2, 2}], {0, 0}},
  "derived first moment * u^tau^3 r tau == 1d_viscous_HQ_cilindrical_BG.jl, raw, read from the .jl"];
exportSys["bg2", derBG, {aL, nuR, tau, R, uR, dtT, drT, drur, dtur, n, dmn, dtn, tauDiff, kappa}, 24];

(* ================================================== M-V7m: medium block *)
Print["== M-V7m  medium block of the 7-field one_d_viscous_matrix_fugacity =="];
vars7 = Join[vars5, {al[tau, rr], nur[tau, rr]}];
sh7 = Values[juliaMatrices[FileNameJoin[{MAT, "1d_viscous_HQ_cilindrical_fugacity.jl"}], "one_d_viscous_matrix_fugacity"]];
jl7 = Join[jl5, {jl["u"][6] -> aL, jl["u"][7] -> nuR, jl["n"] -> n, jl["dtn"] -> dtn, jl["dmn"] -> dmn,
  jl["tauDiff"] -> tauDiff, jl["kappa"] -> kappa, jl["dmp"] -> 0, jl["dtdmp"] -> 0, jl["dmdmp"] -> 0}];
sh7 = sh7 /. jl7;
blk[s_] := {s[[1, 1 ;; 5, 1 ;; 5]], s[[2, 1 ;; 5, 1 ;; 5]], s[[3, 1 ;; 5]]};
dppSol = Solve[Simplify[der5d[[1, 3, 2]] - sh7[[1, 3, 2]]] == 0, dpp];
Print["      dpp solved from A_t[3,2]: ", dppSol];
gate["M-V7m-dpp", dppSol === {{dpp -> (4 tauS)/3}}, "the medium carries dpp = 4/3 tauS (DNMR delta_pipi)"];
d7m = Simplify[(der5d /. dpp -> 4 tauS/3) - blk[sh7]];
gate["M-V7m", d7m === {ConstantArray[0, {5, 5}], ConstantArray[0, {5, 5}], ConstantArray[0, 5]},
  "medium rows/cols 1-5 == IS + (4/3) tauS pi theta, RAW entry by entry"];
gate["M-V7m-ctl", Simplify[(der5d /. dpp -> 0) - blk[sh7]] =!= {ConstantArray[0, {5, 5}], ConstantArray[0, {5, 5}], ConstantArray[0, 5]},
  "control: plain IS (dpp = 0) does NOT reproduce the 7-field medium block"];
Print["      medium-block entries that differ from the 5-field one_d_viscous_matrix at dpp = 0: ",
  Count[Flatten[Simplify[blk[sh7] - sh5]], Except[0]]];
gate["M-V5=V7m", Simplify[sh5d[4 tauS/3] - blk[sh7]] === zero5,
  "at the default dpp = 4/3 tauS the 5-field kernel IS the 7-field :MIS medium block, RAW entry by entry (both read from the .jl)"];

(* ===================================================== M-V7c: all rows *)
Print["== M-V7c  7-field system incl. charm rows (solved form, exact) =="];
eqs7 = Join[medium[piC, PiB[tau, rr], 4 tauS/3], charm[nn[TT[tau, rr], al[tau, rr]], al[tau, rr]]];
gate["M-V7q", quasiLinear[eqs7, vars7], "quasi-linear"];
der7 = flat /@ coefs[eqs7, vars7];
w7 = exactSame[der7, sh7, 8];
gate["M-V7c", w7 === 0, "A_t^-1 A_r and A_t^-1 src == derived, EXACTLY (rational arithmetic), 8 random states"];
gate["M-V7e", elimWorst[der7, sh7, 5, 8] === 0,
  "charm rows 6-7 of the .jl == derived rows with d_tau T, d_tau u^r eliminated with rows 1-5, raw, exact, 8 random states"];
gate["M-V7e-ctl", BlockRandom[With[{st = rstate[]}, (der7[[2, 6 ;;]] /. st) =!= (sh7[[2, 6 ;;]] /. st)]],
  "control: without the elimination the charm rows of A_r differ"];
exportSys["visc7", der7, {uT, uR, pFF, pEE, pBB, aL, nuR, tau, R, pres, dtp, dtdtp, zeta, visc, tauS, tauB,
  n, dtn, dmn, tauDiff, kappa}, 24];

(* ======================================================= M-G7: gamma *)
Print["== M-G7  gamma variant one_d_viscous_matrix_gamma (solved form, exact) =="];
varsG = Join[vars5, {gm[tau, rr], nur[tau, rr]}];
eqsG = Join[medium[piC, PiB[tau, rr], 4 tauS/3], charm[gm[tau, rr] n0[TT[tau, rr]], gm[tau, rr]]];
derG = flat /@ coefs[eqsG, varsG];
shG = Values[juliaMatrices[FileNameJoin[{MAT, "1d_viscous_HQ_cilindrical_gamma.jl"}], "one_d_viscous_matrix_gamma"]] /.
  Join[{jl["u"][6] -> gM}, jl7];     (* first: jl7 maps u[6] to alpha *)
wG = exactSame[derG, shG, 8];
gate["M-G7", wG === 0, "n = gamma n0(T), kappa0 d gamma: solved form == derived EXACTLY, 8 random states"];
If[wG =!= 0, Print["      worst |difference| = ", N[wG]]];
gate["M-G7e", elimWorst[derG, shG, 5, 8] === 0,
  "charm rows 6-7 of the .jl == derived rows with d_tau T, d_tau u^r eliminated, raw, exact, 8 random states"];
exportSys["gamma7", derG, {uT, uR, pFF, pEE, pBB, gM, nuR, tau, R, pres, dtp, dtdtp, zeta, visc, tauS, tauB,
  n, dtn, dmn, tauDiff, kappa}, 24];

(* ======================================================= M-I3: ideal *)
Print["== M-I3  ideal medium + passive charm, 3 fields (solved form, exact) =="];
zero4 = ConstantArray[0, {4, 4}];
eqsIn = Join[medium[zero4, 0, 0][[1 ;; 2]], {divV[CYL, nf[tau, rr] uC]}];
varsIn = {TT[tau, rr], ur[tau, rr], nf[tau, rr]};
derIn = flat /@ coefs[eqsIn, varsIn];
shIn = Values[juliaMatrices[FileNameJoin[{MAT, "1d_ideal_HQ_cilindrical_density.jl"}], "one_d_ideal_matrix_ruwen4"]] /.
  {jl["X"][1] -> uT, jl["X"][2] -> uR, jl["X"][3] -> nF, jl["tau"] -> tau, jl["r"] -> R,
   jl["dPUdT"] -> dtp, jl["dPUdTdT"] -> dtdtp};
wIn = exactSame[derIn, shIn, 8];
gate["M-I3n", wIn === 0, "ruwen4 (T, u^r, n): solved form == derived EXACTLY"];
If[wIn =!= 0, Print["      worst |difference| = ", N[wIn]]];
gate["M-I3ne", elimWorst[derIn, shIn, 2, 8] === 0,
  "charm row 3 of the .jl == derived row with d_tau T, d_tau u^r eliminated, raw, exact, 8 random states"];
exportSys["ideal_n", derIn, {uT, uR, nF, tau, R, dtp, dtdtp}, 24];

eqsIa = Join[medium[zero4, 0, 0][[1 ;; 2]], {divV[CYL, nn[TT[tau, rr], al[tau, rr]] uC]}];
varsIa = {TT[tau, rr], ur[tau, rr], al[tau, rr]};
derIa = flat /@ coefs[eqsIa, varsIa];
shIa = Values[juliaMatrices[FileNameJoin[{MAT, "1d_ideal_HQ_cilindrical_fugacity.jl"}], "one_d_ideal_matrix_ruwen5"]] /.
  {jl["X"][1] -> uT, jl["X"][2] -> uR, jl["X"][3] -> aL, jl["tau"] -> tau, jl["r"] -> R,
   jl["dPUdT"] -> dtp, jl["dPUdTdT"] -> dtdtp, jl["n"] -> n, jl["dnUdT"] -> dtn, jl["dnUdalpha"] -> dmn};
wIa = exactSame[derIa, shIa, 8];
gate["M-I3a", wIa === 0, "ruwen5 (T, u^r, alpha): solved form == derived EXACTLY"];
If[wIa =!= 0, Print["      worst |difference| = ", N[wIa]]];
gate["M-I3ae", elimWorst[derIa, shIa, 2, 8] === 0,
  "charm row 3 of the .jl == derived row with d_tau T, d_tau u^r eliminated, raw, exact, 8 random states"];
exportSys["ideal_a", derIa, {uT, uR, aL, tau, R, dtp, dtdtp, n, dtn, dmn}, 24];

(* the 5-field kernel at the default dpp = 4/3 tauS (exported LAST so the states of every system
   above are the ones written before P2b) *)
exportSys["visc5d", der5d /. dpp -> 4 tauS/3, {uT, uR, pFF, pEE, pBB, tau, R, pres, dtp, dtdtp, zeta, visc, tauS, tauB}, 24];

(* the CSV carries the derived solved form, and is written only if every check passed *)
refDir = Environment["FLUIDUM_REFPOINTS_DIR"];
If[StringQ[refDir] && refDir =!= "",
  outcsv = FileNameJoin[{refDir, "matrix_refpoints_1d.csv"}];
  If[fail == 0,
    Export[outcsv, Prepend[rows, {"system", "point", "key", "value"}], "CSV"];
    Print["\n  wrote ", outcsv, "  (", Length[rows], " rows)"],
    Print["\n  not writing ", outcsv, ": a check failed"]]];
Print["\n", pass, " passed, ", fail, " failed  (", pass + fail, " checks evaluated)"];
If[Environment["FLUIDUM_DERIV_RENDER"] === "1",   (* deriv/rendered/: PDF fragment + notebook *)
  Get[FileNameJoin[{DERIV, "render.wl"}]]; renderOneD[]];
Exit[If[fail == 0, 0, 1]];
