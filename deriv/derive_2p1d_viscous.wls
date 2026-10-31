#!/usr/bin/env wolframscript
(* ============================================================================
   derive_2p1d_viscous.wls -- Fluidum's viscous hydro, 1+1D and 2+1D, symbolically.

   Derives from the tensor equations the quasi-linear systems behind

       Matrix/1d_viscous_cilindrical.jl   one_d_viscous_matrix   (5 fields)
       Matrix/2d_viscous.jl               two_d_viscous_matrix   (7 fields)

   and writes the 2+1D one as `two_d_viscous_matrix_derived` into Matrix/2d_viscous.jl.

   Physics.  Israel-Stewart plus the DNMR term delta_pipi (symbol dpp; no other
   second-order term):

       T^{mu nu} = (e + p + Pi) u^mu u^nu + (p + Pi) g^{mu nu} + pi^{mu nu}
       (1)  Nabla_mu T^{mu nu} = 0
       (2)  tauS Delta^{mu nu}_{al be} D pi^{al be} + pi^{mu nu} + 2 eta sigma^{mu nu}
                                                    + dpp pi^{mu nu} theta = 0
       (3)  tauB D Pi + Pi + zeta theta = 0

   The dpp term is linear, so the system is  M = M_IS + dpp M_dpp , and the generated Julia
   carries the two parts separately: it returns M_IS when dpp == 0 and adds M_dpp otherwise.
   `dpp` is the routines' last argument, default `delta_pipi(tauS)` (= 4/3 tauS; the solver
   closures pass 0 for `delta_pipi = false`).  The term is a scalar times the contravariant pi^{mu nu}, lowered
   with the same `lower2` as the rest of (2).  V-G6 checks the split.

   D = u.Nabla, Delta^{mu nu} = g^{mu nu} + u^mu u^nu,
   Delta^{mu nu}_{al be} = (1/2)(Delta^mu_al Delta^nu_be + Delta^mu_be Delta^nu_al)
                           - (1/3) Delta^{mu nu} Delta_{al be},
   sigma^{mu nu} = Delta^{mu nu}_{al be} Nabla^al u^be,  theta = Nabla.u.
   e = T p'(T) - p(T) (zero chemical potential), so only p'(T) = s and p''(T) appear.
   The transport coefficients eta, tauS, zeta, tauB are local coefficients and are not
   differentiated, as in the hand-written matrices.

   Independent shear fields.  pi is symmetric, traceless and u-orthogonal, so in
   a boost-invariant 2+1D flow it has three free components.  Fluidum's choice:

       1-D cylindrical : pi^phi_phi = r^2 pi^{phi phi},  pi^eta_eta = tau^2 pi^{eta eta}
       2-D Cartesian   : pi^{yy}, pi^eta_eta = tau^2 pi^{eta eta}, pi^{xy}

   with pi^{rr} / pi^{xx} and the whole tau-row eliminated by the constraints.

   Two charts.  The cylindrical chart is the control: V-G4 requires the derivation to
   reproduce the hand-written `one_d_viscous_matrix` entry by entry.  The same construction
   then gives the Cartesian system, which replaces the legacy `two_d_viscous_matrix` (ten of
   its entries are wrong; see Matrix/2d_viscous.jl).  The medium block of the 7-field :MIS
   charm kernel `one_d_viscous_matrix_fugacity` is this system at delta_pipi = 4/3 tau_pi
   (V-G4m).

   Checks (printed with the prefix V-)
     V-G0  metric compatibility Nabla g == 0, both charts
     V-G1  divergence check: Nabla.V == X' + X/tau for V = (X(tau),0,0,0)
     V-G2  u.u == -1, u_mu pi^{mu nu} == 0, pi^mu_mu == 0, identically, both charts
     V-G3  quasi-linearity: no second derivatives of the fields survive
     V-G4  1-D: derived == `one_d_viscous_matrix`, entry by entry, symbolically, read out of
           Matrix/1d_viscous_cilindrical.jl (julia_source_reader.wl), with dpp symbolic
           (At + dpp dAt etc. of the file).  V-G4-0 / V-G4-1 repeat it at dpp = 0 and
           dpp = 4/3 tauS.
     V-G4m 1-D at dpp = 4/3 tauS: derived == the medium block (rows/cols 1-5) of the 7-field
           :MIS kernel `one_d_viscous_matrix_fugacity`, read from its .jl
     V-G5  Bjorken: on u = (1,0,0,0) the derived 2-D rows collapse to the textbook
           0+1D equations, tauS pi' + pi + 4 eta/(3 tau) + dpp pi/tau = 0 for pi^eta_eta
     V-G6  the dpp split: (IS with dpp in eq. (2)) == M_IS + dpp M_dpp, both charts

   Run:    wolframscript -file deriv/derive_2p1d_viscous.wls
   Writes: the generated block of src/Matrix/2d_viscous.jl
   ============================================================================ *)

$HistoryLength = 0;
pass = 0; fail = 0; gateLog = {};
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  AppendTo[gateLog, {"V-" <> id, TrueQ[ok], msg}];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], "V-" <> id, "  ", msg]);
PKG = ParentDirectory[DirectoryName[$InputFileName]];   (* the Fluidum.jl package root *)
DERIV = DirectoryName[$InputFileName];
Get[FileNameJoin[{DERIV, "julia_source_reader.wl"}]];

(* ---------------------------------------------------------------- geometry *)
mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam|>];

divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] +
  Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];

metricCompat[geo_] := Simplify @ Table[
  D[geo["g"][[m, n2]], geo["x"][[l]]] - Sum[geo["Gam"][[s, l, m]] geo["g"][[s, n2]], {s, 4}]
  - Sum[geo["Gam"][[s, l, n2]] geo["g"][[m, s]], {s, 4}], {l, 4}, {m, 4}, {n2, 4}];

(* Nabla_l V^m  for a vector *)
gradV[geo_, V_] := Table[D[V[[m]], geo["x"][[l]]] + Sum[geo["Gam"][[m, l, s]] V[[s]], {s, 4}], {l, 4}, {m, 4}];
(* u^l Nabla_l W^m *)
DcovV[geo_, u_, W_] := Table[Sum[u[[l]] gradV[geo, W][[l, m]], {l, 4}], {m, 4}];
(* u^l Nabla_l P^{mn} for a (2,0) tensor *)
DcovT[geo_, u_, P_] := Table[Sum[u[[l]] (D[P[[m, n2]], geo["x"][[l]]]
   + Sum[geo["Gam"][[m, l, s]] P[[s, n2]], {s, 4}]
   + Sum[geo["Gam"][[n2, l, s]] P[[m, s]], {s, 4}]), {l, 4}], {m, 4}, {n2, 4}];
(* Nabla_m P^{m n} *)
divT[geo_, P_] := Table[Sum[D[P[[m, n2]], geo["x"][[m]]], {m, 4}]
   + Sum[geo["Gam"][[m, m, s]] P[[s, n2]], {m, 4}, {s, 4}]
   + Sum[geo["Gam"][[n2, m, s]] P[[m, s]], {m, 4}, {s, 4}], {n2, 4}];

Print["== G0/G1  geometry =="];
CYL = mkGeom[{tau, rr, ph, et}, {-1, 1, rr^2, tau^2}];
CAR = mkGeom[{tau, xx, yy, et}, {-1, 1, 1, tau^2}];
gate["G0cyl", metricCompat[CYL] === ConstantArray[0, {4, 4, 4}], "cylindrical Milne: Nabla g == 0"];
gate["G0car", metricCompat[CAR] === ConstantArray[0, {4, 4, 4}], "Cartesian   Milne: Nabla g == 0"];
gate["G1cyl", Simplify[divV[CYL, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "cylindrical: Nabla.V = X' + X/tau"];
gate["G1car", Simplify[divV[CAR, {XX[tau], 0, 0, 0}] - (XX'[tau] + XX[tau]/tau)] === 0,
  "Cartesian:   Nabla.V = X' + X/tau"];

(* ------------------------------------------------------- the IS system ---- *)
(* pp[T] is the pressure; e = T p' - p.  eta, tauS, zeta, tauB are coefficients. *)
buildIS[geo_, uvec_, pimat_, PiS_, Tf_, dppA_ : 0] := Module[
  {g = geo["g"], gi = geo["ginv"], Dl, Dud, Dttr, th, gradu, sig, ee, pp0, Tmn, cons, Dpi, Esh, Ebu},
  Dl = gi + Outer[Times, uvec, uvec];                   (* Delta^{mu nu} *)
  Dud = Dl . g;                                          (* Delta^mu_nu     *)
  th = divV[geo, uvec];
  gradu = gradV[geo, uvec];                              (* Nabla_l u^m *)
  (* Nabla^al u^be = g^{al l} Nabla_l u^be *)
  With[{gu = gi . gradu},
   Dttr[a_, b_, c_, d_] := (1/2) (Dud[[a, c]] Dud[[b, d]] + Dud[[a, d]] Dud[[b, c]])
        - (1/3) Dl[[a, b]] (g . Dl . g)[[c, d]];
   sig = Table[Sum[Dttr[a, b, c, d] gu[[c, d]], {c, 4}, {d, 4}], {a, 4}, {b, 4}];
  ];
  pp0 = pp[Tf]; ee = Tf pp'[Tf] - pp[Tf];
  Tmn = Table[(ee + pp0 + PiS) uvec[[m]] uvec[[n2]] + (pp0 + PiS) gi[[m, n2]] + pimat[[m, n2]],
        {m, 4}, {n2, 4}];
  cons = divT[geo, Tmn];
  Dpi = DcovT[geo, uvec, pimat];
  Esh = Table[tauS Sum[Dttr[a, b, c, d] Dpi[[c, d]], {c, 4}, {d, 4}]
        + pimat[[a, b]] + 2 visc sig[[a, b]] + dppA pimat[[a, b]] th, {a, 4}, {b, 4}];
  Ebu = tauB Sum[uvec[[l]] D[PiS, geo["x"][[l]]], {l, 4}] + PiS + zeta th;
  <|"cons" -> cons, "Esh" -> Esh, "Ebu" -> Ebu, "sig" -> sig, "theta" -> th, "T" -> Tmn,
    "Dl" -> Dl, "Dud" -> Dud|>];

(* lower the second index: E^a_b = g_{b l} E^{a l} *)
lower2[geo_, M_] := M . geo["g"];

coefs[eqs_, vars_, cs_] := Module[{As, src},
  As = Table[Table[D[eqs[[i]], D[vars[[j]], c]], {i, Length[eqs]}, {j, Length[vars]}], {c, cs}];
  src = eqs - Sum[As[[k]] . Table[D[v, cs[[k]]], {v, vars}], {k, Length[cs]}];
  {Simplify /@ As, Simplify[src]}];

(* ===================================================== CHART 1: cylindrical *)
Print["== V-G2/G3/G4  cylindrical control (5 fields) =="];
utC = Sqrt[1 + ur[tau, rr]^2];
uC  = {utC, ur[tau, rr], 0, 0};
(* independent: pff = pi^phi_phi = rr^2 pi^{phiphi},  pee = pi^eta_eta = tau^2 pi^{etaeta} *)
prrC = -utC^2 (pff[tau, rr] + pee[tau, rr]);
ptrC = (ur[tau, rr]/utC) prrC;
pttC = (ur[tau, rr]/utC) ptrC;
piC = {{pttC, ptrC, 0, 0}, {ptrC, prrC, 0, 0}, {0, 0, pff[tau, rr]/rr^2, 0},
       {0, 0, 0, pee[tau, rr]/tau^2}};
gate["G2cyl-uu", Simplify[uC . CYL["g"] . uC + 1] === 0, "u.u == -1"];
gate["G2cyl-un", Simplify[(CYL["g"] . uC) . piC] === {0, 0, 0, 0}, "u_mu pi^{mu nu} == 0"];
gate["G2cyl-tr", Simplify[Total[CYL["g"] piC, 2]] === 0, "pi^mu_mu == 0"];

ISC = buildIS[CYL, uC, piC, PiB[tau, rr], TT[tau, rr]];
EshC = lower2[CYL, ISC["Esh"]];
eqsC = {ISC["cons"][[1]], ISC["cons"][[2]], EshC[[3, 3]]/rr^2, EshC[[4, 4]]/tau^2, ISC["Ebu"]};
varsC = {TT[tau, rr], ur[tau, rr], pff[tau, rr], pee[tau, rr], PiB[tau, rr]};
sdC = Flatten @ Table[D[eqsC, D[v, c1], D[v2, c2]], {v, varsC}, {v2, varsC}, {c1, {tau, rr}}, {c2, {tau, rr}}];
gate["G3cyl", Simplify[sdC] === ConstantArray[0, Length[sdC]], "cylindrical system is quasi-linear"];
{AsC, srcC} = coefs[eqsC, varsC, {tau, rr}];
{AtC, ArC} = AsC;
(* the delta_pipi part: the coefficient of dpp, i.e. pi^{mn} theta through the same rows *)
dppPart[geo_, uvec_, pimat_] := lower2[geo, pimat divV[geo, uvec]];
DpC = dppPart[CYL, uC, piC];
eqsDC = {0, 0, DpC[[3, 3]]/rr^2, DpC[[4, 4]]/tau^2, 0};
{AsDC, srcDC} = coefs[eqsDC, varsC, {tau, rr}];
{AtDC, ArDC} = AsDC;
ISCd = buildIS[CYL, uC, piC, PiB[tau, rr], TT[tau, rr], dpp];
EshCd = lower2[CYL, ISCd["Esh"]];
eqsCd = {ISCd["cons"][[1]], ISCd["cons"][[2]], EshCd[[3, 3]]/rr^2, EshCd[[4, 4]]/tau^2, ISCd["Ebu"]};
gate["G6cyl", Simplify[eqsCd - eqsC - dpp eqsDC] === ConstantArray[0, 5],
  "cylindrical: (IS + dpp pi theta in eq. 2) == IS + dpp * (the exported dpp part), equation by equation"];
(* V-G4: the control, against the .jl file itself *)
jlC = {jl["u"][1] -> TT[tau, rr], jl["u"][2] -> ur[tau, rr], jl["u"][3] -> pff[tau, rr],
  jl["u"][4] -> pee[tau, rr], jl["u"][5] -> PiB[tau, rr], jl["tau"] -> tau, jl["R"] -> rr,
  jl["dtp"] -> pp'[TT[tau, rr]], jl["dtdtp"] -> pp''[TT[tau, rr]], jl["p"] -> pp[TT[tau, rr]],
  jl["zeta"] -> zeta, jl["visc"] -> visc, jl["tauS"] -> tauS, jl["tauB"] -> tauB};
sh1 = juliaMatrices[FileNameJoin[{PKG, "src", "Matrix", "1d_viscous_cilindrical.jl"}],
  "one_d_viscous_matrix"] /. jlC;
(* the hand-written kernel is (At, Ax, source) + dpp (dAt, dAx, dsource); dpp stays a symbol *)
derC[d_] := {AtC + d AtDC, ArC + d ArDC, srcC + d srcDC};
shC[d_] := {sh1["At"] + d sh1["dAt"], sh1["Ax"] + d sh1["dAx"], sh1["source"] + d sh1["dsource"]};
zero5 = {ConstantArray[0, {5, 5}], ConstantArray[0, {5, 5}], ConstantArray[0, 5]};
gate["G4", Simplify[derC[dpp] - shC[dpp]] === zero5,
  "1-D: derived == shipped one_d_viscous_matrix, A_t, A_r, src entry by entry, read from the .jl, SYMBOLIC dpp  [CONTROL]"];
gate["G4-0", Simplify[derC[0] - shC[0]] === zero5, "1-D at dpp = 0 (plain IS): plain IS == the shipped IS part"];
gate["G4-1", Simplify[derC[4 tauS/3] - shC[4 tauS/3]] === zero5, "1-D at dpp = 4/3 tauS (default)"];
(* V-G4m: at the default dpp the 5-field medium is the 7-field :MIS medium block *)
sh7 = juliaMatrices[FileNameJoin[{PKG, "src", "Matrix", "1d_viscous_HQ_cilindrical_fugacity.jl"}],
  "one_d_viscous_matrix_fugacity"] /. Join[jlC, {jl["p"] -> pp'[TT[tau, rr]], jl["u"][6] -> al, jl["u"][7] -> nur,
  jl["dmp"] -> 0, jl["dtdmp"] -> 0, jl["dmdmp"] -> 0}];
blk7 = {sh7["At"][[1 ;; 5, 1 ;; 5]], sh7["Ax"][[1 ;; 5, 1 ;; 5]], sh7["source"][[1 ;; 5]]};
gate["G4m", Simplify[derC[4 tauS/3] - blk7] === zero5,
  "1-D at dpp = 4/3 tauS == the medium block of the 7-field :MIS one_d_viscous_matrix_fugacity, raw, read from the .jl"];
gate["G4m-ctl", Simplify[derC[0] - blk7] =!= zero5, "control: at dpp = 0 it is NOT the :MIS medium block"];

(* ================================================ CHART 2: Cartesian 2+1D *)
Print["== G2/G3  Cartesian 2+1D (7 fields) =="];
utK = Sqrt[1 + ux[tau, xx, yy]^2 + uy[tau, xx, yy]^2];
uK  = {utK, ux[tau, xx, yy], uy[tau, xx, yy], 0};
(* independent: pyy = pi^{yy}, pee = pi^eta_eta = tau^2 pi^{etaeta}, pxy = pi^{xy} *)
pxxK = (2 ux[tau, xx, yy] uy[tau, xx, yy] pxy[tau, xx, yy]
        - (1 + ux[tau, xx, yy]^2) pyy[tau, xx, yy]
        - utK^2 pee[tau, xx, yy])/(1 + uy[tau, xx, yy]^2);
ptxK = (ux[tau, xx, yy] pxxK + uy[tau, xx, yy] pxy[tau, xx, yy])/utK;
ptyK = (ux[tau, xx, yy] pxy[tau, xx, yy] + uy[tau, xx, yy] pyy[tau, xx, yy])/utK;
pttK = (ux[tau, xx, yy] ptxK + uy[tau, xx, yy] ptyK)/utK;
piK = {{pttK, ptxK, ptyK, 0}, {ptxK, pxxK, pxy[tau, xx, yy], 0},
       {ptyK, pxy[tau, xx, yy], pyy[tau, xx, yy], 0},
       {0, 0, 0, pee[tau, xx, yy]/tau^2}};
gate["G2car-uu", Simplify[uK . CAR["g"] . uK + 1] === 0, "u.u == -1"];
gate["G2car-un", Simplify[(CAR["g"] . uK) . piK] === {0, 0, 0, 0}, "u_mu pi^{mu nu} == 0"];
gate["G2car-tr", Simplify[Total[CAR["g"] piK, 2]] === 0, "pi^mu_mu == 0"];

ISK = buildIS[CAR, uK, piK, PiB[tau, xx, yy], TT[tau, xx, yy]];
EshK = lower2[CAR, ISK["Esh"]];
eqsK = {ISK["cons"][[1]], ISK["cons"][[2]], ISK["cons"][[3]],
        EshK[[3, 3]], EshK[[4, 4]]/tau^2, EshK[[2, 3]], ISK["Ebu"]};
varsK = {TT[tau, xx, yy], ux[tau, xx, yy], uy[tau, xx, yy], pyy[tau, xx, yy],
         pee[tau, xx, yy], pxy[tau, xx, yy], PiB[tau, xx, yy]};
sdK = Flatten @ Table[D[eqsK, D[v, c1], D[v2, c2]], {v, varsK}, {v2, varsK},
        {c1, {tau, xx, yy}}, {c2, {tau, xx, yy}}];
gate["G3car", Simplify[sdK] === ConstantArray[0, Length[sdK]], "Cartesian system is quasi-linear"];
{AsK, srcK} = coefs[eqsK, varsK, {tau, xx, yy}];
{AtK, AxK, AyK} = AsK;
DpK = dppPart[CAR, uK, piK];
eqsDK = {0, 0, 0, DpK[[3, 3]], DpK[[4, 4]]/tau^2, DpK[[2, 3]], 0};
{AsDK, srcDK} = coefs[eqsDK, varsK, {tau, xx, yy}];
{AtDK, AxDK, AyDK} = AsDK;
ISKd = buildIS[CAR, uK, piK, PiB[tau, xx, yy], TT[tau, xx, yy], dpp];
EshKd = lower2[CAR, ISKd["Esh"]];
eqsKd = {ISKd["cons"][[1]], ISKd["cons"][[2]], ISKd["cons"][[3]],
        EshKd[[3, 3]], EshKd[[4, 4]]/tau^2, EshKd[[2, 3]], ISKd["Ebu"]};
gate["G6car", Simplify[eqsKd - eqsK - dpp eqsDK] === ConstantArray[0, 7],
  "Cartesian: (IS + dpp pi theta in eq. 2) == IS + dpp * (the exported dpp part), equation by equation"];

(* -------------------------------------------------------------- G5 Bjorken *)
Print["== G5  Bjorken collapse =="];
bj = {ux[tau, xx, yy] -> 0, uy[tau, xx, yy] -> 0, pxy[tau, xx, yy] -> 0,
      Derivative[a_, b_, c_][ux][tau, xx, yy] :> 0, Derivative[a_, b_, c_][uy][tau, xx, yy] :> 0,
      Derivative[a_, b_, c_][pxy][tau, xx, yy] :> 0,
      Derivative[a_, b_, c_][TT][tau, xx, yy] :> If[b + c > 0, 0, Derivative[a, b, c][TT][tau, xx, yy]],
      Derivative[a_, b_, c_][pyy][tau, xx, yy] :> If[b + c > 0, 0, Derivative[a, b, c][pyy][tau, xx, yy]],
      Derivative[a_, b_, c_][pee][tau, xx, yy] :> If[b + c > 0, 0, Derivative[a, b, c][pee][tau, xx, yy]],
      Derivative[a_, b_, c_][PiB][tau, xx, yy] :> If[b + c > 0, 0, Derivative[a, b, c][PiB][tau, xx, yy]]};
(* on Bjorken flow: pi^{yy} = pi^{xx} = -pi^eta_eta/2 =: -pib/2 with pib = pi^eta_eta *)
bjE = Simplify[eqsKd /. bj];
(* row1: energy.  d_tau e = -(e+p+Pi+pi^eta_eta)/tau   with e = T p' - p *)
(* The enclosing parentheses are needed: without them `a\n + b` is two statements in WL and
   the second line is silently discarded. *)
row1want = (D[TT[tau, xx, yy], tau] TT[tau, xx, yy] pp''[TT[tau, xx, yy]]
   + (TT[tau, xx, yy] pp'[TT[tau, xx, yy]] + PiB[tau, xx, yy] + pee[tau, xx, yy])/tau);
res1 = Simplify[bjE[[1]] - row1want];
gate["G5-e", res1 === 0, "energy row on Bjorken == d_tau e + (e+p+Pi+pi^eta_eta)/tau"];
If[res1 =!= 0, Print["      residual = ", res1]];
(* row5: tau^2 * [ tauS d_tau pi^eta_eta + pi^eta_eta + 2 eta sigma^eta_eta ], sigma^eta_eta = -2/(3 tau) *)
(* sigma^eta_eta = +2/(3 tau) on Bjorken, so with the sign convention
   tauS D pi + pi + 2 eta sigma = 0 pi^eta_eta relaxes to -4 eta/(3 tau) < 0:
   the longitudinal pressure P_L = p + Pi + pi^eta_eta is reduced, as it should be. *)
row5want = (tauS D[pee[tau, xx, yy], tau] + pee[tau, xx, yy] + 4 visc/(3 tau) + dpp pee[tau, xx, yy]/tau)/tau^2;
res5 = Simplify[bjE[[5]] - row5want];
gate["G5-pi", res5 === 0, "pi^eta_eta row on Bjorken == tauS pi' + pi + 4 eta/(3 tau) + dpp pi/tau  (theta = 1/tau)"];
If[res5 =!= 0, Print["      residual = ", res5]];
row7want = tauB D[PiB[tau, xx, yy], tau] + PiB[tau, xx, yy] + zeta/tau;
gate["G5-Pi", Simplify[bjE[[7]] - row7want] === 0, "bulk row on Bjorken == tauB Pi' + Pi + zeta/tau"];
gate["G5-mom", Simplify[{bjE[[2]], bjE[[3]], bjE[[6]]}] === {0, 0, 0},
  "the x, y and xy rows vanish identically on Bjorken flow"];

(* The slab limit (2-D at u^y = 0, d_y = 0  ==  1-D at r -> Infinity) is not checked here:
   Limit[] on expressions with hundreds of terms costs more than the rest of the file.  It is
   a numerical comparison of the two generated routines at matched states and large R. *)

(* ------------------------------------------------------------------ export *)
Print["== export =="];
subC = {TT[tau, rr] -> uT, ur[tau, rr] -> uR, pff[tau, rr] -> pFF, pee[tau, rr] -> pEE,
   PiB[tau, rr] -> pBB, Derivative[2][pp][_] -> dtdtp, Derivative[1][pp][_] -> dtp, pp[_] -> pres,
   rr -> R};
subK = {TT[tau, xx, yy] -> uT, ux[tau, xx, yy] -> uX, uy[tau, xx, yy] -> uY,
   pyy[tau, xx, yy] -> pYY, pee[tau, xx, yy] -> pEE, pxy[tau, xx, yy] -> pXY,
   PiB[tau, xx, yy] -> pBB, Derivative[2][pp][_] -> dtdtp, Derivative[1][pp][_] -> dtp, pp[_] -> pres};
(* Simplify keeps the entries readable, but a few of the 7x7 entries are slow to simplify;
   fall back to Together, which is always correct, after 25 s. *)
norm[e_] := TimeConstrained[Simplify[e], 25, Together[e]];
jstr[e_, sub_] := StringReplace[ToString[InputForm[norm[norm[e] //. sub]]],
   {"Sqrt[" -> "sqrt(", "]" -> ")", "[" -> "("}];
mat2jl[M_, name_, n_, sub_] := "    " <> name <> " = SMatrix{" <> ToString[n] <> "," <> ToString[n] <> "}((\n" <>
   StringRiffle[Flatten @ Table["        " <> jstr[M[[i, j]], sub], {j, n}, {i, n}], ",\n"] <> "))\n";
vec2jl[V_, name_, n_, sub_] := "    " <> name <> " = SVector{" <> ToString[n] <> "}((\n" <>
   StringRiffle[Table["        " <> jstr[V[[i]], sub], {i, n}], ",\n"] <> "))\n";

hdr2 = "\"\"\"
    two_d_viscous_matrix_derived(u, tau, p, dtp, dtdtp, zeta, visc, tauS, tauB, dpp = delta_pipi(tauS))
(A_t, A_x, A_y, src) for phi = (T, u^x, u^y, pi^{yy}, pi^eta_eta, pi^{xy}, Pi), the layout of `two_d_viscous_matrix`.
Israel-Stewart + the DNMR term dpp pi theta, source on the LHS; dpp = 0 is plain Israel-Stewart.
\"\"\"
function two_d_viscous_matrix_derived(u, tau, pres, dtp, dtdtp, zeta, visc, tauS, tauB, dpp = delta_pipi(tauS))
    uT, uX, uY, pYY, pEE, pXY, pBB = u[1], u[2], u[3], u[4], u[5], u[6], u[7]
";
dppHead = "\n    # the coefficient of dpp of  dpp pi^{mn} theta\n";
body2 = mat2jl[AtK, "At", 7, subK] <> "\n" <> mat2jl[AxK, "Ax", 7, subK] <> "\n" <>
        mat2jl[AyK, "Ay", 7, subK] <> "\n" <> vec2jl[srcK, "source", 7, subK] <>
        "\n    iszero(dpp) && return (At, Ax, Ay, source)\n" <> dppHead <>
        mat2jl[AtDK, "dAt", 7, subK] <> "\n" <> mat2jl[AxDK, "dAx", 7, subK] <> "\n" <>
        mat2jl[AyDK, "dAy", 7, subK] <> "\n" <> vec2jl[srcDK, "dsource", 7, subK];
out = hdr2 <> body2 <> "\n    return (At + dpp*dAt, Ax + dpp*dAx, Ay + dpp*dAy, source + dpp*dsource)\nend\n";
If[fail == 0,
  spliceJulia[FileNameJoin[{PKG, "src", "Matrix", "2d_viscous.jl"}],
    juliaNames[out, {"uT" -> "T", "uX" -> "ux", "uY" -> "uy", "pYY" -> "piyy", "pEE" -> "pizz", "pXY" -> "pixy",
      "pBB" -> "piB"}]],
  Print["  NOT writing src/Matrix/2d_viscous.jl: a check failed"]];
Print["\n", pass, " passed, ", fail, " failed  (", pass + fail, " checks evaluated)"];
If[Environment["FLUIDUM_DERIV_RENDER"] === "1",   (* deriv/rendered/: PDF fragment + notebook *)
  Get[FileNameJoin[{DERIV, "render.wl"}]]; renderViscous[]];
Exit[If[fail == 0, 0, 1]];
