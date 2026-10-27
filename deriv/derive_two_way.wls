#!/usr/bin/env wolframscript
(* ============================================================================
   derive_two_way.wls -- one derivation of Fluidum's medium and charm systems in
   three charts, with the EOS kept symbolic, and of the fully coupled two-way
   system (charm in T^{mu nu}_total).

   Systems (all in the source-on-the-LHS convention of src/Matrix:
            A_t d_tau phi + A_x d_x phi (+ A_y d_y phi) + src = 0)
     S1  ideal hydro                       Nabla_mu T_m^{mu nu} = 0, pi = Pi = 0
     S2  viscous hydro (IS, the term set of derive_2p1d_viscous.wls):
           T_m = (e+P+Pi) uu + (P+Pi) g + pi,  e = T P'(T) - P(T)
           tauS Delta^{mn}_{ab} D pi^{ab} + pi^{mn} + 2 visc sigma^{mn} + dpp pi^{mn} theta = 0
           tauB D Pi + Pi + zeta theta + dPP Pi theta = 0
         dpp and dPP stay symbolic.  The 7-field 1-D :MIS kernel has dpp = 4/3 tauS; the
         5-field 1-D and the derived 2-D matrices take dpp as an argument (default 4/3 tauS);
         the legacy 2-D matrix has none; no medium has dPP.
     S3  charm first moment on a prescribed background (consistent closure):
           Nabla.(n u + nu) = 0,  (Ds/T) Delta^mu_l Nabla_m T_Q^{ml} + nu^mu = 0
     S4  S3 + charm second moment on a prescribed background (consistent m2):
           (Ds/2T) Delta^mu_a Delta^nu_b Nabla_l S^{lab} + dM^{mu nu} = 0,
           dM = pi_Q + piB_Q Delta,
           S = I31 (u Delta)_sym + A (u u nu)_sym + B (Delta nu)_sym + C (u dM)_sym
         (traceless channels + trace channel, as in the 1+1D derivation of the
         consistent second moment)
     S5  the two-way system: fluid rows Nabla_mu (T_m + T_Q)^{mu nu} = 0, charm rows
         as S4, medium pi/Pi rows as S2, all fields evolved together.
           T_Q = (n h - P0) uu + P0 Delta + h (u nu + nu u),  P0 = n T
         The second moment dM does not enter T_Q (c_M = 0): it is evolved, but passive.

   EOS, symbolic.  Medium pp[T] (e = T pp' - pp); charm nQ[T,alpha], hQ[T];
   second-moment coefficients i31Q[T,alpha], aQf[T], bQf[T], cQf[T] (Boltzmann:
   I31 = n T h, A = m^2 + 5 T h, B = T h, C = m K4/K3).  Transport coefficients
   visc, tauS, zeta, tauB, Ds are local symbols and are not differentiated, as in
   the hand-written matrices.  Conformal = the specialisation pp -> aSB T^4, Pi absent.

   Charts: 0+1D Bjorken (fields of tau), 1+1D cylindrical Milne (tau,r,phi,eta),
           2+1D Cartesian Milne (tau,x,y,eta).

   Checks (counted; exit 1 on any FAIL):
     G0/G1  metric compatibility; divergence check  Nabla.(X,0,0,0) = X' + X/tau
     G2     algebraic constraints of every constructed tensor
     G3     quasi-linearity (no second derivatives) of every row
     G-reduce   generic-EOS rows == the hand-written and the other generated matrices
     G-chart    2+1D (axisymmetric, on the x axis) == 1+1D; 1+1D (u^r=0) == 0+1D
     G-conserve S5 fluid rows == Nabla.(T_m+T_Q) (independent sqrt(-g) route),
                charm row == Nabla.N, the dM columns of the other rows vanish, n->0 gives S2
     G-bjorken  closed forms (ideal, IS, generic EOS)
     G-twoway   rest/homogeneous: no momentum exchange; drag: Newton's third law
     G-export   every free symbol of the generated Julia is a function argument
   The check functions are in derive_two_way_gates.wls; the export and its own checks
   are in derive_two_way_export.wls.

   Writes the generated blocks of src/Matrix/1d_viscous_HQ_cilindrical_BG.jl,
   1d_viscous_HQ_cilindrical_BG_m2.jl, 2d_viscous_HQ_BG.jl, 2d_viscous_HQ_BG_m2.jl,
   1d_two_way_cilindrical.jl and 2d_two_way.jl, only when every check passed and all three
   charts ran.

   Run:   wolframscript -file deriv/derive_two_way.wls [--charts=bj,cyl,car] [--no-export]
   Env:   TWOWAY_OUT          scratch directory (char speeds, Julia helpers, the export's own check)
                              [default $TemporaryDirectory/derive_two_way_out]
          FLUIDUM_DERIV_PROJECT  Julia project for the Julia-side checks [default the package root]
          TWOWAY_NORM_LOG     optional file that logs the size and time of every simplification
   ============================================================================ *)

$HistoryLength = 0;
SetOptions[$Output, FormatType -> OutputForm];
t0run = AbsoluteTime[];
PKG = ParentDirectory[DirectoryName[$InputFileName]];   (* the Fluidum.jl package root *)
DERIV = DirectoryName[$InputFileName];
envOr[k_, d_] := With[{v = Environment[k]}, If[StringQ[v] && v =!= "", v, d]];
juliaProject = envOr["FLUIDUM_DERIV_PROJECT", PKG];
outDir = envOr["TWOWAY_OUT", FileNameJoin[{$TemporaryDirectory, "derive_two_way_out"}]];
matDir = FileNameJoin[{PKG, "src", "Matrix"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
args = Rest[$ScriptCommandLine];
chartsArg = SelectFirst[args, StringStartsQ[#, "--charts="] &, "--charts=bj,cyl,car"];
charts = StringSplit[StringDrop[chartsArg, 9], ","];
doExport = !MemberQ[args, "--no-export"];
Print["derive_two_way.wls  charts=", charts, "  out=", outDir, "  matrices=", matDir];

pass = 0; fail = 0; gateLog = {};
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  AppendTo[gateLog, {id, TrueQ[ok], msg}];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], id, "  ", msg]);
info[s__] := Print["         ", s];
stamp[] := ToString[Round[AbsoluteTime[] - t0run]] <> " s";
sci[x_] := ToString[CForm[N[x]]];

(* ============================================================ geometry ===== *)
mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam,
    "sg" -> PowerExpand[Sqrt[-Det[g]]]|>];
divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] +
  Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];
(* the independent route: (1/sqrt(-g)) d_m (sqrt(-g) V^m) *)
divValt[geo_, V_] := Sum[D[geo["sg"] V[[a]], geo["x"][[a]]], {a, 4}]/geo["sg"];
metricCompat[geo_] := Simplify @ Table[
  D[geo["g"][[m, n2]], geo["x"][[l]]] - Sum[geo["Gam"][[s, l, m]] geo["g"][[s, n2]], {s, 4}]
  - Sum[geo["Gam"][[s, l, n2]] geo["g"][[m, s]], {s, 4}], {l, 4}, {m, 4}, {n2, 4}];
gradV[geo_, V_] := Table[D[V[[m]], geo["x"][[l]]] + Sum[geo["Gam"][[m, l, s]] V[[s]], {s, 4}], {l, 4}, {m, 4}];
DcovV[geo_, u_, W_] := Table[Sum[u[[l]] gradV[geo, W][[l, m]], {l, 4}], {m, 4}];
DcovT[geo_, u_, P_] := Table[Sum[u[[l]] (D[P[[m, n2]], geo["x"][[l]]]
   + Sum[geo["Gam"][[m, l, s]] P[[s, n2]], {s, 4}]
   + Sum[geo["Gam"][[n2, l, s]] P[[m, s]], {s, 4}]), {l, 4}], {m, 4}, {n2, 4}];
divT2[geo_, P_] := Table[Sum[D[P[[m, n2]], geo["x"][[m]]], {m, 4}]
   + Sum[geo["Gam"][[m, m, s]] P[[s, n2]], {m, 4}, {s, 4}]
   + Sum[geo["Gam"][[n2, m, s]] P[[m, s]], {m, 4}, {s, 4}], {n2, 4}];
(* independent route for a symmetric (2,0) tensor *)
divT2alt[geo_, P_] := Table[Sum[D[geo["sg"] P[[m, n2]], geo["x"][[m]]], {m, 4}]/geo["sg"]
   + Sum[geo["Gam"][[n2, m, s]] P[[m, s]], {m, 4}, {s, 4}], {n2, 4}];
divS3[geo_, S_] := Table[Sum[D[S[[l, b, c]], geo["x"][[l]]], {l, 4}] +
   Sum[geo["Gam"][[l, l, s]] S[[s, b, c]] + geo["Gam"][[b, l, s]] S[[l, s, c]] +
       geo["Gam"][[c, l, s]] S[[l, b, s]], {l, 4}, {s, 4}], {b, 4}, {c, 4}];

(* ================================================== exact zero tests ====== *)
(* Big identities are tested by exact rational arithmetic at random points
   (Schwartz-Zippel): every field value, derivative, EOS value and coefficient is a
   random rational, u is drawn from Pythagorean tuples so that u^tau is rational,
   and the residual must be exactly 0.  Small identities use Simplify. *)
pyth1 = {{3, 4, 5}, {5, 12, 13}, {8, 15, 17}, {7, 24, 25}, {20, 21, 29}, {12, 35, 37}};
pyth2 = {{1, 2, 2, 3}, {2, 3, 6, 7}, {6, 2, 3, 7}, {4, 4, 7, 9}, {7, 4, 4, 9}, {12, 3, 4, 13},
         {8, 1, 4, 9}, {9, 6, 2, 11}};
rr0[] := (2 RandomInteger[{0, 1}] - 1) RandomInteger[{1, 60}]/RandomInteger[{7, 31}];
rpos[] := RandomInteger[{5, 80}]/RandomInteger[{7, 31}];
fieldApps[e_, heads_] := Union @ Cases[e, (x : ((h_Symbol /; MemberQ[heads, h])[___])) |
     (x : (Derivative[__][h_Symbol /; MemberQ[heads, h]][___])) :> x, {0, Infinity}];
eosHeads = {pp, nQ, hQ, i31Q, aQf, bQf, cQf};
(* uheads: {{ur}} (1-D) or {{ux, uy}} (2-D) -- the velocity heads to draw Pythagorean *)
randPoint[e_, fheads_, uheads_, seed_] := Module[{apps, rules = {}, e1, eapps, syms, uv},
  SeedRandom[seed];
  apps = fieldApps[e, fheads];
  Which[
    Length[uheads] == 1,
      With[{p = RandomChoice[pyth1]}, uv = {(2 RandomInteger[{0, 1}] - 1) p[[1]]/p[[2]]}],
    Length[uheads] == 2,
      With[{p = RandomChoice[pyth2]},
        uv = {(2 RandomInteger[{0, 1}] - 1) p[[2]]/p[[1]], (2 RandomInteger[{0, 1}] - 1) p[[3]]/p[[1]]}],
    True, uv = {}];
  rules = Table[a -> If[MatchQ[a, (h_Symbol /; MemberQ[uheads, h])[___]],
        uv[[First@Flatten@Position[uheads, Head[a]]]],
        If[MatchQ[a, Derivative[__][_][___]], rr0[], rpos[]]], {a, apps}];
  e1 = e /. rules;
  eapps = fieldApps[e1, eosHeads];
  e1 = e1 /. Table[a -> If[MatchQ[a, Derivative[__][_][___]], rr0[], rpos[]], {a, eapps}];
  syms = Complement[Union@Cases[e1, _Symbol?(! NumericQ[#] && Context[#] === "Global`" &), {0, Infinity}], {}];
  e1 /. Table[s -> rpos[], {s, syms}]];
isZeroAt[e_, fheads_, uheads_, k_: 3] := AllTrue[Range[k],
  With[{v = Flatten[{Together[randPoint[e, fheads, uheads, 1000 + #]]}]},
    AllTrue[v, (# === 0 || RootReduce[#] === 0) &]] &];
zeroTest[e_, fheads_, uheads_, t_: 20] := Module[{s},
  s = TimeConstrained[Simplify[e], t, $Aborted];
  If[s === $Aborted, {isZeroAt[e, fheads, uheads, 4], "exact at 4 random rational points"},
     {AllTrue[Flatten[{s}], # === 0 &], "symbolic"}]];

(* ===================================================== charts and fields === *)
BJ  = mkGeom[{tau, xx, yy, eta}, {-1, 1, 1, tau^2}];
CYL = mkGeom[{tau, rr, ph, eta}, {-1, 1, rr^2, tau^2}];
CAR = mkGeom[{tau, xx, yy, eta}, {-1, 1, 1, tau^2}];

(* distinct heads per chart (b = Bjorken, c = cylindrical, k = Cartesian) so that the
   chart-reduction substitutions cannot alias *)
stateBJ = Module[{f, pe, qe},
  f[h_] := h[tau];
  pe = f[peeb]; qe = f[qeb];
  <|"name" -> "0p1d", "geo" -> BJ, "cs" -> {tau}, "f" -> f,
    "T" -> f[TTb], "u" -> {1, 0, 0, 0},
    "piM" -> DiagonalMatrix[{0, -pe/2, -pe/2, pe/tau^2}], "PiB" -> f[PiBb],
    "al" -> f[alb], "nu" -> {0, 0, 0, 0},
    "piQ" -> DiagonalMatrix[{0, -qe/2, -qe/2, qe/tau^2}], "PiQ" -> f[PiQb],
    "med" -> {TTb, peeb, PiBb}, "ch" -> {alb}, "m2" -> {qeb, PiQb},
    "uheads" -> {}, "bg" -> {TTb}|>];

stateCYL = Module[{f, ut, ur, prr, ptr, ptt, pl, pf, pe, lv, phv, etv},
  f[h_] := h[tau, rr];
  ur = f[urc]; ut = Sqrt[1 + ur^2];
  prr = -ut^2 (f[pffc] + f[peec]); ptr = (ur/ut) prr; ptt = (ur/ut) ptr;
  (* charm: as the medium, independent pi_Q^phi_phi, pi_Q^eta_eta (orthonormal); the radial one from the trace *)
  pf = f[pffQc]; pe = f[peeQc]; pl = -(pf + pe);
  lv = {ur, ut, 0, 0}; phv = {0, 0, 1/rr, 0}; etv = {0, 0, 0, 1/tau};
  <|"name" -> "1p1d", "geo" -> CYL, "cs" -> {tau, rr}, "f" -> f,
    "T" -> f[TTc], "u" -> {ut, ur, 0, 0},
    "piM" -> {{ptt, ptr, 0, 0}, {ptr, prr, 0, 0}, {0, 0, f[pffc]/rr^2, 0}, {0, 0, 0, f[peec]/tau^2}},
    "PiB" -> f[PiBc], "al" -> f[alc], "nu" -> {ur f[nurc]/ut, f[nurc], 0, 0},
    "piQ" -> pl Outer[Times, lv, lv] + pf Outer[Times, phv, phv] + pe Outer[Times, etv, etv],
    "PiQ" -> f[PiQc], "triad" -> {lv, phv, etv},
    "med" -> {TTc, urc, pffc, peec, PiBc}, "ch" -> {alc, nurc}, "m2" -> {pffQc, peeQc, PiQc},
    "uheads" -> {urc}, "bg" -> {TTc, urc}|>];

stateCAR = Module[{f, ux, uy, ut, pxx, ptx, pty, ptt, qxx, qxy, qyy, qtx, qty, qtt, qet},
  f[h_] := h[tau, xx, yy];
  ux = f[uxk]; uy = f[uyk]; ut = Sqrt[1 + ux^2 + uy^2];
  (* medium: as in derive_2p1d_viscous.wls, independent pi^{yy}, pi^eta_eta, pi^{xy} *)
  pxx = (2 ux uy f[pxyk] - (1 + ux^2) f[pyyk] - ut^2 f[peek])/(1 + uy^2);
  ptx = (ux pxx + uy f[pxyk])/ut; pty = (ux f[pxyk] + uy f[pyyk])/ut; ptt = (ux ptx + uy pty)/ut;
  (* charm: as the medium, independent pi_Q^{yy}, pi_Q^eta_eta, pi_Q^{xy}; pi_Q^{xx} from the trace,
     the tau row from orthogonality *)
  qyy = f[pyyQk]; qet = f[peeQk]; qxy = f[pxyQk];
  qxx = (2 ux uy qxy - (1 + ux^2) qyy - ut^2 qet)/(1 + uy^2);
  qtx = (ux qxx + uy qxy)/ut; qty = (ux qxy + uy qyy)/ut; qtt = (ux qtx + uy qty)/ut;
  <|"name" -> "2p1d", "geo" -> CAR, "cs" -> {tau, xx, yy}, "f" -> f,
    "T" -> f[TTk], "u" -> {ut, ux, uy, 0},
    "piM" -> {{ptt, ptx, pty, 0}, {ptx, pxx, f[pxyk], 0}, {pty, f[pxyk], f[pyyk], 0}, {0, 0, 0, f[peek]/tau^2}},
    "PiB" -> f[PiBk], "al" -> f[alk],
    "nu" -> {(ux f[nuxk] + uy f[nuyk])/ut, f[nuxk], f[nuyk], 0},
    "piQ" -> {{qtt, qtx, qty, 0}, {qtx, qxx, qxy, 0}, {qty, qxy, qyy, 0}, {0, 0, 0, qet/tau^2}},
    "PiQ" -> f[PiQk], "qxx" -> qxx,
    "med" -> {TTk, uxk, uyk, pyyk, peek, pxyk, PiBk}, "ch" -> {alk, nuxk, nuyk},
    "m2" -> {pyyQk, peeQk, pxyQk, PiQk},
    "uheads" -> {uxk, uyk}, "bg" -> {TTk, uxk, uyk}|>];

allHeads[st_] := Join[st["med"], st["ch"], st["m2"]];

(* ======================================================== the physics ====== *)
(* One builder, every chart.  wv is the 2-D corotation switch. *)
build[st_] := Module[{geo = st["geo"], g, gi, u, Tf, Dup, Dmix, Dlo, th, gu, Dttr, sig, Pm, em,
    Tm, Dpi, Esh, Ebu, piM, PiB, n, h, P0, eQ, nu, piQ, PiQ, dM, TQ, consM, consQ, cont, fm,
    i31, aa, bb, cq, Sten, LHS, tr, TL, pre, W = 0, NPu, aV, omg, pgw},
  g = geo["g"]; gi = geo["ginv"]; u = st["u"]; Tf = st["T"];
  Dup = gi + Outer[Times, u, u]; Dmix = Dup . g; Dlo = g . Dup . g;
  th = divV[geo, u];
  gu = gi . gradV[geo, u];                       (* Nabla^a u^b *)
  Dttr[a_, b_, c_, d_] := (1/2) (Dmix[[a, c]] Dmix[[b, d]] + Dmix[[a, d]] Dmix[[b, c]]) -
     (1/3) Dup[[a, b]] Dlo[[c, d]];
  sig = Table[Sum[Dttr[a, b, c, d] gu[[c, d]], {c, 4}, {d, 4}], {a, 4}, {b, 4}];
  piM = st["piM"]; PiB = st["PiB"];
  Pm = pp[Tf]; em = Tf pp'[Tf] - pp[Tf];
  Tm = Table[(em + Pm + PiB) u[[m]] u[[n2]] + (Pm + PiB) gi[[m, n2]] + piM[[m, n2]], {m, 4}, {n2, 4}];
  consM = divT2[geo, Tm];
  Dpi = DcovT[geo, u, piM];
  (* dpp = delta_pipi, dPP = delta_PiPi, both symbolic (see the header) *)
  Esh = Table[tauS Sum[Dttr[a, b, c, d] Dpi[[c, d]], {c, 4}, {d, 4}] + piM[[a, b]] + 2 visc sig[[a, b]] +
      dpp piM[[a, b]] th, {a, 4}, {b, 4}];
  Ebu = tauB Sum[u[[l]] D[PiB, geo["x"][[l]]], {l, 4}] + PiB + zeta th + dPP PiB th;
  (* charm *)
  n = nQ[Tf, st["al"]]; h = hQ[Tf]; P0 = n Tf; eQ = n h - P0;
  nu = st["nu"]; piQ = st["piQ"]; PiQ = st["PiQ"];
  dM = piQ + PiQ Dup;
  TQ = Table[eQ u[[a]] u[[b]] + P0 Dup[[a, b]] + h (u[[a]] nu[[b]] + nu[[a]] u[[b]]), {a, 4}, {b, 4}];
  consQ = divT2[geo, TQ];
  cont = divV[geo, n u + nu];
  fm = Table[(Ds/Tf) Sum[Dmix[[m, l]] consQ[[l]], {l, 4}] + nu[[m]], {m, 4}];
  (* second moment *)
  i31 = i31Q[Tf, st["al"]]; aa = aQf[Tf]; bb = bQf[Tf]; cq = cQf[Tf];
  Sten = Table[
    i31 (u[[a]] Dup[[b, c]] + u[[b]] Dup[[a, c]] + u[[c]] Dup[[a, b]]) +
    aa (u[[a]] u[[b]] nu[[c]] + u[[a]] u[[c]] nu[[b]] + u[[b]] u[[c]] nu[[a]]) +
    bb (Dup[[a, b]] nu[[c]] + Dup[[a, c]] nu[[b]] + Dup[[b, c]] nu[[a]]) +
    cq (u[[a]] dM[[b, c]] + u[[b]] dM[[a, c]] + u[[c]] dM[[a, b]]), {a, 4}, {b, 4}, {c, 4}];
  LHS = Dmix . divS3[geo, Sten] . Transpose[Dmix];
  tr = Sum[Dlo[[a, b]] LHS[[a, b]], {a, 4}, {b, 4}];
  TL = LHS - Dup tr/3;
  pre = Ds/(2 Tf);
  (* the 2-D corotation coupling W = (pi_Q g omega) + transpose.  It is part of the covariant
     row; wv = 0 subtracts it (see rows).  The export sets wv = 1. *)
  If[st["name"] === "2p1d",
    aV = DcovV[geo, u, u];
    NPu = Table[Sum[gi[[a, l]] gradV[geo, u][[l, n2]], {l, 4}] + u[[a]] aV[[n2]], {a, 4}, {n2, 4}];
    omg = (NPu - Transpose[NPu])/2;
    pgw = Table[Sum[piQ[[m, l]] g[[l, l]] omg[[l, n2]], {l, 4}], {m, 4}, {n2, 4}];
    W = pgw + Transpose[pgw]];
  <|"Tm" -> Tm, "TQ" -> TQ, "consM" -> consM, "consQ" -> consQ, "Esh" -> Esh . g, "Ebu" -> Ebu,
    "cont" -> cont, "fm" -> fm, "LHS" -> LHS, "tr" -> tr, "TL" -> TL, "pre" -> pre, "tauM" -> pre cq,
    "W" -> W, "Nvec" -> n u + nu, "sig" -> sig, "theta" -> th|>];

(* row selection: the row definitions of the other generators, chart by chart *)
rows[st_, B_] := Module[{nm = st["name"], E = B["Esh"], cons, ch, m2, f = st["f"], tri, lo},
  cons = B["consM"] + B["consQ"];
  Switch[nm,
   "0p1d",
     <|"medM" -> {B["consM"][[1]], E[[4, 4]]/tau^2, B["Ebu"]},
       "medQ" -> {B["consQ"][[1]], 0, 0},
       "ch" -> {B["cont"]},
       "m2" -> Module[{et = {0, 0, 0, 1/tau}, etl}, etl = st["geo"]["g"] . et;
          {B["pre"] (etl . B["LHS"] . etl - B["tr"]/3) + f[qeb], B["pre"] B["tr"]/3 + f[PiQb]}]|>,
   "1p1d",
     tri = st["triad"]; lo = (st["geo"]["g"] . #) & /@ tri;
     <|"medM" -> {B["consM"][[1]], B["consM"][[2]], E[[3, 3]]/rr^2, E[[4, 4]]/tau^2, B["Ebu"]},
       "medQ" -> {B["consQ"][[1]], B["consQ"][[2]], 0, 0, 0},
       "ch" -> {B["cont"], B["fm"][[2]]},
       "m2" -> {B["pre"] (lo[[2]] . B["LHS"] . lo[[2]] - B["tr"]/3) + f[pffQc],
                B["pre"] (lo[[3]] . B["LHS"] . lo[[3]] - B["tr"]/3) + f[peeQc],
                B["pre"] B["tr"]/3 + f[PiQc]}|>,
   "2p1d",
     <|"medM" -> {B["consM"][[1]], B["consM"][[2]], B["consM"][[3]], E[[3, 3]], E[[4, 4]]/tau^2,
                  E[[2, 3]], B["Ebu"]},
       "medQ" -> {B["consQ"][[1]], B["consQ"][[2]], B["consQ"][[3]], 0, 0, 0, 0},
       "ch" -> {B["cont"], B["fm"][[2]], B["fm"][[3]]},
       "m2" -> {B["pre"] B["TL"][[3, 3]] + f[pyyQk] - (1 - wv) B["tauM"] B["W"][[3, 3]],
                tau^2 (B["pre"] B["TL"][[4, 4]] - (1 - wv) B["tauM"] B["W"][[4, 4]]) + f[peeQk],
                B["pre"] B["TL"][[2, 3]] + f[pxyQk] - (1 - wv) B["tauM"] B["W"][[2, 3]],
                B["pre"] B["tr"]/3 + f[PiQk]}|>]];

(* quasi-linear coefficient extraction *)
coefs[eqs_, varH_, cs_] := Module[{As, src, vars = (# @@ cs) & /@ varH},
  As = Table[Table[D[eqs[[i]], D[vars[[j]], c]], {i, Length[eqs]}, {j, Length[vars]}], {c, cs}];
  (* quasi-linear: the source is the row at zero first derivatives of the evolved fields
     (check Gc-extract re-assembles A.dphi + src and compares with the raw row) *)
  src = eqs /. Flatten[Table[D[v, c] -> 0, {c, cs}, {v, vars}]];
  {As, src}];
(* The simplification rule depends on size only, so a regeneration is byte-identical: an entry of
   LeafCount <= 1e5 is Simplify'd, a larger one Together'd (four 2+1D second-moment entries).  A time
   limit would make the output depend on machine load. *)
normLog = envOr["TWOWAY_NORM_LOG", ""];
norm[e_] := If[e === 0, 0, Module[{t, r, lc = LeafCount[e]},
   {t, r} = AbsoluteTiming[If[lc > 100000, Together[e], Simplify[e]]];
   If[normLog =!= "", PutAppend[{lc, t}, normLog]];
   r]];
secondDerivs[e_, heads_] := Select[fieldApps[e, heads],
   MatchQ[#, Derivative[ks__][_][___] /; Total[{ks}] > 1] &];

(* field-off substitutions *)
offRules[heads_, nargs_] := Table[h -> Function[Evaluate[Table[Unique["z"], {nargs}]], 0], {h, heads}];

(* ======================================================== per-chart run ==== *)
results = <||>;

runChart[st_] := Module[{nm = st["name"], geo = st["geo"], cs = st["cs"], f = st["f"], B, R,
    hds = allHeads[st], nar = Length[st["cs"]], ok, msg, sys = <||>, eq, vars, A, s, med, ch, m2,
    piOff, chOff, m2Off, zt},
  Print[""]; Print["================ chart ", nm, " ================  (", stamp[], ")"];
  (* --- G0/G1 --- *)
  gate["G0-" <> nm, metricCompat[geo] === ConstantArray[0, {4, 4, 4}], "Nabla g == 0"];
  gate["G1-" <> nm, Simplify[divV[geo, {XX[tau, xx], 0, 0, 0}] - (D[XX[tau, xx], tau] + XX[tau, xx]/tau)] === 0 &&
     Simplify[divValt[geo, {XX[tau, xx], 0, 0, 0}] - (D[XX[tau, xx], tau] + XX[tau, xx]/tau)] === 0,
    "trap detector: Nabla.(X(tau,..),0,0,0) == X' + X/tau (Christoffel AND sqrt(-g) routes)"];
  (* --- G2 constraints --- *)
  With[{g = geo["g"], u = st["u"]},
    gate["G2-" <> nm, Simplify[{u . g . u + 1, (g . u) . st["piM"], Total[g st["piM"], 2],
        (g . u) . st["nu"], (g . u) . st["piQ"], Total[g st["piQ"], 2],
        st["piM"] - Transpose[st["piM"]], st["piQ"] - Transpose[st["piQ"]]}] ===
        {0, {0, 0, 0, 0}, 0, 0, {0, 0, 0, 0}, 0, ConstantArray[0, {4, 4}], ConstantArray[0, {4, 4}]},
      "u.u=-1; u.pi=0, tr pi=0 (medium); u.nu=0; u.piQ=0, tr piQ=0; both symmetric"]];
  If[nm === "1p1d",
    With[{g = geo["g"], tri = st["triad"], u = st["u"]},
      gate["G2-triad", Simplify[Flatten[{Outer[(#1 . g . #2) &, tri, tri, 1] - IdentityMatrix[3],
          (# . g . u) & /@ tri}]] === ConstantArray[0, 12], "the (l, phihat, etahat) triad is orthonormal and transverse"]]];
  Print["  building tensors ... (", stamp[], ")"];
  B = build[st];
  R = rows[st, B];
  Print["  rows built (", stamp[], ")"];
  (* A = I41/P0 enters S only through (u u nu)_sym; after the Delta Delta projection it is never
     differentiated.  Its derivative survives in the unsimplified rows and cancels; prove that its
     coefficient vanishes, then drop it, so that the rows carry no spurious aQ'. *)
  Module[{cf = D[R["m2"], Derivative[1][aQf][st["T"]]], z},
    z = zeroTest[cf, hds, st["uheads"], 20];
    gate["G-aQ-" <> nm, z[[1]] && !FreeQ[R["m2"], Derivative[1][aQf]] || FreeQ[R["m2"], Derivative[1][aQf]],
      "aQ' has identically vanishing coefficient in the second-moment rows (" <> z[[2]] <> "); dropped"]];
  R = R /. Derivative[1][aQf][_] -> 0;
  (* --- G3 quasi-linearity --- *)
  ok = Length[secondDerivs[Join[R["medM"], R["medQ"], R["ch"], R["m2"]], hds]] == 0;
  gate["G3-" <> nm, ok, "every row (medium, charm-in-T, first moment, second moment) is quasi-linear"];

  med = st["med"]; ch = st["ch"]; m2 = st["m2"];
  piOff = offRules[Rest[med] /. {urc -> Nothing, uxk -> Nothing, uyk -> Nothing}, nar];
  chOff = offRules[Rest[ch], nar];   (* nu off; alpha stays *)
  m2Off = offRules[m2, nar];

  (* --- assemble the five systems --- *)
  (* S1: ideal = the conservation rows of S2 with pi = Pi = 0 *)
  Module[{idx = Switch[nm, "0p1d", {1}, "1p1d", {1, 2}, "2p1d", {1, 2, 3}],
          vv = Switch[nm, "0p1d", {TTb}, "1p1d", {TTc, urc}, "2p1d", {TTk, uxk, uyk}]},
    sys["S1"] = <|"eqs" -> (R["medM"][[idx]] /. piOff), "vars" -> vv|>];
  sys["S2"] = <|"eqs" -> R["medM"], "vars" -> med|>;
  sys["S3"] = <|"eqs" -> (R["ch"] /. m2Off), "vars" -> ch|>;
  sys["S4"] = <|"eqs" -> Join[R["ch"], R["m2"]], "vars" -> Join[ch, m2]|>;
  sys["S5"] = <|"eqs" -> Join[R["medM"] + R["medQ"], R["ch"], R["m2"]], "vars" -> Join[med, ch, m2]|>;
  (* in 0+1D there is no transverse current: S3 is the continuity row alone *)
  gate["S3-free-" <> nm, FreeQ[sys["S3"]["eqs"], Alternatives @@ m2],
    "S3 rows carry no second-moment field"];
  gate["S34-nomedpi-" <> nm, FreeQ[Join[R["ch"], R["m2"]], Alternatives @@ (Rest[med] /. {urc -> Nothing, uxk -> Nothing, uyk -> Nothing})],
    "charm rows never see the medium pi or Pi (so S3/S4 need only T, u as background)"];

  Do[
    Print["  extracting ", k, " (", Length[sys[k]["vars"]], " fields) ... (", stamp[], ")"];
    {A, s} = coefs[sys[k]["eqs"], sys[k]["vars"], cs];
    sys[k] = Join[sys[k], <|"A" -> A, "src" -> s|>];
    , {k, {"S1", "S2", "S3", "S4", "S5"}}];
  Print["  extracted (", stamp[], ")"];

  (* --- G-conserve --- *)
  Module[{S5 = sys["S5"], nmed = Length[med], nch = Length[ch], nm2 = Length[m2], dphi, rebuilt,
      Ttot, alt, offIdx, Anon, srcnon, z, offQ, S2m, flu},
    dphi = Table[D[v @@ cs, c], {c, cs}, {v, S5["vars"]}];
    rebuilt = Sum[S5["A"][[k]] . dphi[[k]], {k, Length[cs]}] + S5["src"];
    z = zeroTest[rebuilt - S5["eqs"], hds, st["uheads"], 5];
    gate["Gc-extract-" <> nm, z[[1]], "A_t dphi + A_x dphi (+A_y dphi) + src == the raw rows (" <> z[[2]] <> ")"];
    gate["Gc-ctrl-" <> nm, !isZeroAt[rebuilt[[1]] - S5["eqs"][[1]] + 10^-9 S5["eqs"][[2]], hds, st["uheads"], 2] &&
        !isZeroAt[S5["eqs"][[1]], hds, st["uheads"], 2],
      "control: the exact random-point zero test REJECTS a 1e-9-perturbed identity and a nonzero row"];
    Ttot = B["Tm"] + B["TQ"];
    alt = divT2alt[geo, Ttot];
    flu = Switch[nm, "0p1d", {1}, "1p1d", {1, 2}, "2p1d", {1, 2, 3}];
    z = zeroTest[S5["eqs"][[flu]] - alt[[flu]], hds, st["uheads"], 20];
    gate["Gc-total-" <> nm, z[[1]], "S5 fluid rows == Nabla_mu (T_m + T_Q)^{mu nu}, sqrt(-g) route (" <> z[[2]] <> ")"];
    z = zeroTest[S5["eqs"][[nmed + 1]] - divValt[geo, B["Nvec"]], hds, st["uheads"], 20];
    gate["Gc-number-" <> nm, z[[1]], "charm row == Nabla.N, sqrt(-g) route (" <> z[[2]] <> ")"];
    (* dM is not in T_Q: the dM columns of every non-m2 row vanish, and no dM field survives in their source *)
    Anon = Table[S5["A"][[k]][[1 ;; nmed + nch, nmed + nch + 1 ;; nmed + nch + nm2]], {k, Length[cs]}];
    srcnon = S5["src"][[1 ;; nmed + nch]];
    gate["Gc-dM-" <> nm, Simplify[Anon] === ConstantArray[0, Dimensions[Anon]] && FreeQ[srcnon, Alternatives @@ m2],
      "the second moment is passive: dM columns of the fluid + first-moment rows vanish identically, and their sources carry no dM"];
    (* charm off: n -> 0 (so i31 -> 0), nu -> 0, dM -> 0 ==> S5 fluid block == S2, charm cols of fluid rows = 0 *)
    offQ = Join[{nQ -> Function[{a1, a2}, 0], i31Q -> Function[{a1, a2}, 0]}, chOff, m2Off];
    S2m = sys["S2"];
    gate["Gc-off-" <> nm,
      Simplify[(Table[S5["A"][[k]][[1 ;; nmed, 1 ;; nmed]], {k, Length[cs]}] /. offQ) - S2m["A"]] ===
          ConstantArray[0, Dimensions[S2m["A"]]] &&
      Simplify[(Table[S5["A"][[k]][[nmed + 1 ;;, 1 ;; nmed]], {k, Length[cs]}] /. offQ)] ===
          ConstantArray[0, {Length[cs], nch + nm2, nmed}] &&
      Simplify[(S5["src"][[nmed + 1 ;;]] /. offQ)] === ConstantArray[0, nch + nm2] &&
      Simplify[(S5["src"][[1 ;; nmed]] /. offQ) - S2m["src"]] === ConstantArray[0, nmed],
      "charm off (n, nu, dM -> 0): S5's medium block and medium sources == S2 exactly; the charm rows' " <>
      "medium columns and sources vanish, so charm = 0 is an invariant subspace on which S5 IS S2"];
  ];
  results[nm] = <|"st" -> st, "B" -> B, "R" -> R, "sys" -> sys|>;
  Print["  chart ", nm, " done (", stamp[], ")"];
];

If[MemberQ[charts, "bj"], runChart[stateBJ]];
If[MemberQ[charts, "cyl"], runChart[stateCYL]];
If[MemberQ[charts, "car"], runChart[stateCAR]];

Get[FileNameJoin[{DirectoryName[$InputFileName], "derive_two_way_gates.wls"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "derive_two_way_export.wls"}]];

Print[""]; Print["================ summary ================  (", stamp[], ")"];
Print["  ", pass, " passed, ", fail, " failed  (", pass + fail, " checks evaluated)"];
If[fail > 0, Print["  FAILED: ", Select[gateLog, ! #[[2]] &][[All, 1]]]];
If[Environment["FLUIDUM_DERIV_RENDER"] === "1",   (* deriv/rendered/: PDF fragment + notebook *)
  Get[FileNameJoin[{DERIV, "render.wl"}]]; renderTwoWay[]];
If[pass + fail < 20, Print["  too few checks evaluated -- treating as FAIL"]; Exit[1]];
Exit[If[fail == 0, 0, 1]];
