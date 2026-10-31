#!/usr/bin/env wolframscript
(* ============================================================================
   derive_drivers.wls -- what the drivers in src/Matrix add to the matrices: the solved form, kappa and
   tau_n, the charm EOS closures, the vacuum regulators, the cold edge and the corona update.

   Symbolic here (prefix R-).  The code itself is compared in audit_matrices.jl (prefix A-, D-), which
   this script runs: it rebuilds every system numerically with ForwardDiff, calls each driver at
   physical states and evaluates the corona update on a grid.  It shares no code with the Wolfram scripts.

   Run:    wolframscript -file deriv/derive_drivers.wls      (~3 min, most of it Julia)
   Env:    FLUIDUM_DERIV_PROJECT  Julia project with Fluidum, StaticArrays, ForwardDiff [default the package root]
   Writes: nothing into src/
   ============================================================================ *)

$HistoryLength = 0;
pass = 0; fail = 0; gateLog = {};
gate[id_, ok_, msg_] := (If[TrueQ[ok], pass++, fail++];
  AppendTo[gateLog, {id, TrueQ[ok], msg}];
  Print[If[TrueQ[ok], "  [PASS] ", "  [FAIL] "], id, "  ", msg]);
PKG = ParentDirectory[DirectoryName[$InputFileName]];
DERIV = DirectoryName[$InputFileName];
envOr[k_, d_] := With[{v = Environment[k]}, If[StringQ[v] && v =!= "", v, d]];
juliaProject = envOr["FLUIDUM_DERIV_PROJECT", PKG];

mkGeom[coords_, gdiag_] := Module[{g, ginv, Gam, n = 4},
  g = DiagonalMatrix[gdiag]; ginv = Inverse[g];
  Gam = Table[Simplify[(1/2) Sum[ginv[[a, d]] (D[g[[d, b]], coords[[c]]] +
        D[g[[d, c]], coords[[b]]] - D[g[[b, c]], coords[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  <|"x" -> coords, "g" -> g, "ginv" -> ginv, "Gam" -> Gam|>];
divV[geo_, V_] := Sum[D[V[[a]], geo["x"][[a]]], {a, 4}] + Sum[geo["Gam"][[a, a, b]] V[[b]], {a, 4}, {b, 4}];
divT[geo_, P_] := Table[Sum[D[P[[m, n2]], geo["x"][[m]]], {m, 4}]
   + Sum[geo["Gam"][[m, m, s]] P[[s, n2]], {m, 4}, {s, 4}]
   + Sum[geo["Gam"][[n2, m, s]] P[[m, s]], {m, 4}, {s, 4}], {n2, 4}];
CYL = mkGeom[{tau, rr, ph, et}, {-1, 1, rr^2, tau^2}];

(* ---------------------------------------------- kappa, tau_n: the first moment *)
Print["== kappa and tau_n from the consistent first moment =="];
utC = Sqrt[1 + ur[tau, rr]^2];
uC = {utC, ur[tau, rr], 0, 0};
nuC = {ur[tau, rr] nur[tau, rr]/utC, nur[tau, rr], 0, 0};
Tf = TT[tau, rr]; af = al[tau, rr];
n0 = nn[Tf, af]; h0 = hh[Tf]; P0 = n0 Tf;
Dup = CYL["ginv"] + Outer[Times, uC, uC]; Dmix = Dup . CYL["g"];
TQ = Table[(n0 h0 - P0) uC[[a]] uC[[b]] + P0 Dup[[a, b]] + h0 (uC[[a]] nuC[[b]] + nuC[[a]] uC[[b]]), {a, 4}, {b, 4}];
Enu = (Ds/Tf) Sum[Dmix[[2, l]] divT[CYL, TQ][[l]], {l, 4}] + nuC[[2]];
atRest = {ur[tau, rr] -> 0, Derivative[1, 0][ur][tau, rr] -> 0, Derivative[0, 1][ur][tau, rr] -> 0,
  nur[tau, rr] -> 0, Derivative[0, 1][nur][tau, rr] -> dnu};
restRow = Simplify[Enu /. atRest];
gate["R-fick", Simplify[(restRow /. Derivative[1, 0][nur][tau, rr] -> dtnu /. dnu -> 0) -
      (Ds h0/Tf dtnu + Ds Derivative[0, 1][nn][Tf, af] D[af, rr] + (Ds/Tf) (n0 + Tf Derivative[1, 0][nn][Tf, af]) D[Tf, rr])] === 0,
  "at rest the row is tau_n d_tau nu^r + kappa d_r alpha + (D_s/T)(n + T dn/dT) d_r T"];
gate["R-taun", Simplify[D[Enu, Derivative[1, 0][nur][tau, rr]]/utC - Ds h0/Tf] === 0, "the coefficient of D nu^r is tau_n = D_s h/T"];
gate["R-kappa", Simplify[(D[Enu, Derivative[0, 1][al][tau, rr]]) /. atRest /. Derivative[0, 1][nn][__] -> dna] === Ds dna, "the coefficient of d_r alpha at rest is kappa = D_s dn/dalpha"];
(* the file's tau_diffusion: the Bessel combination *)
bk[k_] := BesselK[k, z];
gate["R-bess", FullSimplify[2 bk[1] - 3 bk[3] + bk[5] - 48 bk[3]/z^2, z > 0] === 0, "2 K_1 - 3 K_3 + K_5 = 48 K_3/z^2, so the code's tau_n is D_s h/T"];
gate["R-tausp", Simplify[(1/(192 Pi^3 T0^3)) (2 Pi DT) q^2 g m^5 (48 BK3/(m/T0)^2) ex /
     (q^2 g m^2 T0 BK2 ex/(2 Pi^2)) - (DT/T0) (m BK3/BK2)/T0] === 0,
  "tau_diffusion_hadron, one species: the prefactor of the code gives (D_s/T) m K_3/K_2"];

(* ------------------------------------------------------ the charm closures *)
Print["== the closures =="];
(* I_nq = (1/(2q+1)!!) Int d^3p/((2 Pi)^3 E) E^(n-2q) p^(2q) f, f = g e^alpha e^(-E/T) *)
inq[nnq_, q_, zz_] := Module[{v},   (* T = 1, m = zz; the factor g e^alpha/(2 Pi^2) dropped *)
  v = NIntegrate[p^2/Sqrt[p^2 + zz^2] Sqrt[p^2 + zz^2]^(nnq - 2 q) p^(2 q) Exp[-Sqrt[p^2 + zz^2]],
    {p, 0, Infinity}, WorkingPrecision -> 34, PrecisionGoal -> 24, MaxRecursion -> 30];
  v/(2 q + 1)!!];
clos = Table[Module[{n1 = inq[1, 0, zz], h1 = zz BesselK[3, zz]/BesselK[2, zz], c1 = zz BesselK[4, zz]/BesselK[3, zz]},
    {n1 - zz^2 BesselK[2, zz], inq[2, 1, zz]/n1 - 1, inq[3, 1, zz]/n1 - h1, inq[4, 1, zz]/n1 - (zz^2 + 5 h1),
     inq[4, 2, zz]/n1 - h1, inq[5, 2, zz]/inq[4, 2, zz] - c1}], {zz, {3, 15/2, 19/2}}];
gate["R-clos", Max[Abs[clos]] < 10^-18, "momentum integrals at z = 3, 15/2, 19/2: n, I21 = P0, I31 = P0 h, I41 = P0 A, I42 = P0 B, I52 = I42 C"];
(* matching: delta f in the nu and the pi channel *)
(* matching in the local rest frame, u = (1,0,0,0): the angular averages of p^a p^b ... /p^n over the sphere, then the
   moments of delta f_nu = f nu_a p^a/P0 and delta f_pi = f pi_ab p^a p^b/(2 I42).  dP = d^3p/((2 Pi)^3 E). *)
sph[e_] := Integrate[(e /. {px -> Sin[th] Cos[fi], py -> Sin[th] Sin[fi], pz -> Cos[th]}) Sin[th]/(4 Pi), {th, 0, Pi}, {fi, 0, 2 Pi}];
avg = {sph[px^2], sph[px^4], sph[px^2 py^2]};
gate["R-ang", avg === {1/3, 1/5, 1/15}, "angular averages: p_x^2 -> p^2/3, p_x^4 -> p^4/5, p_x^2 p_y^2 -> p^4/15"];
(* with Int dP E^(n-2q) p^(2q) f = (2q+1)!! I_nq, so each moment below is a number times I_nq *)
mom[n_, q_] := (2 q + 1)!! Inq[n, q];
gate["R-match",
  ({nuX/P0s (1/3) mom[2, 1], nuX/P0s (1/3) mom[4, 1], nuX/P0s (1/5) mom[4, 2],
    (2 piXY)/(2 Inq[4, 2]) (1/15) mom[4, 2], (2 piXY)/(2 Inq[4, 2]) (1/15) mom[5, 2]} /. {Inq[2, 1] -> P0s}) ===
  {nuX, nuX Inq[4, 1]/P0s, 3 nuX Inq[4, 2]/P0s, piXY, piXY Inq[5, 2]/Inq[4, 2]},
  "Int dP p_x df_nu = nu^x; S^{tau tau x} = (I41/P0) nu^x; S^{xxx} = 3 (I42/P0) nu^x; Int dP p_x p_y df_pi = pi^{xy}; S^{tau xy} = (I52/I42) pi^{xy}"];
(* the T derivatives the code writes *)
hOfT = m BesselK[3, m/T0]/BesselK[2, m/T0];
cOfT = m BesselK[4, m/T0]/BesselK[3, m/T0];
gate["R-dh", FullSimplify[D[hOfT, T0] - (-(m/T0^2) m (-(BesselK[2, m/T0] + BesselK[4, m/T0])/2 BesselK[2, m/T0] -
       BesselK[3, m/T0] (-(BesselK[1, m/T0] + BesselK[3, m/T0])/2))/BesselK[2, m/T0]^2), {m > 0, T0 > 0}] === 0, "dh/dT of charm_eos"];
gate["R-dC", FullSimplify[D[cOfT, T0] - (-(m/T0^2) m (-(BesselK[3, m/T0] + BesselK[5, m/T0])/2 BesselK[3, m/T0] -
       BesselK[4, m/T0] (-(BesselK[2, m/T0] + BesselK[4, m/T0])/2))/BesselK[3, m/T0]^2), {m > 0, T0 > 0}] === 0, "dC/dT of charm_eos"];
gate["R-dn", FullSimplify[D[g m^2 T0 BesselK[2, m/T0] Exp[q al0], T0] -
     g m^2 Exp[q al0] (m BesselK[1, m/T0] + 2 T0 BesselK[2, m/T0] + m BesselK[3, m/T0])/(2 T0), {m > 0, T0 > 0}] === 0, "dn/dT of charm_eos, per species"];
gate["R-I31", {D[n1s[T0, al0] T0 h1s[T0], T0], D[n1s[T0, al0] T0 h1s[T0], al0]} === 
    {n1s[T0, al0] h1s[T0] + T0 h1s[T0] Derivative[1, 0][n1s][T0, al0] + n1s[T0, al0] T0 h1s'[T0], T0 h1s[T0] Derivative[0, 1][n1s][T0, al0]},
  "dI31/dT and dI31/dalpha of I31 = n T h"];

(* ---------------------------------------------------- the vacuum regulators *)
Print["== the regulators =="];
Nv = n0 uC + nuC;
Ealpha = divV[CYL, Nv];
regRow = Ealpha + eps (uC[[1]] D[af, tau] + uC[[2]] D[af, rr]);
gate["R-reg", Simplify[{D[regRow, D[af, tau]] - D[Ealpha, D[af, tau]], D[regRow, D[af, rr]] - D[Ealpha, D[af, rr]]}] === Simplify[eps {uC[[1]], uC[[2]]}],
  "two-way: eps u.d alpha adds eps u^tau, eps u^r to the alpha column of the number row"];
shift = Derivative[0, 1][nn][Tf, af] -> Derivative[0, 1][nn][Tf, af] + eps;
gate["R-regbg", Simplify[(Ealpha /. shift) - regRow] === 0,
  "background drivers: dn/dalpha -> dn/dalpha + eps in the number row is the same term"];
gate["R-regk", Simplify[D[Enu /. shift, eps] /. atRest] === Ds D[af, rr],
  "and through kappa = D_s dn/dalpha it adds D_s eps d alpha to the consistent current row; the two-way drivers avoid this"];
gate["R-leak", Simplify[(1/(tau rr)) (D[tau rr Nv[[1]], tau] + D[tau rr Nv[[2]], rr]) - Ealpha] === 0,
  "tau r N is conserved by the number row: the regulated row gives d_tau(tau r N^tau) + d_r(tau r N^r) = -tau r eps u.d alpha"];

(* ------------------------------------------------------------- the cold edge *)
Print["== the cold edge =="];
sm[x_] := x^3 (10 - 15 x + 6 x^2);
gate["R-step", {sm[0], sm[1], sm'[0], sm'[1], sm''[0], sm''[1]} === {0, 1, 0, 0, 0, 0}, "s(x) = x^3 (10 - 15 x + 6 x^2) is C2: s, s', s'' fixed at both ends"];
gate["R-stepw", Simplify[(1 - (fl + (1 - fl) sm[x]))/(1 - fl) - (1 - sm[x])] === 0, "f = floor + (1 - floor) s, w = 1 - s = (1 - f)/(1 - floor)"];

(* ------------------------------------------------------- the corona update *)
Print["== the corona update =="];
(* d(tau vol N^tau)/dtau with the rates as symbols, solved for the alpha rate *)
utS[uu_] := Sqrt[1 + uu^2];
Ntau[s_] := nn[T0 + s dT0, a0 + s da0] utS[u0 + s du0] + (u0 + s du0) (v0 + s dv0)/utS[u0 + s du0];
rate = D[(tau + s) vol Ntau[s], s] /. s -> 0;
sol = Solve[rate == target, da0];
aRate = da0 /. First[sol];
Nt0 = Ntau[0];
codeU = u0 du0/utS[u0];
codeNu = (du0 v0 + u0 dv0)/utS[u0] - u0 v0 codeU/utS[u0]^2;
codeRate = ((target - vol Nt0)/(tau vol) - Derivative[1, 0][nn][T0, a0] dT0 utS[u0] - nn[T0, a0] codeU - codeNu)/(utS[u0] Derivative[0, 1][nn][T0, a0]);
gate["R-cor", Length[sol] == 1 && Simplify[aRate - codeRate] === 0, "corona_alpha_rate == the solution of d(tau vol N^tau)/dtau = target, with its u^tau and nu^tau rates, symbolically"];
tele = Sum[-(Gg[i + 1] - Gg[i])/dx dx, {i, 1, 4}];
gate["R-tele", Simplify[tele + (Gg[5] - Gg[1])] === 0, "the flux differences telescope"];
llf[UL_, UR_, FL_, FR_, a_] := (FL + FR)/2 - a (UR - UL)/2;
gate["R-llf", Simplify[{llf[U1, U1, F1, F1, a], llf[U1, U2, F1, F2, a] - llf[U2, U1, F2, F1, -a]}] === {F1, 0}, "the LLF flux is consistent, and antisymmetric under exchanging the two sides with a"];
gate["R-pos", Simplify[llf[0, UR, 0, FR, a] /. FR -> a UR] === 0 && Simplify[(llf[UL, UR, FL, FR, a] /. {FL -> a UL, FR -> a UR}) - a UL] === 0,
  "with a >= |N^r|/N^tau the flux is upwind where the drift equals a: the outflow of a vacuum cell vanishes"];

(* -------------------------------------------------- the code, in Julia *)
Print["== the code: audit_matrices.jl =="];
run = RunProcess[{"julia", "--startup-file=no", "--project=" <> juliaProject, FileNameJoin[{DERIV, "audit_matrices.jl"}]}];
auditLines = Select[StringSplit[run["StandardOutput"], "\n"], StringStartsQ[#, "CHK "] &];
If[auditLines === {}, Print[run["StandardError"]]];
gate["R-audit", run["ExitCode"] === 0 && Length[auditLines] >= 50, "audit_matrices.jl ran and every check passed: " <> ToString[Length[auditLines]] <> " checks"];
auditRows = Map[Function[l, Module[{w = StringSplit[l, " ", 5]}, {w[[2]], w[[3]] === "PASS", w[[5]], w[[4]]}]], auditLines];
Scan[gate[#[[1]], #[[2]], #[[3]]] &, auditRows];

Print["\n", pass, " passed, ", fail, " failed  (", pass + fail, " checks evaluated)"];
If[Environment["FLUIDUM_DERIV_RENDER"] === "1",
  Get[FileNameJoin[{DERIV, "render.wl"}]]; renderDrivers[]];
Exit[If[fail == 0, 0, 1]];
