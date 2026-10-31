(* ============================================================================
   render.wl -- the human-readable form of a derivation: a LaTeX fragment for
   deriv/rendered/derivations.tex and a notebook deriv/rendered/<name>.nb.

   Loaded by the derive_*.wls scripts only when FLUIDUM_DERIV_RENDER=1, after the
   derivation and its checks have run.  It reads the script's own objects (matrices,
   sources, field vectors, the check log) and writes nothing into src/.

   Every entry in a Result section is an expression the script computed, renamed with
   the script's own export rules, simplified for reading, and then compared with the
   unsimplified expression at random exact points.  An entry broken over several lines
   is re-assembled from its pieces and compared again.  A Steps section shows the script's
   intermediate objects the same way, and checks at random exact points that they assemble
   into the rows.  Only the Setup equations are
   written here, transcribed from the scripts' header comments.  Problems (an unnamed
   symbol, a mismatch) are printed with the prefix [render] and listed at the end;
   they never change the script's own check counts.
   ============================================================================ *)

rdOutDir = FileNameJoin[{DERIV, "rendered"}];
If[!DirectoryQ[rdOutDir], CreateDirectory[rdOutDir, CreateIntermediateDirectories -> True]];
rdProblems = {};
rdNote[s__] := (AppendTo[rdProblems, StringJoin[ToString /@ {s}]]; Print["  [render] ", s]);
rdStats = <|"entries" -> 0, "verified" -> 0, "broken" -> 0|>;
rdW = 70;   (* line width for a displayed entry, in the units of rdBW (about one letter) *)

(* ------------------------------------------------------------------------------
   1. LaTeX (a small math subset) -> notebook boxes.  Used for the Setup equations,
      for the symbol names inside result boxes, and to estimate typeset widths.
   ------------------------------------------------------------------------------ *)
rdGreek = <|
  "alpha" -> "\[Alpha]", "beta" -> "\[Beta]", "gamma" -> "\[Gamma]", "delta" -> "\[Delta]",
  "epsilon" -> "\[Epsilon]", "zeta" -> "\[Zeta]", "eta" -> "\[Eta]", "theta" -> "\[Theta]",
  "kappa" -> "\[Kappa]", "lambda" -> "\[Lambda]", "mu" -> "\[Mu]", "nu" -> "\[Nu]", "xi" -> "\[Xi]",
  "pi" -> "\[Pi]", "rho" -> "\[Rho]", "sigma" -> "\[Sigma]", "tau" -> "\[Tau]", "phi" -> "\[Phi]",
  "varphi" -> "\[CurlyPhi]", "chi" -> "\[Chi]", "psi" -> "\[Psi]", "omega" -> "\[Omega]",
  "Gamma" -> "\[CapitalGamma]", "Delta" -> "\[CapitalDelta]", "Theta" -> "\[CapitalTheta]",
  "Lambda" -> "\[CapitalLambda]", "Pi" -> "\[CapitalPi]", "Sigma" -> "\[CapitalSigma]",
  "Phi" -> "\[CapitalPhi]", "Omega" -> "\[CapitalOmega]", "ell" -> "\[ScriptL]",
  "nabla" -> "\[Del]", "partial" -> "\[PartialD]", "cdot" -> "\[CenterDot]", "times" -> "\[Times]",
  "to" -> "\[RightArrow]", "equiv" -> "\[Congruent]", "le" -> "\[LessEqual]", "ge" -> "\[GreaterEqual]",
  "ne" -> "\[NotEqual]", "pm" -> "\[PlusMinus]", "infty" -> "\[Infinity]", "langle" -> "\[LeftAngleBracket]",
  "rangle" -> "\[RightAngleBracket]", "ldots" -> "\[Ellipsis]", "dots" -> "\[Ellipsis]",
  "cdots" -> "\[CenterEllipsis]", "prime" -> "\[Prime]", "perp" -> "\[UpTee]", "bullet" -> "\[Bullet]",
  "Leftrightarrow" -> "\[DoubleLeftRightArrow]", "Rightarrow" -> "\[DoubleRightArrow]",
  "approx" -> "\[TildeTilde]", "in" -> "\[Element]", "sum" -> "\[Sum]", "int" -> "\[Integral]", "max" -> "max", "min" -> "min"|>;
rdScript[l_String] := ToExpression["\"\\[ScriptCapital" <> l <> "]\""];

rdTB[s_String] := Module[{toks, pos = 1, n, peek, next, row, rawArg, arg, seq, atom, ws},
  toks = StringCases[s, RegularExpression["\\\\[a-zA-Z]+|\\\\.|[{}^_]|\\s+|[0-9]+|."]];
  n = Length[toks];
  peek[] := If[pos <= n, toks[[pos]], None];
  next[] := (pos++; toks[[pos - 1]]);
  ws[] := While[pos <= n && StringMatchQ[toks[[pos]], WhitespaceCharacter ..], pos++];
  row[l_] := Which[l === {}, "", Length[l] == 1, First[l], True, RowBox[l]];
  rawArg[] := Module[{d, out = "", t},
    ws[];
    If[peek[] =!= "{", Return[next[]]];
    next[]; d = 1;
    While[pos <= n, t = next[];
      If[t === "{", d++]; If[t === "}", d--; If[d == 0, Break[]]];
      out = out <> t];
    out];
  arg[] := (ws[]; If[peek[] === "{", next[]; seq["}"], atom[]]);
  seq[stop_] := Module[{out = {}, a, sup, sub, t},
    While[pos <= n && peek[] =!= stop,
      t = peek[];
      If[StringMatchQ[t, WhitespaceCharacter ..], next[]; Continue[]];
      a = If[t === "^" || t === "_", "", atom[]];
      sup = None; sub = None;
      While[peek[] === "^" || peek[] === "_", If[next[] === "^", sup = arg[], sub = arg[]]];
      AppendTo[out, Which[sup =!= None && sub =!= None, SubsuperscriptBox[a, sub, sup],
        sup =!= None, SuperscriptBox[a, sup], sub =!= None, SubscriptBox[a, sub], True, a]]];
    If[stop =!= None && peek[] === stop, next[]];
    row[out]];
  atom[] := Module[{t = next[], c},
    Which[
      t === "{", seq["}"],
      t === "'", "\[Prime]",
      t === "&", "",
      StringStartsQ[t, "\\"],
        c = StringDrop[t, 1];
        Switch[c,
          "frac" | "tfrac" | "dfrac", FractionBox[arg[], arg[]],
          "sqrt", SqrtBox[arg[]],
          "text" | "mathrm" | "textrm" | "operatorname", StyleBox[rawArg[], FontSlant -> "Plain"],
          "texttt", StyleBox[rawArg[], FontFamily -> "Courier"],
          "mathcal", rdScript[rawArg[]],
          "hat", OverscriptBox[arg[], "^"],
          "dot", OverscriptBox[arg[], "."],
          "bar" | "overline", OverscriptBox[arg[], "_"],
          "tilde", OverscriptBox[arg[], "~"],
          "boldsymbol" | "mathbf", StyleBox[arg[], FontWeight -> Bold],
          "left" | "right" | "bigl" | "bigr" | "Bigl" | "Bigr" | "big" | "Big" | "biggl" | "biggr",
            With[{d = next[]}, Switch[d, ".", "", "\\{", "{", "\\}", "}", "\\langle", "\[LeftAngleBracket]",
              "\\rangle", "\[RightAngleBracket]", _, d]],
          "," | ";" | ":" | " " | "quad" | "qquad", " ",
          "!", "",
          "{", "{", "}", "}",
          "\\", "",
          _, If[KeyExistsQ[rdGreek, c], rdGreek[c], rdNote["TeX->boxes: unknown command ", t]; c]],
      True, t]];
  seq[None]];

(* typeset width of a box, roughly in letters: a fraction is as wide as its wider part *)
rdBW[b_] := Switch[Head[b],
  String, Which[b === "", 0, MemberQ[{"+", "-", "="}, b], 2.6, StringMatchQ[b, WhitespaceCharacter ..], 0.5,
    True, StringLength[b]],
  RowBox, Total[rdBW /@ First[b]],
  FractionBox, Max[rdBW[b[[1]]], rdBW[b[[2]]]] + 0.4,
  SuperscriptBox | SubscriptBox, rdBW[b[[1]]] + 0.7 rdBW[b[[2]]],
  SubsuperscriptBox, rdBW[b[[1]]] + 0.7 Max[rdBW[b[[2]]], rdBW[b[[3]]]],
  SqrtBox, rdBW[b[[1]]] + 1.2,
  StyleBox | OverscriptBox, rdBW[b[[1]]],
  _, 1];
rdWidth[tex_String] := rdBW[rdTB[tex]];

(* ------------------------------------------------------------------------------
   2. A computed expression -> LaTeX and -> boxes.  Every Global symbol is mapped
      through rdSymTeX (symbol name -> LaTeX), set by the caller from the script's own
      export names.  A symbol without a name is reported, never guessed.
   ------------------------------------------------------------------------------ *)
rdSymTeX = <||>;
rdLetters[i_] := StringJoin[FromCharacterCode /@ (97 + IntegerDigits[i, 26, 3])];
rdSyms[e_] := Union @ Cases[e, s_Symbol /; Context[s] === "Global`", {0, Infinity}, Heads -> False];
rdAtomicQ[t_String] := StringMatchQ[t, RegularExpression["[A-Za-z0-9]|\\\\[a-zA-Z]+|\\\\mathcal\\{[A-Z]\\}|\\(.*\\)"]];
rdTeXOf[sym_Symbol] := With[{nm = SymbolName[sym]},
  If[KeyExistsQ[rdSymTeX, nm], rdSymTeX[nm],
    rdNote["no LaTeX name for the symbol ", nm]; "\\mathrm{" <> StringReplace[nm, "_" -> "\\_"] <> "}"]];

rdTeX[e_] := Module[{syms = rdSyms[e], phs, s, reps},
  phs = Table["zqz" <> rdLetters[i], {i, Length[syms]}];
  s = ToString[TeXForm[e /. Thread[syms -> (Symbol /@ phs)]]];
  reps = Flatten @ Table[With[{t = rdTeXOf[syms[[i]]]},
      {"\\text{" <> phs[[i]] <> "}^" -> If[rdAtomicQ[t], "{" <> t <> "}^", "\\left(" <> t <> "\\right)^"],
       "\\text{" <> phs[[i]] <> "}" -> "{" <> t <> "}"}], {i, Length[syms]}];
  s = StringReplace[s, reps];
  If[StringContainsQ[s, "zqz"], rdNote["placeholder left in LaTeX: ", s]];
  s];

rdBoxes[e_] := Module[{syms = rdSyms[e], phs, b, map, atom},
  phs = Table["zqz" <> rdLetters[i], {i, Length[syms]}];
  b = ToBoxes[e /. Thread[syms -> (Symbol /@ phs)], TraditionalForm];
  If[MatchQ[b, FormBox[_, TraditionalForm]], b = First[b]];
  map = Association @ Table[phs[[i]] -> rdTB[rdTeXOf[syms[[i]]]], {i, Length[syms]}];
  atom = Association @ Table[phs[[i]] -> rdAtomicQ[rdTeXOf[syms[[i]]]], {i, Length[syms]}];
  b /. {SuperscriptBox[p_String /; KeyExistsQ[map, p], x_] :>
          SuperscriptBox[If[atom[p], map[p], RowBox[{"(", map[p], ")"}]], x],
        p_String /; KeyExistsQ[map, p] :> map[p]}];

(* equal at random exact points: every symbol a random rational, the velocities drawn from
   Pythagorean tuples so that u^tau is rational too; the residual must be exactly 0 *)
rdPy1 = {{3, 4}, {5, 12}, {8, 15}, {7, 24}, {20, 21}};          (* 1 + (a/b)^2 is a square *)
rdPy2 = {{1, 2, 2}, {2, 3, 6}, {6, 2, 3}, {4, 4, 7}, {12, 3, 4}}; (* 1 + (b/a)^2 + (c/a)^2 is a square *)
rdVerify[a_, b_, what_] := Module[{syms = Union[rdSyms[a], rdSyms[b]], ok = True, rule, v, p1, p2},
  rdStats["verified"]++;
  Do[SeedRandom[4242 + k];
    p1 = RandomChoice[rdPy1]; p2 = RandomChoice[rdPy2];
    rule = Table[s -> Switch[SymbolName[s],
        "uR" | "ur", p1[[1]]/p1[[2]], "uX" | "ux", p2[[2]]/p2[[1]], "uY" | "uy", -p2[[3]]/p2[[1]],
        _, RandomInteger[{1, 97}]/RandomInteger[{3, 41}]], {s, syms}];
    v = Together[(a - b) /. rule];
    If[v =!= 0, ok = False], {k, 3}];
  If[!ok, rdNote["MISMATCH (", what, "): displayed form != derived expression"]];
  ok];

(* the script's expression -> the displayed one: rename (the script's rules), simplify, and
   write sqrt(1 + u_perp^2) as u^tau; then check that nothing changed *)
rdNestQ[e_] := AnyTrue[Cases[e, _Times, {0, Infinity}, Heads -> False],
  With[{a = Numerator[#], b = Denominator[#]}, b =!= 1 && (!FreeQ[a, Power[_, _?Negative]] || !FreeQ[b, Power[_, _?Negative]])] &];
rdFlatten[d_] := If[rdNestQ[d], With[{f = TimeConstrained[Together[d], 5, d]}, If[rdNestQ[f], d, f]], d];
rdPrep[e_, rules_, ut_, what_ : ""] := Module[{o = e //. rules, d},
  d = TimeConstrained[Simplify[o], 10, o];
  d = rdFlatten[d /. ut[[1]]];
  rdVerify[d /. ut[[2]], o, what];
  d];
rdNoUt = {{}, {}};

(* ------------------------------------------------------------------------------
   3. Breaking a long entry over lines.  Returns {lines, the expression re-assembled
      from the pieces}; the caller checks the re-assembly against the entry.
   ------------------------------------------------------------------------------ *)
rdBreak[e_, w_] := Module[{s = rdTeX[e]},
  If[rdWidth[s] <= w, Return[{{s}, e}]];
  Which[
    Head[e] === Plus, rdBreakSum[List @@ e, w],
    Head[e] === Times || MatchQ[e, Power[_, _?Negative]], rdBreakProd[e, w],
    True, {{s}, e}]];

rdTermLines[t_, w_] := Module[{neg, tp, ls, r},
  neg = MatchQ[t, _?Negative] || (Head[t] === Times && MatchQ[First[t], _?Negative]);
  tp = If[neg, -t, t];
  {ls, r} = rdBreak[tp, w - 3];
  ls[[1]] = If[neg, "- ", "+ "] <> ls[[1]];
  {ls, If[neg, -r, r]}];

rdBreakSum[terms_, w_] := Module[{lines = {}, cur = "", curW = 0, rebuilt = 0, tl, tr, tw},
  Do[
    {tl, tr} = rdTermLines[t, w];
    rebuilt += tr;
    If[Length[tl] == 1,
      tw = rdWidth[tl[[1]]];
      If[cur =!= "" && curW + tw > w, AppendTo[lines, cur]; cur = ""; curW = 0];
      cur = If[cur === "", tl[[1]], cur <> " " <> tl[[1]]]; curW += tw,
      If[cur =!= "", AppendTo[lines, cur]];
      lines = Join[lines, Most[tl]]; cur = Last[tl]; curW = rdWidth[cur]],
    {t, terms}];
  If[cur =!= "", AppendTo[lines, cur]];
  lines[[1]] = StringReplace[lines[[1]], StartOfString ~~ "+ " -> ""];
  {lines, rebuilt}];

rdBreakProd[e_, w_] := Module[{num = Numerator[e], den = Denominator[e], fs, c, sgn = 1, cand, big,
    rest, pre, post = "", inner, il, ir, lines},
  fs = If[Head[num] === Times, List @@ num, {num}];
  c = Times @@ Select[fs, NumericQ]; fs = Select[fs, !NumericQ[#] &];
  If[c < 0, sgn = -1; c = -c];
  den = den Denominator[c]; c = Numerator[c];
  cand = Select[fs, Head[#] === Plus || MatchQ[#, Power[_Plus, _Integer?Positive]] &];
  If[cand === {}, Return[{{rdTeX[e]}, e}]];
  big = First @ MaximalBy[cand, StringLength[rdTeX[#]] &];
  rest = Times @@ Delete[fs, First @ FirstPosition[fs, big, {1}, Heads -> False]];
  pre = If[sgn < 0, "-", ""] <>
    If[den =!= 1, "\\frac{" <> rdTeX[c] <> "}{" <> rdTeX[den] <> "}", If[c =!= 1, rdTeX[c], ""]] <>
    If[rest =!= 1, rdTeX[rest], ""];
  If[Head[big] === Power, inner = big[[1]]; post = "^{" <> ToString[big[[2]]] <> "}", inner = big];
  {il, ir} = rdBreak[inner, Max[w - rdWidth[pre] - 2, 40]];
  lines = il;
  lines[[1]] = pre <> "\\Bigl(" <> lines[[1]];
  lines[[-1]] = lines[[-1]] <> "\\Bigr)" <> post;
  {lines, sgn c rest If[Head[big] === Power, ir^big[[2]], ir]/den}];

(* ------------------------------------------------------------------------------
   4. Items -> LaTeX.  An item is one of
        {"sec", t} {"sub", t} {"subsub", t}   headings (t may hold $math$)
        {"text", t}                           a line of prose with $math$
        {"eq", {tex, ...}}                    displayed equations (Setup)
        {"entries", {{lhsTeX, expr}, ...}}    computed entries, one per row
        {"sparsity", {{labelTeX, boolMatrix}, ...}}
        {"checks", log, skipped}              log = {{id, ok, msg}, ...}, skipped = {{id, why}, ...}
        {"source", {file, ...}}               notebook only: the script text
   ------------------------------------------------------------------------------ *)
rdEsc[s_String] := StringReplace[s, {"\\" -> "\\textbackslash{}", "{" -> "\\{", "}" -> "\\}",
  "_" -> "\\_\\allowbreak{}", "^" -> "\\textasciicircum{}", "#" -> "\\#", "%" -> "\\%", "&" -> "\\&", "$" -> "\\$",
  "~" -> "\\~{}", "<" -> "\\textless{}", ">" -> "\\textgreater{}", "/" -> "/\\allowbreak{}",
  "=" -> "=\\allowbreak{}", "," -> ",\\allowbreak{}"}];
rdEscId[s_String] := StringReplace[s, {"_" -> "\\_\\allowbreak{}", "-" -> "-\\allowbreak{}", "#" -> "\\#",
  "&" -> "\\&", "%" -> "\\%", "$" -> "\\$"}];

(* a check's message, typeset: its ASCII symbols -> math, the rest escaped *)
rdMsgRules = {
  "trap detector: " -> "",
  "Nabla_mu (T_m + T_Q + dM)^{mu nu}" -> "\\nabla_{\\mu}(T_{m} + T_{Q} + \\delta M)^{\\mu\\nu}", "Nabla.(T_Q + dM)^i = -(T/Ds) nu^i" -> "\\nabla_{\\mu}(T_{Q} + \\delta M)^{\\mu i} = -(T/D_{s})\\nu^{i}",
  "(T/Ds) nu^i" -> "(T/D_{s})\\nu^{i}", "Nabla_m dM^{mr} = d_r (pi_rr + Pi_Q) + (pi_rr - pi^phi_phi)/r" -> "\\nabla_{\\mu}\\delta M^{\\mu r} = \\partial_{r}(\\pi_{rr} + \\Pi_{Q}) + (\\pi_{rr} - \\pi^{\\phi}{}_{\\phi})/r",
  "d_tau(e + e_Q) + (e+P+Pi+pi + n h + piQ^eta_eta + Pi_Q)/tau" -> "\\partial_{\\tau}(e + e_{Q}) + (e + P + \\Pi + \\pi + nh + \\pi_{Q}{}^{\\eta}{}_{\\eta} + \\Pi_{Q})/\\tau",
  "tau_pipi, lambda_1" -> "\\tau_{\\pi\\pi}, \\lambda_{1}", "e + P = T d_T P_tot + alpha d_alpha P_tot" -> "e + P = T\\partial_{T}P_{\\mathrm{tot}} + \\alpha\\partial_{\\alpha}P_{\\mathrm{tot}}", "e + P = T d_T P_tot" -> "e + P = T\\partial_{T}P_{\\mathrm{tot}}", "tau_n Delta D nu + nu + Ds n Delta grad(alpha) = 0" -> "\\tau_{n}\\Delta D\\nu + \\nu + D_{s}n\\,\\Delta\\nabla\\alpha = 0",
  "d(tau r N^tau)/dtau == -(G_{i+1/2} - G_{i-1/2})/dr" -> "\\partial_{\\tau}(\\tau rN^{\\tau}) = -(G_{i+1/2} - G_{i-1/2})/\\Delta r",
  "d(tau N^tau)/dtau" -> "\\partial_{\\tau}(\\tau N^{\\tau})", "d(tau vol N^tau)/dtau" -> "\\partial_{\\tau}(\\tau VN^{\\tau})",
  "d_tau(tau r N^tau) + d_r(tau r N^r) = -tau r eps u.d alpha" -> "\\partial_{\\tau}(\\tau rN^{\\tau}) + \\partial_{r}(\\tau rN^{r}) = -\\tau r\\epsilon\\,u\\cdot\\partial\\alpha",
  "tau_n d_tau nu^r + kappa d_r alpha + (D_s/T)(n + T dn/dT) d_r T" -> "\\tau_{n}\\partial_{\\tau}\\nu^{r} + \\kappa\\partial_{r}\\alpha + (D_{s}/T)(n + T\\partial_{T}n)\\partial_{r}T",
  "tau_n = D_s h/T" -> "\\tau_{n} = D_{s}h/T", "kappa = D_s dn/dalpha" -> "\\kappa = D_{s}\\partial_{\\alpha}n", "kappa = Ds n" -> "\\kappa = D_{s}n",
  "tau_diffusion_hadron == D_s h/T" -> "\\tau_{n}\\text{ of the code} = D_{s}h/T",
  "diffusion_hadron == D_s dn/dalpha" -> "\\kappa\\text{ of the code} = D_{s}\\partial_{\\alpha}n",
  "D_s -> f(T) D_s" -> "D_{s}\\to f(T)D_{s}", "dn/dalpha -> dn/dalpha + eps" -> "\\partial_{\\alpha}n\\to\\partial_{\\alpha}n + \\epsilon",
  "2 K_1 - 3 K_3 + K_5 = 48 K_3/z^2" -> "2K_{1} - 3K_{3} + K_{5} = 48K_{3}/z^{2}", "(D_s/T) m K_3/K_2" -> "(D_{s}/T)mK_{3}/K_{2}",
  "D_s h/T" -> "D_{s}h/T", "D_s eps d alpha" -> "D_{s}\\epsilon\\,\\partial\\alpha", "D nu^r" -> "D\\nu^{r}", "d_r alpha" -> "\\partial_{r}\\alpha",
  "eps u^tau, eps u^r" -> "\\epsilon u^{\\tau}, \\epsilon u^{r}", "eps u.d alpha" -> "\\epsilon\\,u\\cdot\\partial\\alpha", "u.d alpha" -> "u\\cdot\\partial\\alpha",
  "dmn_eps = 1e-3" -> "\\epsilon = 10^{-3}", "dmn_eps = 0" -> "\\epsilon = 0",
  "e + P_tot = T d_T P_tot + alpha d_alpha P_tot" -> "e + P_{\\mathrm{tot}} = T\\partial_{T}P_{\\mathrm{tot}} + \\alpha\\partial_{\\alpha}P_{\\mathrm{tot}}",
  "e + P_tot = T d_T P_tot" -> "e + P_{\\mathrm{tot}} = T\\partial_{T}P_{\\mathrm{tot}}", "n = d_alpha P_h / T" -> "n = \\partial_{\\alpha}P_{h}/T",
  "d_T alpha" -> "\\partial_{T}\\alpha", "d_n alpha" -> "\\partial_{n}\\alpha", "dn/dalpha" -> "\\partial_{\\alpha}n", "dn/dT" -> "\\partial_{T}n", "dh/dT" -> "\\partial_{T}h",
  "dC/dT" -> "\\partial_{T}\\mathcal{C}", "dI31/dT" -> "\\partial_{T}I_{31}", "dI31/dalpha" -> "\\partial_{\\alpha}I_{31}", "I31 = n T h" -> "I_{31} = nTh",
  "s(x) = x^3 (10 - 15 x + 6 x^2)" -> "s(x) = x^{3}(10 - 15x + 6x^{2})", "s, s', s''" -> "s, s', s''", "f = floor + (1 - floor) s, w = 1 - s = (1 - f)/(1 - floor)" -> "f = f_{0} + (1 - f_{0})s,\\ w = 1 - s = (1 - f)/(1 - f_{0})",
  "corona_weight = (1 - f)/(1 - floor)" -> "w = (1 - f)/(1 - f_{0})",   "T_edge - width" -> "T_{\\mathrm{edge}} - w_{T}", "T_edge" -> "T_{\\mathrm{edge}}", "f = floor" -> "f = f_{0}", "du_alpha = (1 - w) quasi-linear + w flux form" -> "\\dot\\alpha = (1 - w)\\dot\\alpha_{\\mathrm{ql}} + w\\dot\\alpha_{\\mathrm{flux}}",
  "delta f_nu = f nu.p/P0, delta f_pi = f pi_ab p^a p^b/(2 I42)" -> "\\delta f_{\\nu} = f\\nu\\cdot p/P_{0},\\ \\delta f_{\\pi} = f\\pi_{ab}p^{a}p^{b}/(2I_{42})",
  "Int dP p_x df_nu = nu^x; S^{tau tau x} = (I41/P0) nu^x; S^{xxx} = 3 (I42/P0) nu^x; Int dP p_x p_y df_pi = pi^{xy}; S^{tau xy} = (I52/I42) pi^{xy}" ->
    "\\int dP\\,p_{x}\\delta f_{\\nu} = \\nu^{x},\\ S^{\\tau\\tau x} = (I_{41}/P_{0})\\nu^{x},\\ S^{xxx} = 3(I_{42}/P_{0})\\nu^{x},\\ \\int dP\\,p_{x}p_{y}\\delta f_{\\pi} = \\pi^{xy},\\ S^{\\tau xy} = (I_{52}/I_{42})\\pi^{xy}",
  "p_x^2 -> p^2/3, p_x^4 -> p^4/5, p_x^2 p_y^2 -> p^4/15" -> "p_{x}^{2}\\to p^{2}/3,\\ p_{x}^{4}\\to p^{4}/5,\\ p_{x}^{2}p_{y}^{2}\\to p^{4}/15",
  "I31/n, I41/P0, I42/P0, I52/I42" -> "I_{31}/n,\\ I_{41}/P_{0},\\ I_{42}/P_{0},\\ I_{52}/I_{42}", "T h, A, B, C" -> "Th,\\ \\mathcal{A},\\ \\mathcal{B},\\ \\mathcal{C}", "B = T h" -> "\\mathcal{B} = Th",
  "n, I21 = P0, I31 = P0 h, I41 = P0 A, I42 = P0 B, I52 = I42 C" -> "n,\\ I_{21} = P_{0},\\ I_{31} = P_{0}h,\\ I_{41} = P_{0}\\mathcal{A},\\ I_{42} = P_{0}\\mathcal{B},\\ I_{52} = I_{42}\\mathcal{C}",
  "z = 3, 15/2, 19/2" -> "z = 3,\\ 15/2,\\ 19/2", "A_t of the file is the derived A_r and A_r of the file is the derived A_t" -> "A_{t}$ of the file is the derived $A_{r}$ and $A_{r}$ of the file is the derived $A_{t}",
  "eta, tau_pi, zeta, tau_Pi" -> "\\eta, \\tau_{\\pi}, \\zeta, \\tau_{\\Pi}", "kappa, tau_n" -> "\\kappa, \\tau_{n}", "tau_n" -> "\\tau_{n}", "n, kappa" -> "n, \\kappa", "alpha = 0" -> "\\alpha = 0", "kappa" -> "\\kappa", "D_s" -> "D_{s}", "pi^{yy} and pi^{xy}" -> "\\pi^{yy}$ and $\\pi^{xy}",
  "tau r N" -> "\\tau rN", "Ds" -> "D_{s}", " READ FROM THE FILE" -> " read from the file", "RAW" -> "raw", "SHIPPED " -> "",
  "EXACTLY" -> "exactly", "SYMBOLIC" -> "symbolic", "  [CONTROL]" -> " (control)", "REJECTS" -> "rejects",
  "NOT" -> "not", " IS " -> " is ",
  "Nabla_mu (T_m + T_Q)^{mu nu}" -> "\\nabla_{\\mu}(T_{m} + T_{Q})^{\\mu\\nu}",
  "Nabla.(X(tau,..),0,0,0) == X' + X/tau" -> "\\nabla_{\\mu}(X, 0, 0, 0)^{\\mu} = \\partial_{\\tau}X + X/\\tau",
  "Nabla.V = X' + X/tau" -> "\\nabla_{\\mu}V^{\\mu} = \\partial_{\\tau}X + X/\\tau",
  "V = (X(tau),0,0,0)" -> "V^{\\mu} = (X(\\tau), 0, 0, 0)", "V=(X(tau),0,0,0)" -> "V^{\\mu} = (X(\\tau), 0, 0, 0)",
  "Nabla.T_Q^i" -> "\\nabla_{\\mu}T_{Q}^{\\mu i}", "Nabla.T_m^i" -> "\\nabla_{\\mu}T_{m}^{\\mu i}",
  "Delta.Nabla.T_Q" -> "\\Delta^{\\mu}{}_{\\nu}\\nabla_{\\lambda}T_{Q}^{\\lambda\\nu}",
  "Nabla.T_m + Nabla.T_Q" -> "\\nabla_{\\mu}T_{m}^{\\mu\\nu} + \\nabla_{\\mu}T_{Q}^{\\mu\\nu}",
  "Nabla.N" -> "\\nabla_{\\mu}N^{\\mu}", "Nabla g" -> "\\nabla_{\\lambda}g_{\\mu\\nu}",
  "sqrt(-g)" -> "\\sqrt{-g}", "u^tau^3 r tau" -> "(u^{\\tau})^{3} r\\tau", "u^tau" -> "u^{\\tau}",
  "A_t dphi + A_x dphi (+A_y dphi) + src" -> "A_{t}\\partial_{\\tau}\\phi + A_{x}\\partial_{x}\\phi + A_{y}\\partial_{y}\\phi + S",
  "A_t, A_r, src" -> "A_{t}, A_{r}, S", "A_t, A_x, src" -> "A_{t}, A_{x}, S", "A_t^-1 A_r" -> "A_{t}^{-1}A_{r}",
  "A_t^-1 src" -> "A_{t}^{-1}S", "A_t[3,2]" -> "(A_{t})_{32}", "A_t" -> "A_{t}", "A_r" -> "A_{r}", "A_x" -> "A_{x}",
  " src " -> " S ", "src_live" -> "S^{L}", "src_background" -> "S", "src  " -> "S ",
  "u_mu pi^{mu nu}" -> "u_{\\mu}\\pi^{\\mu\\nu}", "pi^mu_mu" -> "\\pi^{\\mu}{}_{\\mu}", "pi^{mu nu}" -> "\\pi^{\\mu\\nu}",
  "u.u=-1" -> "u\\cdot u = -1", "u.u = -1" -> "u\\cdot u = -1", "u.u" -> "u\\cdot u", "u.nu=0" -> "u\\cdot\\nu = 0",
  "u.nu" -> "u\\cdot\\nu", "u.piQ=0" -> "u_{\\mu}\\pi_{Q}^{\\mu\\nu} = 0", "tr piQ=0" -> "\\pi_{Q}{}^{\\mu}{}_{\\mu} = 0",
  "u.pi=0" -> "u_{\\mu}\\pi^{\\mu\\nu} = 0", "tr pi=0" -> "\\pi^{\\mu}{}_{\\mu} = 0", "piQ" -> "\\pi_{Q}",
  "pi^eta_eta" -> "\\pi^{\\eta}{}_{\\eta}", "E^{phiphi}" -> "E^{\\phi\\phi}", "E^{etaeta}" -> "E^{\\eta\\eta}",
  "E^{yy}" -> "E^{yy}", "E^{xy}" -> "E^{xy}", "^{yy}" -> "^{yy}", "^{xy}" -> "^{xy}", "(phi phi)" -> "(\\phi\\phi)",
  "(eta eta)" -> "(\\eta\\eta)", "^eta_eta" -> "{}^{\\eta}{}_{\\eta}", "nu^r" -> "\\nu^{r}", "nu^x" -> "\\nu^{x}", "nu^y" -> "\\nu^{y}",
  "u^r = 0" -> "u^{r} = 0", "u^r" -> "u^{r}", "u^y -> 0" -> "u^{y}\\to 0", "d_y -> 0" -> "\\partial_{y}\\to 0",
  "(alpha, nu^x, nu^y)" -> "(\\alpha, \\nu^{x}, \\nu^{y})", "(n, nu, dM -> 0)" -> "(n, \\nu, \\delta M\\to 0)",
  "nu = 0, dM = 0, grad alpha = 0" -> "\\nu = 0, \\delta M = 0, \\nabla\\alpha = 0", "dM" -> "\\delta M",
  "d_tau(e + e_Q) + (e+P+Pi+pi + n h)/tau" -> "\\partial_{\\tau}(e + e_{Q}) + (e + P + \\Pi + \\pi + nh)/\\tau",
  "d_tau e + (e+p+Pi+pi^eta_eta)/tau" -> "\\partial_{\\tau}e + (e + p + \\Pi + \\pi^{\\eta}{}_{\\eta})/\\tau",
  "tauS pi' + pi + 4 eta/(3 tau) + dpp pi/tau" -> "\\tau_{\\pi}\\partial_{\\tau}\\pi + \\pi + 4\\eta/(3\\tau) + \\delta_{\\pi\\pi}\\pi/\\tau",
  "(theta = 1/tau)" -> "(\\theta = 1/\\tau)", "tauB Pi' + Pi + zeta/tau" -> "\\tau_{\\Pi}\\partial_{\\tau}\\Pi + \\Pi + \\zeta/\\tau",
  "dT/dtau == -T/(3 tau)" -> "\\partial_{\\tau}T = -T/(3\\tau)",
  "dT/dtau == -P'/(P'' tau) == -c_s^2 T/tau, c_s^2 = P'/(T P'')" -> "\\partial_{\\tau}T = -P'/(P''\\tau) = -c_{s}^{2}T/\\tau,\\ c_{s}^{2} = P'/(TP'')",
  "T' = -(T P' + Pi + pi)/(tau T P''), tauS pi' = -pi - 4 visc/(3 tau) - dpp pi/tau, tauB Pi' = -Pi - zeta/tau - dPP Pi/tau" ->
    "\\partial_{\\tau}T = -(TP' + \\Pi + \\pi)/(\\tau TP''),\\ \\tau_{\\pi}\\partial_{\\tau}\\pi = -\\pi - 4\\eta/(3\\tau) - \\delta_{\\pi\\pi}\\pi/\\tau,\\ \\tau_{\\Pi}\\partial_{\\tau}\\Pi = -\\Pi - \\zeta/\\tau - \\delta_{\\Pi\\Pi}\\Pi/\\tau",
  "dT/dtau = -T/(3tau) - pi/(12 a T^3 tau), tauS dpi/dtau = -pi - 4 visc/(3 tau) - dpp pi/tau" ->
    "\\partial_{\\tau}T = -T/(3\\tau) - \\pi/(12aT^{3}\\tau),\\ \\tau_{\\pi}\\partial_{\\tau}\\pi = -\\pi - 4\\eta/(3\\tau) - \\delta_{\\pi\\pi}\\pi/\\tau",
  "pi = -4 visc/(3 tau): dT/dtau = -T/(3tau) + (4/9)(eta/s)/tau^2, s = 4 a T^3" ->
    "\\pi = -4\\eta/(3\\tau)$: $\\partial_{\\tau}T = -T/(3\\tau) + (4/9)(\\eta/s)/\\tau^{2},\\ s = 4aT^{3}",
  "(1/tau) d_tau (n tau) = 0" -> "\\tau^{-1}\\partial_{\\tau}(n\\tau) = 0",
  "Nabla.T_Q^i = -(T/Ds) nu^i" -> "\\nabla_{\\mu}T_{Q}^{\\mu i} = -(T/D_{s})\\nu^{i}",
  "Nabla.T_m^i = +(T/Ds) nu^i" -> "\\nabla_{\\mu}T_{m}^{\\mu i} = (T/D_{s})\\nu^{i}",
  "dpp = 4/3 tauS" -> "\\delta_{\\pi\\pi} = \\tfrac{4}{3}\\tau_{\\pi}", "dpp = dPP = 0" -> "\\delta_{\\pi\\pi} = \\delta_{\\Pi\\Pi} = 0",
  "dPP = 0" -> "\\delta_{\\Pi\\Pi} = 0", "dpp = 0" -> "\\delta_{\\pi\\pi} = 0", "dpp pi theta" -> "\\delta_{\\pi\\pi}\\pi\\theta",
  "(4/3) tauS pi theta" -> "\\tfrac{4}{3}\\tau_{\\pi}\\pi\\theta", "dpp *" -> "\\delta_{\\pi\\pi}\\times", "dpp" -> "\\delta_{\\pi\\pi}",
  "kappa0 d gamma" -> "\\kappa_{0}\\partial\\gamma", "n = gamma n0(T)" -> "n = \\gamma n_{0}(T)", "(T, u^r, n)" -> "(T, u^{r}, n)",
  "(T, u^r, alpha)" -> "(T, u^{r}, \\alpha)", "d_tau T, d_tau u^r" -> "\\partial_{\\tau}T, \\partial_{\\tau}u^{r}", "(l, phihat, etahat)" -> "(\\ell, \\hat\\phi, \\hat\\eta)",
  "delta_pipi" -> "\\delta_{\\pi\\pi}", " * " -> "\\times ", "aQ'" -> "\\mathcal{A}'", " == " -> " = ", "pi or Pi" -> "\\pi$ or $\\Pi", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}",
  "T_Q" -> "T_{Q}", "T_m" -> "T_{m}", "(alpha, " -> "(\\alpha, ", "alpha" -> "\\alpha", " AND " -> " and ", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}", "T_Q" -> "T_{Q}", "T_m" -> "T_{m}"};
rdMsg[s_String] := Module[{parts = StringSplit[StringReplace[s, (First[#] -> If[StringFreeQ[Last[#], "\\"] && !StringMatchQ[Last[#], ___ ~~ ("^" | "_") ~~ ___], Last[#], "\[Section]" <> Last[#] <> "\[Section]"]) & /@ rdMsgRules], "\[Section]", All]},
  StringReplace[StringJoin[MapIndexed[If[OddQ[First[#2]], rdEsc[#1], "$" <> #1 <> "$"] &, parts]], "$$" -> ""]];

rdEntryTeX[{lhs_, d_}] := Module[{lines, rb},
  rdStats["entries"]++;
  {lines, rb} = rdBreak[d, rdW];
  If[Length[lines] > 1, rdStats["broken"]++; rdVerify[rb, d, "line break of " <> lhs]];
  lhs <> " &= " <> lines[[1]] <> StringJoin[("\\\\\n  &\\quad " <> #) & /@ Rest[lines]]];

rdSparsTeX[{lab_, m_}] := lab <> " = \\begin{pmatrix}" <>
  StringRiffle[StringRiffle[If[#, "\\bullet", "\\cdot"] & /@ #, " & "] & /@ m, " \\\\ "] <> "\\end{pmatrix}";

(* a heading with math gets a plain bookmark string *)
rdHead[t_String] := StringReplace[t, "$" ~~ m : Except["$"] .. ~~ "$" :>
  "\\texorpdfstring{$" <> m <> "$}{" <> StringReplace[m, {"\\delta" -> "\\unichar{\"03B4}",
    "\\pi" -> "\\unichar{\"03C0}", "\\alpha" -> "\\unichar{\"03B1}", "\\nabla" -> "\\unichar{\"2207}",
    "_" -> "", "{" -> "", "}" -> "", "\\" -> ""}] <> "}"];
rdItemTeX[it_] := Switch[First[it],
  "sec", "\\section{" <> rdHead[it[[2]]] <> "}\n\n",
  "sub", "\\subsection{" <> rdHead[it[[2]]] <> "}\n\n",
  "subsub", "\\subsubsection{" <> rdHead[it[[2]]] <> "}\n\n",
  "text", it[[2]] <> "\n\n",
  "eq", If[Length[it[[2]]] == 1, "\\[\n" <> it[[2, 1]] <> "\n\\]\n",
    "\\begin{gather*}\n" <> StringRiffle[it[[2]], " \\\\\n"] <> "\n\\end{gather*}\n"],
  "entries", If[it[[2]] === {}, "",
    "\\begin{align*}\n" <> StringRiffle[rdEntryTeX /@ it[[2]], " \\\\\n"] <> "\n\\end{align*}\n"],
  "sparsity", Module[{ps = it[[2]], per},
    per = With[{nf = Max[Length /@ ps[[All, 2]]]}, Which[nf <= 5, 4, nf <= 9, 3, True, 2]];
    StringJoin[("\\[\n\\arraycolsep=2.5pt\n" <> StringRiffle[rdSparsTeX /@ #, ",\\qquad "] <> "\n\\]\n") & /@ Partition[ps, UpTo[per]]]],
  "checks", Module[{log = it[[2]], skip = it[[3]], np, nf},
    np = Count[log, {_, True, _}]; nf = Count[log, {_, False, _}];
    "{\\small\n\\begin{longtable}{@{}>{\\raggedright\\arraybackslash}p{0.27\\textwidth}" <>
      ">{\\raggedright\\arraybackslash}p{0.60\\textwidth}l@{}}\n\\toprule\nid & what it tests & \\\\\n\\midrule\n\\endhead\n" <>
      StringJoin[("\\texttt{" <> rdEscId[#[[1]]] <> "} & " <> rdMsg[#[[3]]] <> " & " <>
        If[#[[2]], "PASS", "\\textbf{FAIL}"] <> " \\\\\n") & /@ log] <>
      "\\bottomrule\n\\end{longtable}\n}\n" <>
      ToString[np] <> " passed, " <> ToString[nf] <> " failed (" <> ToString[np + nf] <> " checks evaluated)." <>
      If[skip === {}, "\n\n",
        " Skipped: " <> StringRiffle[("\\texttt{" <> rdEscId[#[[1]]] <> "}") & /@ skip, ", "] <> ". " <>
        StringRiffle[DeleteDuplicates[skip[[All, 2]]], " "] <> "\n\n"]],
  "source", "Re-run with \\texttt{wolframscript -file deriv/" <> rdEscId[FileNameTake[it[[2, 1]]]] <> "}.\n\n",
  _, ""];

(* ------------------------------------------------------------------------------
   5. Items -> notebook cells
   ------------------------------------------------------------------------------ *)
rdPlain[t_String] := StringReplace[t, {"\\texttt{" ~~ x : Except["}"] .. ~~ "}" :> StringReplace[x, "\\_" -> "_"],
  "---" -> "\[LongDash]", "--" -> "\[Dash]", "~" -> " ", "\\," -> " ", "\\%" -> "%", "\\&" -> "&",
  "\\_" -> "_", "\\ " -> " ", "``" -> "\[OpenCurlyDoubleQuote]", "''" -> "\[CloseCurlyDoubleQuote]"}];
rdInline[s_String] := Module[{parts = StringSplit[s, "$", All]},
  Flatten @ MapIndexed[If[OddQ[First[#2]], If[#1 === "", {}, rdPlain[#1]],
      Cell[BoxData[FormBox[rdTB[#1], TraditionalForm]]]] &, parts]];
rdDF[b_] := Cell[BoxData[FormBox[b, TraditionalForm]], "DisplayFormula"];

rdItemNB[it_] := Switch[First[it],
  "sec", {Cell[TextData[rdInline[it[[2]]]], "Title"]},
  "sub", {Cell[TextData[rdInline[it[[2]]]], "Section"]},
  "subsub", {Cell[TextData[rdInline[it[[2]]]], "Subsection"]},
  "text", {Cell[TextData[rdInline[it[[2]]]], "Text"]},
  "eq", rdDF[rdTB[#]] & /@ it[[2]],
  "entries", rdDF[RowBox[{rdTB[#[[1]]], "=", rdBoxes[#[[2]]]}]] & /@ it[[2]],
  "sparsity", {rdDF[RowBox[Riffle[(RowBox[{rdTB[#[[1]]], "=", RowBox[{"(", GridBox[
        Map[If[#, "\[Bullet]", "\[CenterDot]"] &, #[[2]], {2}]], ")"}]}]) & /@ it[[2]], "   "]]]},
  "checks", Module[{log = it[[2]], skip = it[[3]]},
    Join[{Cell[BoxData[ToBoxes[Grid[Prepend[{#[[1]], #[[3]], If[#[[2]], "PASS", Style["FAIL", Bold, Red]]} & /@ log,
          Style[#, Bold] & /@ {"id", "what it tests", ""}], Frame -> All, Alignment -> Left,
          ItemSize -> {{Automatic, 45, Automatic}}], StandardForm]], "Output"],
      Cell[ToString[Count[log, {_, True, _}]] <> " passed, " <> ToString[Count[log, {_, False, _}]] <>
        " failed (" <> ToString[Length[log]] <> " checks evaluated).", "Text"]},
      If[skip === {}, {}, {Cell["Skipped: " <> StringRiffle[skip[[All, 1]], ", "] <> ". " <>
        StringRiffle[DeleteDuplicates[skip[[All, 2]]], " "], "Text"]}]]],
  "source", Join[{Cell["Source", "Section"],
      Cell[TextData[{"The derivation is re-run with ",
        StyleBox["wolframscript -file deriv/" <> FileNameTake[it[[2, 1]]], FontFamily -> "Courier"],
        If[Length[it[[2]]] > 1, " (which loads the files below it).", "."]}], "Text"]},
    Flatten[{Cell[FileNameTake[#], "Subsection"], Cell[Import[#, "Text"], "Program"]} & /@ it[[2]]]],
  _, {}];

(* ------------------------------------------------------------------------------
   6. Writing a derivation
   ------------------------------------------------------------------------------ *)
rdWrite[name_String, title_String, items_List] := Module[{tex, nb, f, cells, t0 = AbsoluteTime[]},
  tex = "% GENERATED by deriv/render.wl from the derivation script.  Do not edit.\n\n" <>
    StringJoin[rdItemTeX /@ items];
  f = FileNameJoin[{rdOutDir, name <> ".tex"}];
  Export[f, tex, "Text"];
  Print["  [render] wrote ", f, "  (", StringLength[tex], " characters)"];
  cells = Join[{Cell[title, "Title"]}, Flatten[rdItemNB /@ Select[items, First[#] =!= "sec" &]]];
  nb = Notebook[cells, StyleDefinitions -> "Default.nb"];
  f = FileNameJoin[{rdOutDir, name <> ".nb"}];
  Export[f, nb, "NB"];
  Print["  [render] wrote ", f, "  (", FileByteCount[f], " bytes; loads as ", Head[Get[f]], ")"];
  Print["  [render] ", rdStats["entries"], " entries typeset, ", rdStats["verified"],
    " exact-point comparisons, ", rdStats["broken"], " entries broken over lines, ",
    Length[rdProblems], " problems  (", Round[AbsoluteTime[] - t0], " s to write)"];
  If[rdProblems =!= {}, Print["  [render] PROBLEMS:"]; Scan[Print["    ", #] &, rdProblems]];
  ];

(* entries of a matrix / vector: nonzero only, labelled (A_t)_{ij} *)
rdIdx[i_, j_, n_] := If[n >= 10, ToString[i] <> "," <> ToString[j], ToString[i] <> ToString[j]];
rdMat[M_, lab_, prep_, cols_ : All] := Module[{n = Length[M], js},
  js = If[cols === All, Range[Length[First[M]]], cols];
  Flatten[Table[With[{d = prep[M[[i, j]]]},
      If[d === 0, Nothing, {"(" <> lab <> ")_{" <> rdIdx[i, j, Max[n, Length[First[M]]]] <> "}", d}]],
    {i, n}, {j, js}], 1]];
rdVec[V_, lab_, prep_] := Table[With[{d = prep[V[[i]]]}, If[d === 0, Nothing, {lab <> "_{" <> ToString[i] <> "}", d}]],
  {i, Length[V]}];
rdNZ[M_] := Map[# =!= 0 &, M, {2}];
rdTuple[l_List] := "(" <> StringRiffle[l, ", "] <> ")";
rdDiag[geo_] := "x^{\\mu} = " <> rdTuple[rdTeX /@ geo["x"]] <> ",\\qquad g_{\\mu\\nu} = \\mathrm{diag}" <>
  rdTuple[rdTeX /@ Diagonal[geo["g"]]];
rdChecks[] := {"checks", gateLog, {}};

(* symbols common to every derivation *)
rdBaseTeX = <|"tau" -> "\\tau", "utau" -> "u^{\\tau}", "xx" -> "x", "yy" -> "y", "rr" -> "r", "r" -> "r",
  "ph" -> "\\phi", "phi" -> "\\phi", "et" -> "\\eta", "eta" -> "\\eta"|>;
rdDeriv[pre_, f_] := "(\\partial_{" <> pre <> "} " <> f <> ")";

(* ------------------------------------------------------------------------------
   7. Steps: the script's intermediate objects.  Field heads get a symbol and a LaTeX
      name for the field and for each first derivative.
   ------------------------------------------------------------------------------ *)
rdFieldRules[fs_, cs_] := Flatten @ Table[
  With[{h = f[[1]], nm = f[[2]], tx = f[[3]]},
    rdSymTeX[nm] = tx;
    Append[Table[With[{dn = "d" <> ToString[k] <> nm, ord = UnitVector[Length[cs], k]},
        rdSymTeX[dn] = rdDeriv[rdTeX[cs[[k]]], tx];
        (Derivative @@ ord)[h] @@ cs -> Symbol[dn]], {k, Length[cs]}],
      h @@ cs -> Symbol[nm]]], {f, fs}];
rdGammaEq[geo_] := With[{x = rdTeX /@ geo["x"], G = geo["Gam"]},
  StringRiffle[Flatten @ Table[If[G[[a, b, c]] === 0 || b > c, Nothing,
      "\\Gamma^{" <> x[[a]] <> "}{}_{" <> x[[b]] <> x[[c]] <> "} = " <> rdTeX[G[[a, b, c]]]],
    {a, 4}, {b, 4}, {c, 4}], ",\\quad "]];
rdPrs[lab_, e_, pr_] := With[{d = pr[e]}, If[d === 0, Nothing, {lab, d}]];
(* the nonzero components a <= b of a symmetric tensor *)
rdComps[M_, lab_, x_, pr_] := Flatten[Table[If[a > b, Nothing,
    rdPrs[lab <> "^{" <> x[[a]] <> x[[b]] <> "}", M[[a, b]], pr]], {a, 4}, {b, 4}], 1];
(* Nabla_m T^{mn} written out for a generic symmetric tensor with the nonzero pattern of M, with the
   script's own divergence: {entries} for the rows `idx`.  sym gives the component's LaTeX name. *)
rdDivForm[geo_, cs_, M_, idx_, sym_, divFn_] := Module[{x = rdTeX /@ geo["x"], gen, fs, rl, d, tag = ToString[Unique["zT"]]},
  fs = Flatten[Table[If[a > b || M[[a, b]] === 0, Nothing,
      {Symbol[tag <> ToString[a] <> ToString[b]], tag <> "c" <> ToString[a] <> ToString[b], sym <> "^{" <> x[[a]] <> x[[b]] <> "}"}],
      {a, 4}, {b, 4}], 1];
  gen = Table[If[M[[a, b]] === 0, 0, Symbol[tag <> ToString[Min[a, b]] <> ToString[Max[a, b]]] @@ cs], {a, 4}, {b, 4}];
  rl = rdFieldRules[fs, cs];
  d = divFn[geo, gen];
  Table[{"\\nabla_{\\mu}" <> sym <> "^{\\mu " <> x[[n]] <> "}", (d[[n]] //. rl)}, {n, idx}]];
(* Nabla_l S^{l a b} written out for a generic totally symmetric tensor with the nonzero pattern of S *)
rdDivS[geo_, cs_, S_, idx_] := Module[{x = rdTeX /@ geo["x"], gen, fs, rl, d, tag = ToString[Unique["zS"]], key},
  key[a_, b_, c_] := StringJoin[ToString /@ Sort[{a, b, c}]];
  fs = DeleteDuplicates @ Flatten[Table[If[S[[a, b, c]] === 0, Nothing,
      With[{s = Sort[{a, b, c}]}, {Symbol[tag <> key[a, b, c]], tag <> "c" <> key[a, b, c], "S^{" <> x[[s[[1]]]] <> x[[s[[2]]]] <> x[[s[[3]]]] <> "}"}]],
      {a, 4}, {b, 4}, {c, 4}], 2];
  gen = Table[If[S[[a, b, c]] === 0, 0, Symbol[tag <> key[a, b, c]] @@ cs], {a, 4}, {b, 4}, {c, 4}];
  rl = rdFieldRules[fs, cs];
  d = divS3[geo, gen];
  Table[{"\\nabla_{\\lambda}S^{\\lambda " <> x[[ab[[1]]]] <> x[[ab[[2]]]] <> "}", d[[ab[[1]], ab[[2]]]] //. rl}, {ab, idx}]];
rdExtractEq[dirs_] := "A^{(k)}_{ij} = \\frac{\\partial E_{i}}{\\partial(\\partial_{k}\\phi_{j})},\\qquad S_{i} = E_{i}\\bigl|_{\\partial\\phi = 0},\\qquad k = " <> dirs;
(* A.dphi + S against the raw rows, at random exact points *)
rdRebuild[As_, src_, vars_, cs_, eqs_, rules_, what_] := Module[{r = Sum[As[[k]] . (D[#, cs[[k]]] & /@ vars), {k, Length[cs]}] + src},
  Do[rdVerify[r[[i]] //. rules, eqs[[i]] //. rules, what <> ", row " <> ToString[i]], {i, Length[eqs]}]];
(* denominators too wide for a line get a name: {rules, definitions}.  The caller prints the definitions
   and applies the rules to its entries; the renaming is checked by putting them back. *)
rdDenNames[ents_] := Module[{bases, big, nms},
  bases = If[Head[#] === Power, First[#], #] & /@ Flatten[With[{d = Denominator[#[[2]]]}, If[Head[d] === Times, List @@ d, {d}]] & /@ ents];
  big = Select[DeleteDuplicates[bases], Head[#] === Plus && rdWidth[rdTeX[#]] > rdW - 8 &];
  nms = Table[With[{nm = "zDen" <> ToString[k]}, rdSymTeX[nm] = "\\mathcal{D}_{" <> ToString[k] <> "}"; Symbol[nm]], {k, Length[big]}];
  {Thread[big -> nms], Thread[{("\\mathcal{D}_{" <> ToString[#] <> "}") & /@ Range[Length[big]], big}]}];
rdDenApply[ents_, rules_] := ({#[[1]], With[{e = #[[2]] /. rules},
    rdVerify[e /. (Reverse /@ rules), #[[2]], "denominator name in " <> #[[1]]]; e]} & /@ ents);
(* a long field list over several display lines *)
rdListLines[pre_String, l_List] := Module[{lines = {}, cur = {}, w = rdWidth[pre]},
  Do[If[w + rdWidth[s] + 2 > rdW && cur =!= {},
      AppendTo[lines, cur]; cur = {}; w = 0];
    AppendTo[cur, s]; w += rdWidth[s] + 2, {s, l}];
  AppendTo[lines, cur];
  MapIndexed[If[First[#2] == 1, pre <> "(", ""] <> StringRiffle[#1, ", "] <>
    If[First[#2] == Length[lines], ")", ","] &, lines]];

(* ============================================================================
   derive_2p1d_viscous.wls
   ============================================================================ *)
renderViscous[] := Module[{ut1, ut2, p1, p2, items, prepNo, xC, xK, rsC, rsK, s1, s2, DpiC, ProjC, ProjF, DlC, DudC, DloC, DpiK, ProjK},
  rdSymTeX = Join[rdBaseTeX, <|"uT" -> "T", "uR" -> "u^{r}", "uX" -> "u^{x}", "uY" -> "u^{y}",
    "pFF" -> "\\pi^{\\phi}{}_{\\phi}", "pEE" -> "\\pi^{\\eta}{}_{\\eta}", "pBB" -> "\\Pi",
    "pYY" -> "\\pi^{yy}", "pXY" -> "\\pi^{xy}", "pres" -> "p", "dtp" -> "p'(T)", "dtdtp" -> "p''(T)",
    "R" -> "r", "visc" -> "\\eta", "zeta" -> "\\zeta", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}",
    "dpp" -> "\\delta_{\\pi\\pi}"|>];
  ut1 = {{Power[1 + uR^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + uR^2 :> utau^2}, {utau -> Sqrt[1 + uR^2]}};
  ut2 = {{Power[1 + uX^2 + uY^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + uX^2 + uY^2 :> utau^2}, {utau -> Sqrt[1 + uX^2 + uY^2]}};
  p1 = rdPrep[#, subC, ut1, "1+1D"] &;
  p2 = rdPrep[#, subK, ut2, "2+1D"] &;
  rsC = Join[rdFieldRules[{{TT, "uT", "T"}, {ur, "uR", "u^{r}"}, {pff, "pFF", "\\pi^{\\phi}{}_{\\phi}"},
      {pee, "pEE", "\\pi^{\\eta}{}_{\\eta}"}, {PiB, "pBB", "\\Pi"}}, {tau, rr}], subC];
  rsK = Join[rdFieldRules[{{TT, "uT", "T"}, {ux, "uX", "u^{x}"}, {uy, "uY", "u^{y}"}, {pyy, "pYY", "\\pi^{yy}"},
      {pee, "pEE", "\\pi^{\\eta}{}_{\\eta}"}, {pxy, "pXY", "\\pi^{xy}"}, {PiB, "pBB", "\\Pi"}}, {tau, xx, yy}], subK];
  s1 = rdPrep[#, rsC, ut1, "steps 1+1D"] &;
  s2 = rdPrep[#, rsK, ut2, "steps 2+1D"] &;
  xC = rdTeX /@ CYL["x"]; xK = rdTeX /@ CAR["x"];
  DlC = ISC["Dl"]; DudC = ISC["Dud"]; DloC = CYL["g"] . DlC . CYL["g"];
  DpiC = DcovT[CYL, uC, piC];
  ProjF[a_, b_] := Sum[((1/2) (DudC[[a, c]] DudC[[b, d]] + DudC[[a, d]] DudC[[b, c]]) - (1/3) DlC[[a, b]] DloC[[c, d]]) DpiC[[c, d]], {c, 4}, {d, 4}];
  ProjC = Table[ProjF[a, b], {a, 4}, {b, 4}];
  (* the rows from these pieces, and from the extracted matrices *)
  rdVerify[(tauS ProjC[[3, 3]] + piC[[3, 3]] + 2 visc ISC["sig"][[3, 3]]) //. rsC, ISC["Esh"][[3, 3]] //. rsC, "steps: E^{phi phi}"];
  rdVerify[(tauS ProjC[[4, 4]] + piC[[4, 4]] + 2 visc ISC["sig"][[4, 4]]) //. rsC, ISC["Esh"][[4, 4]] //. rsC, "steps: E^{eta eta}"];
  rdRebuild[AsC, srcC, varsC, {tau, rr}, eqsC, rsC, "steps: A dphi + S == E (1+1D)"];
  rdRebuild[AsK, srcK, varsK, {tau, xx, yy}, eqsK, rsK, "steps: A dphi + S == E (2+1D)"];
  DpiK = DcovT[CAR, uK, piK];
  ProjK = With[{Dl = ISK["Dl"], Du = ISK["Dud"], Dlo = CAR["g"] . ISK["Dl"] . CAR["g"]},
    Table[Sum[((1/2) (Du[[a, c]] Du[[b, d]] + Du[[a, d]] Du[[b, c]]) - (1/3) Dl[[a, b]] Dlo[[c, d]]) DpiK[[c, d]], {c, 4}, {d, 4}], {a, 4}, {b, 4}]];
  Do[rdVerify[(tauS ProjK[[ab[[1]], ab[[2]]]] + piK[[ab[[1]], ab[[2]]]] + 2 visc ISK["sig"][[ab[[1]], ab[[2]]]]) //. rsK,
      ISK["Esh"][[ab[[1]], ab[[2]]]] //. rsK, "steps: 2+1D shear row " <> ToString[ab]], {ab, {{3, 3}, {4, 4}, {2, 3}}}];
  items = {
    {"sec", "Viscous hydrodynamics: Israel--Stewart with $\\delta_{\\pi\\pi}$"},
    {"text", "Script: \\texttt{deriv/derive\\_2p1d\\_viscous.wls}"}, {"text", "Writes: \\texttt{two\\_d\\_viscous\\_matrix\\_derived} in \\texttt{src/Matrix/2d\\_viscous.jl}; the 1+1D system is the control, compared entry by entry with \\texttt{one\\_d\\_viscous\\_matrix} (V-G4)."},
    {"sub", "Setup"},
    {"text", "Milne coordinates, boost invariant. Two charts:"},
    {"eq", {rdDiag[CYL], rdDiag[CAR]}},
    {"eq", {
      "T^{\\mu\\nu} = (e + p + \\Pi)\\,u^{\\mu}u^{\\nu} + (p + \\Pi)\\,g^{\\mu\\nu} + \\pi^{\\mu\\nu},\\qquad e = T\\,p'(T) - p(T)",
      "\\nabla_{\\mu}T^{\\mu\\nu} = 0",
      "\\tau_{\\pi}\\,\\Delta^{\\mu\\nu}_{\\alpha\\beta}\\,D\\pi^{\\alpha\\beta} + \\pi^{\\mu\\nu} + 2\\eta\\,\\sigma^{\\mu\\nu} + \\delta_{\\pi\\pi}\\,\\pi^{\\mu\\nu}\\theta = 0",
      "\\tau_{\\Pi}\\,D\\Pi + \\Pi + \\zeta\\,\\theta = 0"}},
    {"text", "with"},
    {"eq", {
      "D = u^{\\mu}\\nabla_{\\mu},\\qquad \\theta = \\nabla_{\\mu}u^{\\mu},\\qquad \\Delta^{\\mu\\nu} = g^{\\mu\\nu} + u^{\\mu}u^{\\nu}",
      "\\sigma^{\\mu\\nu} = \\Delta^{\\mu\\nu}_{\\alpha\\beta}\\nabla^{\\alpha}u^{\\beta},\\qquad \\Delta^{\\mu\\nu}_{\\alpha\\beta} = \\tfrac{1}{2}\\bigl(\\Delta^{\\mu}_{\\alpha}\\Delta^{\\nu}_{\\beta} + \\Delta^{\\mu}_{\\beta}\\Delta^{\\nu}_{\\alpha}\\bigr) - \\tfrac{1}{3}\\Delta^{\\mu\\nu}\\Delta_{\\alpha\\beta}"}},
    {"text", "$\\eta$, $\\tau_{\\pi}$, $\\zeta$, $\\tau_{\\Pi}$ are not differentiated. Quasi-linear form, source on the left:"},
    {"eq", {"A_t\\,\\partial_{\\tau}\\phi + A_r\\,\\partial_{r}\\phi + S = 0\\quad\\text{(1+1D)}", "A_t\\,\\partial_{\\tau}\\phi + A_x\\,\\partial_{x}\\phi + A_y\\,\\partial_{y}\\phi + S = 0\\quad\\text{(2+1D)}"}},
    {"text", "The $\\delta_{\\pi\\pi}$ term is linear, so each matrix and the source split as"},
    {"eq", {"A = A^{\\mathrm{IS}} + \\delta_{\\pi\\pi}\\,A^{(\\delta)},\\qquad S = S^{\\mathrm{IS}} + \\delta_{\\pi\\pi}\\,S^{(\\delta)}"}},
    {"sub", "Steps"},
    {"subsub", "1+1D: Christoffels and constraints"},
    {"eq", {rdGammaEq[CYL]}},
    {"entries", {rdPrs["u^{\\tau}", uC[[1]], (# //. rsC &)]}},
    {"text", "$u_{\\mu}\\pi^{\\mu\\nu} = 0$, $\\pi^{\\mu}{}_{\\mu} = 0$:"},
    {"entries", {rdPrs["\\pi^{rr}", prrC, s1], rdPrs["\\pi^{\\tau r}", ptrC, s1], rdPrs["\\pi^{\\tau\\tau}", pttC, s1],
      rdPrs["\\pi^{\\phi\\phi}", piC[[3, 3]], s1], rdPrs["\\pi^{\\eta\\eta}", piC[[4, 4]], s1]}},
    {"subsub", "1+1D: velocity gradients"},
    {"entries", Join[{rdPrs["\\theta", ISC["theta"], s1]},
      Flatten[Table[rdPrs["\\nabla_{" <> xC[[l]] <> "}u^{" <> xC[[m]] <> "}", gradV[CYL, uC][[l, m]], s1], {l, 4}, {m, 4}], 1],
      rdComps[ISC["sig"], "\\sigma", xC, s1]]},
    {"subsub", "1+1D: conservation"},
    {"entries", rdComps[ISC["T"], "T", xC, s1]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivForm[CYL, {tau, rr}, ISC["T"], {1, 2}, "T", divT]},
    {"subsub", "1+1D: shear and bulk"},
    {"entries", rdComps[DpiC, "(D\\pi)", xC, s1]},
    {"entries", {rdPrs["\\Delta^{\\phi\\phi}_{\\alpha\\beta}D\\pi^{\\alpha\\beta}", ProjC[[3, 3]], s1],
      rdPrs["\\Delta^{\\eta\\eta}_{\\alpha\\beta}D\\pi^{\\alpha\\beta}", ProjC[[4, 4]], s1]}},
    {"eq", {"E^{\\phi\\phi} = \\tau_{\\pi}\\Delta^{\\phi\\phi}_{\\alpha\\beta}D\\pi^{\\alpha\\beta} + \\pi^{\\phi\\phi} + 2\\eta\\sigma^{\\phi\\phi} + \\delta_{\\pi\\pi}\\pi^{\\phi\\phi}\\theta,\\quad \\phi\\to\\eta",
      "E_{\\Pi} = \\tau_{\\Pi}(u^{\\tau}\\partial_{\\tau} + u^{r}\\partial_{r})\\Pi + \\Pi + \\zeta\\theta"}},
    {"subsub", "1+1D: from the rows to the matrices"},
    {"eq", {rdExtractEq["\\tau, r"]}},
    {"entries", {rdPrs["E_{1}", eqsC[[1]], s1], rdPrs["E_{5}", eqsC[[5]], s1]}},
    {"subsub", "2+1D"},
    {"eq", {rdGammaEq[CAR]}},
    {"text", "$u_{\\mu}\\pi^{\\mu\\nu} = 0$, $\\pi^{\\mu}{}_{\\mu} = 0$:"},
    {"entries", {rdPrs["\\pi^{xx}", pxxK, s2], rdPrs["\\pi^{\\tau x}", ptxK, s2], rdPrs["\\pi^{\\tau y}", ptyK, s2], rdPrs["\\pi^{\\tau\\tau}", pttK, s2]}},
    {"entries", Join[{rdPrs["\\theta", ISK["theta"], s2]},
      Flatten[Table[rdPrs["\\nabla_{" <> xK[[l]] <> "}u^{" <> xK[[m]] <> "}", gradV[CAR, uK][[l, m]], s2], {l, 4}, {m, 4}], 1],
      rdComps[ISK["sig"], "\\sigma", xK, s2]]},
    {"entries", rdComps[ISK["T"], "T", xK, s2]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivForm[CAR, {tau, xx, yy}, ISK["T"], {1, 2, 3}, "T", divT]},
    {"entries", rdComps[DpiK, "(D\\pi)", xK, s2]},
    {"entries", {rdPrs["\\Delta^{yy}_{\\alpha\\beta}D\\pi^{\\alpha\\beta}", ProjK[[3, 3]], s2], rdPrs["\\Delta^{\\eta\\eta}_{\\alpha\\beta}D\\pi^{\\alpha\\beta}", ProjK[[4, 4]], s2],
      rdPrs["\\Delta^{xy}_{\\alpha\\beta}D\\pi^{\\alpha\\beta}", ProjK[[2, 3]], s2]}},
    {"eq", {"E^{ab} = \\tau_{\\pi}\\Delta^{ab}_{\\alpha\\beta}D\\pi^{\\alpha\\beta} + \\pi^{ab} + 2\\eta\\sigma^{ab} + \\delta_{\\pi\\pi}\\pi^{ab}\\theta,\\qquad ab = yy,\\ \\eta\\eta,\\ xy",
      "E_{\\Pi} = \\tau_{\\Pi}(u^{\\tau}\\partial_{\\tau} + u^{x}\\partial_{x} + u^{y}\\partial_{y})\\Pi + \\Pi + \\zeta\\theta", rdExtractEq["\\tau, x, y"]}},
    {"entries", {rdPrs["E_{1}", eqsK[[1]], s2], rdPrs["E_{7}", eqsK[[7]], s2]}},
    {"sub", "1+1D cylindrical"},
    {"text", "Fields, with $\\pi^{rr}$ and the $\\tau$ row of $\\pi^{\\mu\\nu}$ fixed by $u_{\\mu}\\pi^{\\mu\\nu} = 0$, $\\pi^{\\mu}{}_{\\mu} = 0$:"},
    {"eq", {"\\phi = " <> rdTuple[rdTeX /@ (varsC //. subC)] <> ",\\qquad \\pi^{\\phi}{}_{\\phi} = r^{2}\\pi^{\\phi\\phi},\\quad \\pi^{\\eta}{}_{\\eta} = \\tau^{2}\\pi^{\\eta\\eta}",
      "\\pi^{rr} = " <> rdTeX[p1[prrC]] <> ",\\qquad u^{\\tau} = " <> rdTeX[Sqrt[1 + uR^2]]}},
    {"text", "Rows: $\\nabla_{\\mu}T^{\\mu\\tau}$, $\\nabla_{\\mu}T^{\\mu r}$, $E^{\\phi\\phi}$, $E^{\\eta\\eta}$, $E_{\\Pi}$, where $E^{\\mu\\nu}$ and $E_{\\Pi}$ are the left-hand sides of the shear and bulk equations."},
    {"subsub", "IS part"},
    {"entries", rdMat[AtC, "A^{\\mathrm{IS}}_t", p1]},
    {"entries", rdMat[ArC, "A^{\\mathrm{IS}}_r", p1]},
    {"entries", rdVec[srcC, "S^{\\mathrm{IS}}", p1]},
    {"subsub", "$\\delta_{\\pi\\pi}$ part"},
    {"entries", rdMat[AtDC, "A^{(\\delta)}_t", p1]},
    {"entries", rdMat[ArDC, "A^{(\\delta)}_r", p1]},
    {"entries", rdVec[srcDC, "S^{(\\delta)}", p1]},
    {"sub", "2+1D Cartesian"},
    {"text", "Fields, with $\\pi^{xx}$ and the $\\tau$ row fixed by the same constraints:"},
    {"eq", {"\\phi = " <> rdTuple[rdTeX /@ (varsK //. subK)],
      "\\pi^{xx} = " <> rdTeX[p2[pxxK]] <> ",\\qquad u^{\\tau} = " <> rdTeX[Sqrt[1 + uX^2 + uY^2]]}},
    {"text", "Rows: $\\nabla_{\\mu}T^{\\mu\\tau}$, $\\nabla_{\\mu}T^{\\mu x}$, $\\nabla_{\\mu}T^{\\mu y}$, $E^{yy}$, $E^{\\eta\\eta}$, $E^{xy}$, $E_{\\Pi}$."},
    {"subsub", "IS part"},
    {"entries", rdMat[AtK, "A^{\\mathrm{IS}}_t", p2]},
    {"entries", rdMat[AxK, "A^{\\mathrm{IS}}_x", p2]},
    {"entries", rdMat[AyK, "A^{\\mathrm{IS}}_y", p2]},
    {"entries", rdVec[srcK, "S^{\\mathrm{IS}}", p2]},
    {"subsub", "$\\delta_{\\pi\\pi}$ part"},
    {"entries", rdMat[AtDK, "A^{(\\delta)}_t", p2]},
    {"entries", rdMat[AxDK, "A^{(\\delta)}_x", p2]},
    {"entries", rdMat[AyDK, "A^{(\\delta)}_y", p2]},
    {"entries", rdVec[srcDK, "S^{(\\delta)}", p2]},
    {"sub", "Checks"},
    {"checks", gateLog, {}},
    {"text", "Not checked in this script: the slab limit (2+1D at $u^{y} = 0$, $\\partial_{y} = 0$ against 1+1D at $r \\to \\infty$); it is a numerical comparison of the two generated routines."},
    {"source", {FileNameJoin[{DERIV, "derive_2p1d_viscous.wls"}]}}};
  rdWrite["viscous", "Viscous hydrodynamics: Israel\[Dash]Stewart with \[Delta]\[Pi]\[Pi]", items]];

(* ============================================================================
   derive_hq_2p1d.wls and derive_hq_consistent_2p1d.wls
   ============================================================================ *)
rdHQTeX[] := Join[rdBaseTeX, <|"ux" -> "u^{x}", "uy" -> "u^{y}", "T" -> "T", "al" -> "\\alpha",
  "nuxF" -> "\\nu^{x}", "nuyF" -> "\\nu^{y}", "n" -> "n", "rdDnda" -> "(\\partial_{\\alpha} n)",
  "rdDndT" -> "(\\partial_{T} n)", "taun" -> "\\tau_{n}", "kappa" -> "\\kappa", "h" -> "h", "hp" -> "h'(T)",
  "Ds" -> "D_{s}",
  "dtT" -> rdDeriv["\\tau", "T"], "dxT" -> rdDeriv["x", "T"], "dyT" -> rdDeriv["y", "T"],
  "dtux" -> rdDeriv["\\tau", "u^{x}"], "dxux" -> rdDeriv["x", "u^{x}"], "dyux" -> rdDeriv["y", "u^{x}"],
  "dtuy" -> rdDeriv["\\tau", "u^{y}"], "dxuy" -> rdDeriv["x", "u^{y}"], "dyuy" -> rdDeriv["y", "u^{y}"]|>];
rdUtCar = {{Power[1 + ux^2 + uy^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + ux^2 + uy^2 :> utau^2}, {utau -> Sqrt[1 + ux^2 + uy^2]}};
rdHQSetupCommon[] := {
  {"text", "Cartesian Milne coordinates, boost invariant, on a prescribed background $T$, $u^{x}$, $u^{y}$:"},
  {"eq", {rdDiag[CAR] <> ",\\qquad u^{\\mu} = (u^{\\tau}, u^{x}, u^{y}, 0),\\quad u^{\\tau} = " <>
    rdTeX[Sqrt[1 + ux^2 + uy^2]]}}};

renderHQ[] := Module[{jr, pr, items, rs, s2, xK, accK, DnuK, ndaK},
  rdSymTeX = rdHQTeX[];
  jr = jrepl /. {Verbatim[dn_da] -> rdDnda, Verbatim[dn_dT] -> rdDndT};
  pr = rdPrep[#, jr, rdUtCar, "2+1D"] &;
  rs = Join[rdFieldRules[{{al2, "al", "\\alpha"}, {nux, "nuxF", "\\nu^{x}"}, {nuy, "nuyF", "\\nu^{y}"}}, {tau, xx, yy}], jr];
  s2 = rdPrep[#, rs, rdUtCar, "steps 2+1D"] &;
  xK = rdTeX /@ CAR["x"];
  accK = DcovV[CAR, ucar, ucar]; DnuK = DcovV[CAR, ucar, nucar];
  ndaK = Sum[(CAR["g"] . accK)[[m]] nucar[[m]], {m, 4}];
  rdVerify[(DnuK[[2]] - ucar[[2]] ndaK) taun + nucar[[2]] + kappa Sum[(CAR["ginv"] + Outer[Times, ucar, ucar])[[2, b]] D[al2[tau, xx, yy], CAR["x"][[b]]], {b, 4}] //. rs,
    e2k[[2]] //. rs, "steps: nu^x row"];
  rdRebuild[AsK, srcK, varsK, {tau, xx, yy}, eqsK, rs, "steps: A dphi + S == E"];
  items = {
    {"sec", "Charm first moment on a 2+1D background ($\\nabla\\alpha$ closure)"},
    {"text", "Script: \\texttt{deriv/derive\\_hq\\_2p1d.wls}"}, {"text", "Writes: \\texttt{two\\_d\\_HQ\\_BG\\_matrix} and \\texttt{two\\_d\\_HQ\\_live\\_rows} in \\texttt{src/Matrix/2d\\_viscous\\_HQ\\_BG.jl}"},
    {"sub", "Setup"},
    Splice[rdHQSetupCommon[]],
    {"text", "Charm is a passive tracer:"},
    {"eq", {"N^{\\mu} = n\\,u^{\\mu} + \\nu^{\\mu},\\qquad \\nabla_{\\mu}N^{\\mu} = 0",
      "\\tau_{n}\\,\\Delta^{\\mu}{}_{\\beta}\\,D\\nu^{\\beta} + \\nu^{\\mu} + \\kappa\\,\\Delta^{\\mu\\beta}\\partial_{\\beta}\\alpha = 0",
      "n = n(T, \\alpha),\\qquad \\kappa = D_{s}\\,\\partial_{\\alpha}n,\\qquad D = u^{\\mu}\\nabla_{\\mu},\\qquad \\Delta^{\\mu\\nu} = g^{\\mu\\nu} + u^{\\mu}u^{\\nu}"}},
    {"text", "Navier--Stokes limit: Fick's law $\\nu^{\\mu} = -\\kappa\\,\\nabla^{\\mu}_{\\perp}\\alpha$. Fields, with $\\nu^{\\tau}$ from $u_{\\mu}\\nu^{\\mu} = 0$:"},
    {"eq", {"\\phi = " <> rdTuple[rdTeX /@ (varsK //. jr)] <> ",\\qquad \\nu^{\\tau} = " <> rdTeX[pr[nucar[[1]]]]}},
    {"text", "Rows: $\\nabla_{\\mu}N^{\\mu}$ and the $x$, $y$ components of the current equation, in"},
    {"eq", {"A_t\\,\\partial_{\\tau}\\phi + A_x\\,\\partial_{x}\\phi + A_y\\,\\partial_{y}\\phi + S = 0"}},
    {"sub", "Steps"},
    {"eq", {rdGammaEq[CAR]}},
    {"entries", Join[
      Table[rdPrs["a^{" <> xK[[m]] <> "}", accK[[m]], s2], {m, 4}],
      {rdPrs["\\nu_{\\lambda}a^{\\lambda}", ndaK, s2]},
      Table[rdPrs["D\\nu^{" <> xK[[m]] <> "}", DnuK[[m]], s2], {m, 4}]]},
    {"text", "$\\Delta^{\\mu}{}_{\\beta}D\\nu^{\\beta} = D\\nu^{\\mu} - u^{\\mu}\\nu_{\\lambda}a^{\\lambda}$, because $u_{\\beta}\\nu^{\\beta} = 0$. Rows:"},
    {"entries", {rdPrs["E_{1} = \\nabla_{\\mu}N^{\\mu}", eqsK[[1]], s2], rdPrs["E_{2}", eqsK[[2]], s2], rdPrs["E_{3}", eqsK[[3]], s2]}},
    {"eq", {rdExtractEq["\\tau, x, y"]}},
    {"sub", "Result"},
    {"entries", rdMat[AtK, "A_t", pr]},
    {"entries", rdMat[AxK, "A_x", pr]},
    {"entries", rdMat[AyK, "A_y", pr]},
    {"entries", rdVec[srcK, "S", pr]},
    {"subsub", "Live background"},
    {"text", "The same three rows with the medium promoted to fields, $\\phi_{L} = " <> rdTuple[rdTeX /@ (varsL //. jr)] <>
      "$. Columns 4--6 are the matrices above (H-G6a); columns 1--3 and the remaining source:"},
    {"entries", rdMat[AtL, "A^{L}_t", pr, {1, 2, 3}]},
    {"entries", rdMat[AxL, "A^{L}_x", pr, {1, 2, 3}]},
    {"entries", rdMat[AyL, "A^{L}_y", pr, {1, 2, 3}]},
    {"entries", rdVec[srcL, "S^{L}", pr]},
    {"sub", "Checks"},
    {"checks", gateLog, {}},
    {"source", {FileNameJoin[{DERIV, "derive_hq_2p1d.wls"}]}}};
  rdWrite["hq_2p1d", "Charm first moment on a 2+1D background (\[Del]\[Alpha] closure)", items]];

(* ============================================================================
   derive_two_way.wls (with _gates and _export)
   ============================================================================ *)
rdTwoWayTeX[] := Module[{fld, m},
  fld = <|"T" -> "T", "ur" -> "u^{r}", "ux" -> "u^{x}", "uy" -> "u^{y}",
    "piphi" -> "\\pi^{\\phi}{}_{\\phi}", "pieta" -> "\\pi^{\\eta}{}_{\\eta}", "piyy" -> "\\pi^{yy}",
    "pixy" -> "\\pi^{xy}", "PiB" -> "\\Pi", "al" -> "\\alpha", "nur" -> "\\nu^{r}", "nux" -> "\\nu^{x}",
    "nuy" -> "\\nu^{y}", "piphiphiQ" -> "\\pi_{Q}{}^{\\phi}{}_{\\phi}",
    "pietaetaQ" -> "\\pi_{Q}{}^{\\eta}{}_{\\eta}", "piBQ" -> "\\Pi_{Q}", "piyyQ" -> "\\pi_{Q}^{yy}",
    "pizzQ" -> "\\pi_{Q}{}^{\\eta}{}_{\\eta}", "pixyQ" -> "\\pi_{Q}^{xy}"|>;
  m = Join[rdBaseTeX, fld, <|"pres" -> "P", "dtp" -> "P'(T)", "dtdtp" -> "P''(T)", "aSB" -> "a_{\\mathrm{SB}}",
    "zeta" -> "\\zeta", "visc" -> "\\eta", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}",
    "dpp" -> "\\delta_{\\pi\\pi}", "dPP" -> "\\delta_{\\Pi\\Pi}", "n" -> "n", "dndT" -> "(\\partial_{T} n)",
    "dnda" -> "(\\partial_{\\alpha} n)", "h" -> "h", "dhdT" -> "h'(T)", "I31" -> "I_{31}",
    "dI31dT" -> "(\\partial_{T} I_{31})", "dI31da" -> "(\\partial_{\\alpha} I_{31})", "aQ" -> "\\mathcal{A}",
    "bQ" -> "\\mathcal{B}", "dbQdT" -> "\\mathcal{B}'(T)", "cQ" -> "\\mathcal{C}", "dcQdT" -> "\\mathcal{C}'(T)",
    "Ds" -> "D_{s}", "wv" -> "w_{\\omega}"|>];
  (* derivatives named by the export: dt/dr/dx/dy + field name *)
  Do[m[d[[1]] <> f] = rdDeriv[d[[2]], fld[f]],
    {d, {{"dt", "\\tau"}, {"dr", "r"}, {"dx", "x"}, {"dy", "y"}}}, {f, Keys[fld]}];
  m];

renderTwoWay[] := Module[{items, full, sysLab, chartLab, dirLab, fnItems, fnPrep, stJ, fl, conv, st1, r1, B1, R1, sys1, x1, g1, s1, ch, Dmx,
    lo1, lhs1, tr1, shown, appItems, medL, stK, sK, utT, rK, BK, RK, sysK, gK, xK, utK, DmxK},
  rdSymTeX = rdTwoWayTeX[];
  stJ[st_] := (toJ[#, st] &);
  fl[st_, key_] := rdTuple[rdTeX /@ (fieldNames /@ st[key] /. s_String :> Symbol[s])];
  sysLab = <|"S1" -> "S1, ideal hydro", "S2" -> "S2, viscous hydro", "S3" -> "S3, charm first moment on a background",
    "S4" -> "S4, charm first and second moment on a background", "S5" -> "S5, two-way", "S6" -> "S6, fully coupled"|>;
  chartLab = <|"0p1d" -> "0+1D Bjorken", "1p1d" -> "1+1D cylindrical", "2p1d" -> "2+1D Cartesian"|>;
  dirLab = <|"0p1d" -> {"A_t"}, "1p1d" -> {"A_t", "A_r"}, "2p1d" -> {"A_t", "A_x", "A_y"}|>;
  (* typeset in full: the 1+1D charm systems and the 1+1D conformal two-way system *)
  full[r_] := r["chart"] === "1p1d" && (!MemberQ[{"S5", "S6"}, r["sys"]] || r["conf"]);
  (* the exported entries, with 1 + u_perp^2 written (u^tau)^2 where it stands alone *)
  fnPrep[nm_] := With[{q = If[nm === "1p1d", 1 + ur^2, 1 + ux^2 + uy^2], rt = If[nm === "1p1d", Sqrt[1 + ur^2], Sqrt[1 + ux^2 + uy^2]]},
    Function[e, With[{d = rdFlatten[e /. {x_Plus /; x === q :> utau^2}]},
      If[d =!= e, rdVerify[d /. utau -> rt, e /. utau -> rt, "u^tau in an exported entry"]]; d]]];
  fnItems[fn_, r_, wantFull_] := Module[{nm = r["chart"], labs = dirLab[r["chart"]], head, ents},
    head = {{"subsub", sysLab[r["sys"]] <> If[r["conf"], ", conformal", ""]},
      {"text", "\\texttt{" <> rdEscId[fn] <> "} in \\texttt{src/Matrix/" <> rdEscId[r["file"]] <> "}, fields:"},
      {"eq", rdListLines["\\phi = ", rdTeX /@ (Symbol /@ r["fields"])]}};
    If[wantFull,
      ents = Join[Flatten[Table[rdMat[r["A"][[ic]], labs[[ic]], fnPrep[nm]], {ic, Length[labs]}], 1],
        rdVec[r["src"], "S", fnPrep[nm]]];
      Join[head, {{"entries", ents}}],
      Join[head, {{"sparsity", Append[Table[{labs[[ic]], rdNZ[r["A"][[ic]]]}, {ic, Length[labs]}],
        {"S", List /@ (# =!= 0 & /@ r["src"])}]}}]]];
  (* the steps, 1+1D: the script's tensors in the export's names *)
  r1 = results["1p1d"]; st1 = r1["st"]; B1 = r1["B"]; R1 = r1["R"]; sys1 = r1["sys"]; g1 = st1["geo"]["g"];
  x1 = rdTeX /@ st1["geo"]["x"];
  utT = {{Power[1 + ur^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + ur^2 :> utau^2}, {utau -> Sqrt[1 + ur^2]}};
  s1 = rdPrep[toJ[# /. Derivative[1][aQf][_] -> 0, st1] /. {dPP -> 0}, {}, utT, "steps 1+1D"] &;
  Dmx = (st1["geo"]["ginv"] + Outer[Times, st1["u"], st1["u"]]) . g1;
  lo1 = (g1 . #) & /@ st1["triad"]; lhs1 = B1["LHS"]; tr1 = B1["tr"];
  (* the pieces assemble into the rows, and the extracted matrices back into the rows *)
  rdVerify[toJ[Sum[Dmx[[2, l]] B1["consQ"][[l]], {l, 4}] (Ds/st1["T"]) + st1["nu"][[2]], st1], toJ[R1["ch"][[2]], st1], "steps: nu^r row"];
  rdVerify[toJ[(B1["pre"] (lo1[[2]] . lhs1 . lo1[[2]] - tr1/3) + st1["f"][pffQc]) /. Derivative[1][aQf][_] -> 0, st1], toJ[R1["m2"][[1]], st1], "steps: phi phi row"];
  rdVerify[toJ[(B1["pre"] (lo1[[3]] . lhs1 . lo1[[3]] - tr1/3) + st1["f"][peeQc]) /. Derivative[1][aQf][_] -> 0, st1], toJ[R1["m2"][[2]], st1], "steps: eta eta row"];
  Do[rdRebuild[sys1[k]["A"], sys1[k]["src"], (# @@ st1["cs"]) & /@ sys1[k]["vars"], st1["cs"], sys1[k]["eqs"], jRules[st1],
      "steps: A dphi + S == E, " <> k], {k, {"S3", "S4", "S5", "S6"}}];
  rdVerify[toJ[(Sum[Dmx[[2, l]] B1["consD"][[l]], {l, 4}] (Ds/st1["T"])), st1], toJ[R1["chD"][[2]], st1], "steps: c_M part of the current row"];
  (* the same, 2+1D *)
  rK = results["2p1d"]; stK = rK["st"]; BK = rK["B"]; RK = rK["R"]; sysK = rK["sys"]; gK = stK["geo"]["g"];
  xK = rdTeX /@ stK["geo"]["x"];
  utK = {{Power[1 + ux^2 + uy^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + ux^2 + uy^2 :> utau^2}, {utau -> Sqrt[1 + ux^2 + uy^2]}};
  sK = rdPrep[toJ[# /. Derivative[1][aQf][_] -> 0, stK] /. {dPP -> 0, wv -> 1}, {}, utK, "steps 2+1D"] &;
  DmxK = (stK["geo"]["ginv"] + Outer[Times, stK["u"], stK["u"]]) . gK;
  Do[rdVerify[toJ[Sum[DmxK[[m, l]] BK["consQ"][[l]], {l, 4}] (Ds/stK["T"]) + stK["nu"][[m]], stK], toJ[RK["ch"][[m]], stK], "steps: nu row 2+1D"], {m, {2, 3}}];
  rdVerify[toJ[(BK["pre"] BK["TL"][[3, 3]] + stK["f"][pyyQk]) /. Derivative[1][aQf][_] -> 0, stK], toJ[RK["m2"][[1]], stK] /. wv -> 1, "steps: yy row"];
  rdVerify[toJ[(tau^2 BK["pre"] BK["TL"][[4, 4]] + stK["f"][peeQk]) /. Derivative[1][aQf][_] -> 0, stK], toJ[RK["m2"][[2]], stK] /. wv -> 1, "steps: eta eta row"];
  rdVerify[toJ[(BK["pre"] BK["TL"][[2, 3]] + stK["f"][pxyQk]) /. Derivative[1][aQf][_] -> 0, stK], toJ[RK["m2"][[3]], stK] /. wv -> 1, "steps: xy row"];
  Do[rdRebuild[sysK[k]["A"], sysK[k]["src"], (# @@ stK["cs"]) & /@ sysK[k]["vars"], stK["cs"], sysK[k]["eqs"] /. wv -> 1, Append[jRules[stK], wv -> 1],
      "steps: A dphi + S == E, 2+1D " <> k], {k, {"S3", "S4", "S5", "S6"}}];
  items = {
    {"sec", "Medium, charm and the two-way coupled system"},
    {"text", "Script: \\texttt{deriv/derive\\_two\\_way.wls}, with \\texttt{derive\\_two\\_way\\_gates.wls} and \\texttt{derive\\_two\\_way\\_export.wls}"}, {"text", "Writes: the charm and two-way matrices of \\texttt{src/Matrix} (one function per system, listed below)"},
    {"sub", "Setup"},
    {"text", "Three charts, boost invariant: 0+1D Bjorken (fields of $\\tau$ only), 1+1D cylindrical, 2+1D Cartesian Milne:"},
    {"eq", {rdDiag[BJ] <> ",\\quad u^{\\mu} = " <> rdTuple[rdTeX /@ stateBJ["u"]], rdDiag[CYL], rdDiag[CAR]}},
    {"text", "Medium, with the EOS $P(T)$ kept symbolic:"},
    {"eq", {"T_{m}^{\\mu\\nu} = (e + P + \\Pi)\\,u^{\\mu}u^{\\nu} + (P + \\Pi)\\,g^{\\mu\\nu} + \\pi^{\\mu\\nu},\\qquad e = T P'(T) - P(T)",
      "\\tau_{\\pi}\\,\\Delta^{\\mu\\nu}_{\\alpha\\beta}D\\pi^{\\alpha\\beta} + \\pi^{\\mu\\nu} + 2\\eta\\,\\sigma^{\\mu\\nu} + \\delta_{\\pi\\pi}\\,\\pi^{\\mu\\nu}\\theta = 0,\\qquad \\tau_{\\Pi}D\\Pi + \\Pi + \\zeta\\theta + \\delta_{\\Pi\\Pi}\\Pi\\theta = 0"}},
    {"text", "($D$, $\\theta$, $\\Delta^{\\mu\\nu}$, $\\Delta^{\\mu\\nu}_{\\alpha\\beta}$, $\\sigma^{\\mu\\nu}$ as in the viscous section.) Charm, with $n = n(T, \\alpha)$, $h = h(T)$, $P_{0} = nT$ and $\\delta M^{\\mu\\nu} = \\pi_{Q}^{\\mu\\nu} + \\Pi_{Q}\\Delta^{\\mu\\nu}$:"},
    {"eq", {"\\nabla_{\\mu}(n\\,u^{\\mu} + \\nu^{\\mu}) = 0,\\qquad \\frac{D_{s}}{T}\\,\\Delta^{\\mu}{}_{\\lambda}\\nabla_{\\nu}T_{Q}^{\\nu\\lambda} + \\nu^{\\mu} = 0",
      "T_{Q}^{\\mu\\nu} = (n h - P_{0})\\,u^{\\mu}u^{\\nu} + P_{0}\\,\\Delta^{\\mu\\nu} + h\\,(u^{\\mu}\\nu^{\\nu} + \\nu^{\\mu}u^{\\nu})",
      "S^{\\lambda\\alpha\\beta} = I_{31}\\,(u^{\\lambda}\\Delta^{\\alpha\\beta} + 2\\ \\text{perm.}) + \\mathcal{A}\\,(u^{\\lambda}u^{\\alpha}\\nu^{\\beta} + 2\\ \\text{perm.}) + \\mathcal{B}\\,(\\Delta^{\\lambda\\alpha}\\nu^{\\beta} + 2\\ \\text{perm.}) + \\mathcal{C}\\,(u^{\\lambda}\\delta M^{\\alpha\\beta} + 2\\ \\text{perm.})",
      "\\mathcal{E}^{\\mu\\nu} = \\Delta^{\\mu}{}_{\\alpha}\\Delta^{\\nu}{}_{\\beta}\\nabla_{\\lambda}S^{\\lambda\\alpha\\beta}:\\qquad \\frac{D_{s}}{2T}\\Bigl(\\mathcal{E}^{\\mu\\nu} - \\tfrac{1}{3}\\Delta^{\\mu\\nu}\\mathcal{E}^{\\lambda}{}_{\\lambda}\\Bigr) + \\pi_{Q}^{\\mu\\nu} = 0,\\qquad \\frac{D_{s}}{2T}\\,\\tfrac{1}{3}\\mathcal{E}^{\\lambda}{}_{\\lambda} + \\Pi_{Q} = 0"}},
    {"text", "``2 perm.'' are the two other placements of the lone index. $I_{31}(T, \\alpha)$, $\\mathcal{A}(T)$, $\\mathcal{B}(T)$, $\\mathcal{C}(T)$ stay symbolic (Boltzmann: $I_{31} = nTh$, $\\mathcal{A} = m^{2} + 5Th$, $\\mathcal{B} = Th$, $\\mathcal{C} = m K_{4}/K_{3}$). In S1--S5 the second moment $\\delta M$ is evolved but does not enter $T_{Q}$. In 2+1D the traceless rows also carry"},
    {"eq", {"-(1 - w_{\\omega})\\,\\tau_{M}W^{\\mu\\nu},\\qquad \\tau_{M} = \\frac{D_{s}\\mathcal{C}}{2T},\\qquad W^{\\mu\\nu} = \\pi_{Q}^{\\mu}{}_{\\lambda}\\omega^{\\lambda\\nu} + \\pi_{Q}^{\\nu}{}_{\\lambda}\\omega^{\\lambda\\mu},\\qquad \\omega^{\\mu\\nu} = \\tfrac{1}{2}(\\Delta^{\\mu\\lambda}\\nabla_{\\lambda}u^{\\nu} - \\Delta^{\\nu\\lambda}\\nabla_{\\lambda}u^{\\mu})"}},
    {"text", "so $w_{\\omega} = 1$ is the covariant row, which is exported; $w_{\\omega} = 0$ removes the corotation coupling. The five systems (S1, S2 and every 0+1D system are derived for the checks; S3, S4, S5 in 1+1D and 2+1D are exported, with $\\delta_{\\Pi\\Pi} = 0$):"},
    {"eq", {"\\text{S1}:\\ \\nabla_{\\mu}T_{m}^{\\mu\\nu} = 0,\\ \\pi = \\Pi = 0\\qquad \\text{S2}:\\ \\text{S1} + \\text{the }\\pi, \\Pi\\text{ equations}",
      "\\text{S3}:\\ \\text{charm number} + \\text{first moment, background }T, u\\qquad \\text{S4}:\\ \\text{S3} + \\text{second moment, background }T, u",
      "\\text{S5}:\\ \\nabla_{\\mu}(T_{m} + T_{Q})^{\\mu\\nu} = 0,\\ \\text{the }\\pi, \\Pi\\text{ equations, the charm rows of S4; all fields evolved}",
      "\\text{S6}:\\ \\nabla_{\\mu}(T_{m} + T_{Q} + \\delta M)^{\\mu\\nu} = 0,\\ \\frac{D_{s}}{T}\\Delta^{\\mu}{}_{\\lambda}\\nabla_{\\nu}(T_{Q} + \\delta M)^{\\nu\\lambda} + \\nu^{\\mu} = 0,\\ \\text{the rest as S5}"}},
    {"text", "Conformal: $P = a_{\\mathrm{SB}}T^{4}$, no bulk. Source on the left, $A_t\\,\\partial_{\\tau}\\phi + A_r\\,\\partial_{r}\\phi + S = 0$ in 1+1D and $A_t\\,\\partial_{\\tau}\\phi + A_x\\,\\partial_{x}\\phi + A_y\\,\\partial_{y}\\phi + S = 0$ in 2+1D."},
    {"subsub", "Fields"},
    {"text", "Medium, charm first moment, charm second moment, chart by chart:"},
    {"eq", {"\\text{0+1D}:\\ " <> fl[stateBJ, "med"] <> ",\\ " <> fl[stateBJ, "ch"] <> ",\\ " <> fl[stateBJ, "m2"],
      "\\text{1+1D}:\\ " <> fl[stateCYL, "med"] <> ",\\ " <> fl[stateCYL, "ch"] <> ",\\ " <> fl[stateCYL, "m2"],
      "\\text{2+1D}:\\ " <> fl[stateCAR, "med"] <> ",\\ " <> fl[stateCAR, "ch"],
      fl[stateCAR, "m2"]}},
    {"text", "Rows, in the same order: the fluid rows $\\nabla_{\\mu}T^{\\mu\\tau}$, $\\nabla_{\\mu}T^{\\mu i}$ ($T_{m} + T_{Q}$ in S5); the shear rows $E^{\\eta\\eta}$ (0+1D), $E^{\\phi\\phi}$, $E^{\\eta\\eta}$ (1+1D), $E^{yy}$, $E^{\\eta\\eta}$, $E^{xy}$ (2+1D); the bulk row; $\\nabla_{\\mu}N^{\\mu}$; the current components; the traceless second-moment rows, components $\\eta\\eta$ (0+1D), $\\phi\\phi$, $\\eta\\eta$ (1+1D), $yy$, $\\eta\\eta$, $xy$ (2+1D), as the medium's shear rows; the trace row."},
    {"text", "In 0+1D the tensors are diagonal and $\\nu = 0$:"},
    {"eq", {"\\pi^{\\mu\\nu} = \\mathrm{diag}" <> rdTuple[rdTeX /@ Diagonal[toJ[stateBJ["piM"], stateBJ]]] <>
      ",\\qquad \\pi_{Q}^{\\mu\\nu} = \\mathrm{diag}" <> rdTuple[rdTeX /@ Diagonal[toJ[stateBJ["piQ"], stateBJ]]]}},
    {"text", "In 1+1D the medium $\\pi^{\\mu\\nu}$ is as in the viscous section, and so is the charm stress: $\\pi_{Q}{}^{\\phi}{}_{\\phi}$ and $\\pi_{Q}{}^{\\eta}{}_{\\eta}$ are independent, the radial component follows from the trace. With the orthonormal triad transverse to $u$,"},
    {"eq", {"\\pi_{Q}^{\\mu\\nu} = -(\\pi_{Q}{}^{\\phi}{}_{\\phi} + \\pi_{Q}{}^{\\eta}{}_{\\eta})\\,\\ell^{\\mu}\\ell^{\\nu} + \\pi_{Q}{}^{\\phi}{}_{\\phi}\\,\\hat\\phi^{\\mu}\\hat\\phi^{\\nu} + \\pi_{Q}{}^{\\eta}{}_{\\eta}\\,\\hat\\eta^{\\mu}\\hat\\eta^{\\nu}",
      "\\ell^{\\mu} = " <> rdTuple[rdTeX /@ (toJ[stateCYL["triad"][[1]], stateCYL] /. Sqrt[1 + ur^2] -> utau)] <>
      ",\\quad \\hat\\phi^{\\mu} = " <> rdTuple[rdTeX /@ stateCYL["triad"][[2]]] <> ",\\quad \\hat\\eta^{\\mu} = " <>
      rdTuple[rdTeX /@ stateCYL["triad"][[3]]]}},
    {"text", "In 2+1D the medium $\\pi^{\\mu\\nu}$ is as in the viscous section, and so is the charm stress: $\\pi_{Q}^{yy}$, $\\pi_{Q}{}^{\\eta}{}_{\\eta}$, $\\pi_{Q}^{xy}$ are independent, the $\\tau$ row follows from $u_{\\mu}\\pi_{Q}^{\\mu\\nu} = 0$ and"},
    {"eq", {"\\pi_{Q}^{xx} = " <> rdTeX[TimeConstrained[Simplify[toJ[stateCAR["qxx"], stateCAR]], 20,
        toJ[stateCAR["qxx"], stateCAR]] /. rdUtTW2[[1]]]}},
    {"sub", "Steps"},
    {"subsub", "S3: charm first moment"},
    {"eq", {rdGammaEq[st1["geo"]]}},
    {"entries", Append[Table[rdPrs["N^{" <> x1[[m]] <> "}", B1["Nvec"][[m]], s1], {m, 4}],
      rdPrs["E_{\\alpha} = \\nabla_{\\mu}N^{\\mu}", R1["ch"][[1]], s1]]},
    {"entries", rdComps[B1["TQ"], "T_{Q}", x1, s1]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivForm[st1["geo"], st1["cs"], B1["TQ"], {1, 2}, "T_{Q}", divT2]},
    {"entries", {rdPrs["\\nabla_{\\mu}T_{Q}^{\\mu\\tau}", B1["consQ"][[1]], s1], rdPrs["\\nabla_{\\mu}T_{Q}^{\\mu r}", B1["consQ"][[2]], s1],
      rdPrs["\\Delta^{r}{}_{\\tau}", Dmx[[2, 1]], s1], rdPrs["\\Delta^{r}{}_{r}", Dmx[[2, 2]], s1]}},
    {"eq", {"E_{\\nu} = \\frac{D_{s}}{T}\\bigl(\\Delta^{r}{}_{\\tau}\\nabla_{\\mu}T_{Q}^{\\mu\\tau} + \\Delta^{r}{}_{r}\\nabla_{\\mu}T_{Q}^{\\mu r}\\bigr) + \\nu^{r}"}},
    {"entries", {rdPrs["E_{\\nu}", R1["ch"][[2]], s1]}},
    {"subsub", "S4: charm second moment"},
    {"text", "$\\mathcal{A}$ enters $S^{\\lambda\\alpha\\beta}$ only through $u^{\\lambda}u^{\\alpha}\\nu^{\\beta}$: after the projection its derivative has a vanishing coefficient (G-aQ)."},
    {"eq", {"\\nabla_{\\lambda}S^{\\lambda\\alpha\\beta} = \\partial_{\\lambda}S^{\\lambda\\alpha\\beta} + \\Gamma^{\\lambda}{}_{\\lambda\\sigma}S^{\\sigma\\alpha\\beta} + \\Gamma^{\\alpha}{}_{\\lambda\\sigma}S^{\\lambda\\sigma\\beta} + \\Gamma^{\\beta}{}_{\\lambda\\sigma}S^{\\lambda\\alpha\\sigma}"}},
    {"entries", {rdPrs["\\mathcal{E}^{\\phi\\phi}", lhs1[[3, 3]], s1], rdPrs["\\mathcal{E}^{\\eta\\eta}", lhs1[[4, 4]], s1],
      rdPrs["\\mathcal{E}^{\\lambda}{}_{\\lambda}", tr1, s1]}},
    {"eq", {"E_{\\phi} = \\frac{D_{s}}{2T}\\Bigl(r^{2}\\mathcal{E}^{\\phi\\phi} - \\tfrac{1}{3}\\mathcal{E}^{\\lambda}{}_{\\lambda}\\Bigr) + \\pi_{Q}{}^{\\phi}{}_{\\phi},\\quad E_{\\eta} = \\frac{D_{s}}{2T}\\Bigl(\\tau^{2}\\mathcal{E}^{\\eta\\eta} - \\tfrac{1}{3}\\mathcal{E}^{\\lambda}{}_{\\lambda}\\Bigr) + \\pi_{Q}{}^{\\eta}{}_{\\eta},\\quad E_{\\Pi_{Q}} = \\frac{D_{s}}{6T}\\mathcal{E}^{\\lambda}{}_{\\lambda} + \\Pi_{Q}"}},
    {"subsub", "S5: two-way"},
    {"entries", rdComps[B1["Tm"], "T_{m}", x1, s1]},
    {"eq", {"E_{\\tau} = \\nabla_{\\mu}T_{m}^{\\mu\\tau} + \\nabla_{\\mu}T_{Q}^{\\mu\\tau},\\qquad E_{r} = \\nabla_{\\mu}T_{m}^{\\mu r} + \\nabla_{\\mu}T_{Q}^{\\mu r}"}},
    {"text", "The shear and bulk rows as in the viscous section, the charm rows as S4."},
    {"entries", {rdPrs["E_{\\tau}", sys1["S5"]["eqs"][[1]], s1]}},
    {"subsub", "S6: fully coupled"},
    {"entries", rdComps[B1["dM"], "\\delta M", x1, s1]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivForm[st1["geo"], st1["cs"], B1["dM"], {1, 2}, "\\delta M", divT2]},
    {"entries", {rdPrs["\\nabla_{\\mu}\\delta M^{\\mu\\tau}", B1["consD"][[1]], s1], rdPrs["\\nabla_{\\mu}\\delta M^{\\mu r}", B1["consD"][[2]], s1]}},
    {"eq", {"E_{\\tau}^{\\mathrm{S6}} = E_{\\tau}^{\\mathrm{S5}} + \\nabla_{\\mu}\\delta M^{\\mu\\tau},\\quad E_{r}^{\\mathrm{S6}} = E_{r}^{\\mathrm{S5}} + \\nabla_{\\mu}\\delta M^{\\mu r},\\quad E_{\\nu}^{\\mathrm{S6}} = E_{\\nu} + c_{M}\\bigl(\\Delta^{r}{}_{\\tau}\\nabla_{\\mu}\\delta M^{\\mu\\tau} + \\Delta^{r}{}_{r}\\nabla_{\\mu}\\delta M^{\\mu r}\\bigr),\\quad c_{M} = \\frac{D_{s}}{T}"}},
    {"text", "Every other row as S5."},
    {"subsub", "From the rows to the matrices"},
    {"eq", {rdExtractEq["\\tau, r"]}},
    {"subsub", "2+1D, S3"},
    {"eq", {rdGammaEq[stK["geo"]]}},
    {"entries", Append[Table[rdPrs["N^{" <> xK[[m]] <> "}", BK["Nvec"][[m]], sK], {m, 4}],
      rdPrs["E_{\\alpha} = \\nabla_{\\mu}N^{\\mu}", RK["ch"][[1]], sK]]},
    {"entries", rdComps[BK["TQ"], "T_{Q}", xK, sK]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivForm[stK["geo"], stK["cs"], BK["TQ"], {1, 2, 3}, "T_{Q}", divT2]},
    {"entries", Table[rdPrs["\\nabla_{\\mu}T_{Q}^{\\mu " <> xK[[m]] <> "}", BK["consQ"][[m]], sK], {m, 3}]},
    {"entries", Flatten[Table[rdPrs["\\Delta^{" <> xK[[m]] <> "}{}_{" <> xK[[l]] <> "}", DmxK[[m, l]], sK], {m, {2, 3}}, {l, 3}], 1]},
    {"eq", {"E_{\\nu^{i}} = \\frac{D_{s}}{T}\\Delta^{i}{}_{\\lambda}\\nabla_{\\mu}T_{Q}^{\\mu\\lambda} + \\nu^{i},\\qquad i = x, y"}},
    {"entries", {rdPrs["E_{\\nu^{x}}", RK["ch"][[2]], sK], rdPrs["E_{\\nu^{y}}", RK["ch"][[3]], sK]}},
    {"subsub", "2+1D, S4"},
    {"entries", Flatten[Table[If[l <= m <= k, rdPrs["S^{" <> xK[[l]] <> xK[[m]] <> xK[[k]] <> "}", BK["S"][[l, m, k]], sK], Nothing], {l, 4}, {m, 4}, {k, 4}], 2]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivS[stK["geo"], stK["cs"], BK["S"], {{2, 2}, {2, 3}, {3, 3}, {4, 4}}]},
    {"eq", {"\\mathcal{E}^{ab} = \\Delta^{a}{}_{\\alpha}\\Delta^{b}{}_{\\beta}\\nabla_{\\lambda}S^{\\lambda\\alpha\\beta},\\qquad \\mathcal{E}^{\\lambda}{}_{\\lambda} = \\Delta_{ab}\\mathcal{E}^{ab}"}},
    {"eq", {"E_{ab} = \\frac{D_{s}}{2T}\\Bigl(\\mathcal{E}^{ab} - \\tfrac{1}{3}\\Delta^{ab}\\mathcal{E}^{\\lambda}{}_{\\lambda}\\Bigr) + \\pi_{Q}^{ab},\\quad ab = yy,\\ xy,\\qquad E_{\\eta} = \\tau^{2}\\frac{D_{s}}{2T}\\Bigl(\\mathcal{E}^{\\eta\\eta} - \\tfrac{1}{3}\\Delta^{\\eta\\eta}\\mathcal{E}^{\\lambda}{}_{\\lambda}\\Bigr) + \\pi_{Q}{}^{\\eta}{}_{\\eta}",
      "E_{\\Pi_{Q}} = \\frac{D_{s}}{6T}\\mathcal{E}^{\\lambda}{}_{\\lambda} + \\Pi_{Q}"}},
    {"subsub", "2+1D, S5"},
    {"entries", rdComps[BK["Tm"], "T_{m}", xK, sK]},
    {"eq", {"E_{\\tau} = \\nabla_{\\mu}(T_{m} + T_{Q})^{\\mu\\tau},\\quad E_{x} = \\nabla_{\\mu}(T_{m} + T_{Q})^{\\mu x},\\quad E_{y} = \\nabla_{\\mu}(T_{m} + T_{Q})^{\\mu y}", rdExtractEq["\\tau, x, y"]}},
    {"entries", {rdPrs["E_{\\tau}", sysK["S5"]["eqs"][[1]], sK]}},
    {"subsub", "2+1D, S6"},
    {"entries", rdComps[BK["dM"], "\\delta M", xK, sK]},
    {"entries", {#[[1]], rdPrep[#[[2]], {}, rdNoUt, "div form"]} & /@ rdDivForm[stK["geo"], stK["cs"], BK["dM"], {1, 2, 3}, "\\delta M", divT2]},
    {"entries", Table[rdPrs["\\nabla_{\\mu}\\delta M^{\\mu " <> xK[[m]] <> "}", BK["consD"][[m]], sK], {m, 3}]},
    {"eq", {"E_{a}^{\\mathrm{S6}} = E_{a}^{\\mathrm{S5}} + \\nabla_{\\mu}\\delta M^{\\mu a},\\quad a = \\tau, x, y,\\qquad E_{\\nu^{i}}^{\\mathrm{S6}} = E_{\\nu^{i}} + c_{M}\\Delta^{i}{}_{\\lambda}\\nabla_{\\mu}\\delta M^{\\mu\\lambda}"}},
    {"text", "Conformal:"},
    {"eq", {"P = a_{\\mathrm{SB}}T^{4},\\quad P' = 4a_{\\mathrm{SB}}T^{3},\\quad P'' = 12a_{\\mathrm{SB}}T^{2},\\quad \\Pi = \\zeta = 0,\\quad\\text{drop the }\\Pi\\text{ row and column}"}},
    {"sub", "Result"},
    {"text", "The 1+1D charm systems and the 1+1D conformal two-way system in full; the generic 1+1D two-way system and every 2+1D system as the pattern of nonzero entries ($\\bullet$), rows and columns in the field order given, with their entries in the appendix."}};
  Do[Module[{fns = Select[exported[[All, 1]], renderSys[#]["chart"] === nm &]},
      AppendTo[items, {"sub", chartLab[nm]}];
      Do[items = Join[items, fnItems[fn, renderSys[fn], full[renderSys[fn]]]], {fn, fns}]],
    {nm, Select[{"1p1d", "2p1d"}, KeyExistsQ[results, #] &]}];
  items = Join[items, {
    {"sub", "Checks"},
    {"checks", gateLog, {}},
    {"text", "Not a check: the characteristic speeds of S5 are computed and reported (\\texttt{two\\_way\\_char\\_speeds.tsv} in the scratch directory)."},
    {"source", {FileNameJoin[{DERIV, "derive_two_way.wls"}], FileNameJoin[{DERIV, "derive_two_way_gates.wls"}],
      FileNameJoin[{DERIV, "derive_two_way_export.wls"}]}}}];
  rdWrite["two_way", "Medium, charm and the two-way coupled system", items];
  appItems = Join[{{"sec", "The entries of the remaining two-way and charm matrices"},
      {"text", "As computed by \\texttt{deriv/derive\\_two\\_way.wls}: the entries of its generated blocks, each compared with the derived expression (E-same). Field names and $u^{\\tau}$ as in the section on the two-way system."}},
    Flatten[Table[Module[{fns = Select[exported[[All, 1]], renderSys[#]["chart"] === nm && !full[renderSys[#]] &]},
        If[fns === {}, {}, Join[{{"sub", chartLab[nm]}}, Flatten[fnItems[#, renderSys[#], True] & /@ fns, 1]]]],
      {nm, Select[{"1p1d", "2p1d"}, KeyExistsQ[results, #] &]}], 1]];
  rdWrite["two_way_full", "The entries of the remaining two-way and charm matrices", appItems]];
rdUtTW2 = {{Power[1 + ux^2 + uy^2, pw_] :> utau^(2 pw)}, {utau -> Sqrt[1 + ux^2 + uy^2]}};

(* ============================================================================
   derive_1d_matrices.wls
   ============================================================================ *)
renderOneD[] := Module[{rs, uts, pr, prs, labs, items, ent, x, acc, Dnu, nda, shE, elimItems, sysItems, stL},
  rdSymTeX = Join[rdBaseTeX, <|"R" -> "r", "uT" -> "T", "uR" -> "u^{r}", "pFF" -> "\\pi^{\\phi}{}_{\\phi}", "pEE" -> "\\pi^{\\eta}{}_{\\eta}",
    "pBB" -> "\\Pi", "aL" -> "\\alpha", "nuR" -> "\\nu^{r}", "nF" -> "n", "gM" -> "\\gamma", "pres" -> "p", "dtp" -> "p'(T)", "dtdtp" -> "p''(T)",
    "visc" -> "\\eta", "zeta" -> "\\zeta", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}", "n" -> "n", "dtn" -> "(\\partial_{T} n)",
    "dmn" -> "(\\partial_{\\alpha} n)", "tauDiff" -> "\\tau_{n}", "kappa" -> "\\kappa", "dtT" -> rdDeriv["\\tau", "T"], "drT" -> rdDeriv["r", "T"],
    "dtur" -> rdDeriv["\\tau", "u^{r}"], "drur" -> rdDeriv["r", "u^{r}"]|>];
  rs = Join[rdFieldRules[{{TT, "uT", "T"}, {ur, "uR", "u^{r}"}, {pff, "pFF", "\\pi^{\\phi}{}_{\\phi}"}, {pee, "pEE", "\\pi^{\\eta}{}_{\\eta}"},
      {PiB, "pBB", "\\Pi"}, {al, "aL", "\\alpha"}, {nur, "nuR", "\\nu^{r}"}, {nf, "nF", "n"}, {gm, "gM", "\\gamma"}}, {tau, rr}], subD];
  uts = {{Power[1 + uR^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + uR^2 :> utau^2}, {utau -> Sqrt[1 + uR^2]}};
  pr = rdPrep[#, {}, uts, "1+1D"] &;
  prs = rdPrep[#, rs, uts, "steps 1+1D"] &;
  x = rdTeX /@ CYL["x"];
  acc = DcovV[CYL, uC, uC];
  Dnu = DcovV[CYL, uC, {ur[tau, rr] nur[tau, rr]/utC, nur[tau, rr], 0, 0}];
  nda = Sum[(CYL["g"] . acc)[[m]] {ur[tau, rr] nur[tau, rr]/utC, nur[tau, rr], 0, 0}[[m]], {m, 4}];
  (* the rows, then the elimination, for one kernel *)
  ent[mat_, lab_, rowsel_, cols_] := Flatten[Table[rdPrs["(" <> lab <> ")_{" <> ToString[i] <> ToString[j] <> "}", mat[[i, j]], pr], {i, rowsel}, {j, cols}], 1];
  elimItems[der_, sh_, nm_, fn_, file_] := Module[{n = Length[der[[1]]], c, L, entL, entS, rules, defs},
    c = Range[nm + 1, n];
    L = elimL[der, nm];
    entL = Flatten[Table[rdPrs["L_{" <> ToString[i] <> ToString[j] <> "}", L[[i - nm, j]], pr], {i, c}, {j, nm}], 1];
    entS = Join[ent[sh[[1]], "A_t", c, Range[n]], ent[sh[[2]], "A_r", c, Range[n]],
      Table[rdPrs["S_{" <> ToString[i] <> "}", sh[[3, i]], pr], {i, c}]];
    {rules, defs} = rdDenNames[Join[entL, entS]];
    {{"subsub", "\\texttt{" <> rdEscId[fn] <> "} in \\texttt{src/Matrix/" <> rdEscId[file] <> "}"},
     {"text", "Derived, charm rows:"},
     {"entries", Join[ent[der[[1]], "A_t", c, Range[n]], ent[der[[2]], "A_r", c, Range[n]],
        Table[rdPrs["S_{" <> ToString[i] <> "}", der[[3, i]], pr], {i, c}]]},
     {"text", "With $c$ the charm rows and $m$ the medium rows, $\\mathrm{row}_{c}\\to\\mathrm{row}_{c} - L\\,\\mathrm{row}_{m}$, $L = A_t^{cm}(A_t^{mm})^{-1}$, removes $\\partial_{\\tau}$ of the medium from the charm rows:"},
     {"entries", defs},
     {"entries", rdDenApply[entL, rules]},
     {"text", "The file, charm rows (M-V7e, M-G7e, M-I3ne, M-I3ae):"},
     {"entries", rdDenApply[entS, rules]}}];
  rdRebuild[{der7[[1]], der7[[2]]}, der7[[3]], vars7, {tau, rr}, eqs7 /. tauDiff -> tauDiff, rs, "steps: A dphi + S == E, 7 fields"];
  items = {
    {"sec", "The hand-written 1+1D matrices"},
    {"text", "Script: \\texttt{deriv/derive\\_1d\\_matrices.wls}"},
    {"text", "Writes nothing into \\texttt{src/Matrix}: it compares the derivation with the hand-written \\texttt{.jl} files, read from the file."},
    {"sub", "Setup"},
    {"eq", {rdDiag[CYL] <> ",\\qquad u^{\\mu} = (u^{\\tau}, u^{r}, 0, 0)"}},
    {"text", "Medium as in the viscous section. Charm:"},
    {"eq", {"\\nabla_{\\mu}(n\\,u^{\\mu} + \\nu^{\\mu}) = 0,\\qquad \\tau_{n}\\,\\Delta^{\\mu}{}_{\\beta}D\\nu^{\\beta} + \\nu^{\\mu} + \\kappa\\,\\Delta^{\\mu\\beta}\\partial_{\\beta}\\alpha = 0,\\qquad \\nu^{\\tau} = \\frac{u^{r}\\nu^{r}}{u^{\\tau}}"}},
    {"text", "Rows $\\nabla_{\\mu}N^{\\mu}$ and the $r$ component of the current equation. Variants:"},
    {"eq", {"\\text{MIS}:\\ n = n(T, \\alpha),\\ \\delta_{\\pi\\pi} = \\tfrac{4}{3}\\tau_{\\pi},\\quad \\phi = (T, u^{r}, \\pi^{\\phi}{}_{\\phi}, \\pi^{\\eta}{}_{\\eta}, \\Pi, \\alpha, \\nu^{r})",
      "\\text{gamma}:\\ n = \\gamma\\,n(T),\\ \\kappa\\partial\\alpha\\to\\kappa\\partial\\gamma,\\ \\delta_{\\pi\\pi} = \\tfrac{4}{3}\\tau_{\\pi},\\quad \\phi = (T, u^{r}, \\pi^{\\phi}{}_{\\phi}, \\pi^{\\eta}{}_{\\eta}, \\Pi, \\gamma, \\nu^{r})",
      "\\text{ideal}:\\ \\pi = \\Pi = 0,\\ \\nu = 0,\\quad \\phi = (T, u^{r}, n)\\ \\text{or}\\ (T, u^{r}, \\alpha)"}},
    {"sub", "Steps"},
    {"eq", {rdGammaEq[CYL]}},
    {"entries", Join[Table[rdPrs["a^{" <> x[[m]] <> "}", acc[[m]], prs], {m, 2}], {rdPrs["\\nu_{\\lambda}a^{\\lambda}", nda, prs]},
      Table[rdPrs["D\\nu^{" <> x[[m]] <> "}", Dnu[[m]], prs], {m, 2}]]},
    {"eq", {"E_{\\nu} = \\tau_{n}\\bigl(D\\nu^{r} - u^{r}\\nu_{\\lambda}a^{\\lambda}\\bigr) + \\nu^{r} + \\kappa\\bigl(\\Delta^{r\\tau}\\partial_{\\tau}\\alpha + \\Delta^{rr}\\partial_{r}\\alpha\\bigr),\\qquad \\Delta^{r\\tau} = u^{r}u^{\\tau},\\ \\Delta^{rr} = 1 + (u^{r})^{2}"}},
    {"entries", {rdPrs["E_{\\alpha}", eqs7[[6]], prs], rdPrs["E_{\\nu}", eqs7[[7]], prs]}},
    {"eq", {rdExtractEq["\\tau, r"]}},
    {"sub", "Result"},
    {"subsub", "\\texttt{one\\_d\\_viscous\\_matrix\\_fugacity\\_BG\\_only} in \\texttt{src/Matrix/1d\\_viscous\\_HQ\\_cilindrical\\_BG.jl}"},
    {"text", "Background $T$, $u^{r}$; $\\phi = (\\alpha, \\nu^{r})$; the file multiplies every row by $(u^{\\tau})^{3}r\\tau$ (M-BG). As derived:"},
    {"entries", Join[rdMat[derBG[[1]], "A_t", pr], rdMat[derBG[[2]], "A_r", pr], rdVec[derBG[[3]], "S", pr]]},
    Splice[elimItems[der7, sh7, 5, "one_d_viscous_matrix_fugacity", "1d_viscous_HQ_cilindrical_fugacity.jl"]],
    {"text", "Rows 1--5 are the viscous section's system at $\\delta_{\\pi\\pi} = \\tfrac{4}{3}\\tau_{\\pi}$ (M-V7m)."},
    Splice[elimItems[derG, shG, 5, "one_d_viscous_matrix_gamma", "1d_viscous_HQ_cilindrical_gamma.jl"]],
    {"text", "$\\partial_{T}n$ above is $\\gamma$-free: $n$ and $\\partial_{T}n$ are those of $n(T)$ at $\\gamma = 1$."},
    {"subsub", "Ideal medium, rows 1--2"},
    {"entries", Join[ent[derIn[[1]], "A_t", {1, 2}, {1, 2}], ent[derIn[[2]], "A_r", {1, 2}, {1, 2}],
        Table[rdPrs["S_{" <> ToString[i] <> "}", derIn[[3, i]], pr], {i, 2}]]},
    Splice[elimItems[derIn, shIn, 2, "one_d_ideal_matrix_ruwen4", "1d_ideal_HQ_cilindrical_density.jl"]],
    Splice[elimItems[derIa, shIa, 2, "one_d_ideal_matrix_ruwen5", "1d_ideal_HQ_cilindrical_fugacity.jl"]],
    {"sub", "Checks"},
    {"checks", gateLog, {}},
    {"source", {FileNameJoin[{DERIV, "derive_1d_matrices.wls"}]}}};
  rdWrite["one_d", "The hand-written 1+1D matrices", items]];

(* ============================================================================
   derive_upstream.wls
   ============================================================================ *)
renderUpstream[] := Module[{items, pr, pc, wrongItems, names, prC},
  rdSymTeX = Join[rdBaseTeX, <|"T" -> "T", "ux" -> "u^{x}", "uy" -> "u^{y}", "pyy" -> "\\pi^{yy}", "pee" -> "\\pi^{\\eta}{}_{\\eta}", "pxy" -> "\\pi^{xy}",
    "PB" -> "\\Pi", "p0" -> "p", "p1" -> "p'(T)", "p2" -> "p''(T)", "zeta" -> "\\zeta", "visc" -> "\\eta", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}",
    "al" -> "\\alpha", "nx" -> "\\nu^{x}", "ny" -> "\\nu^{y}", "n" -> "n", "nT" -> "(\\partial_{T}n)", "na" -> "(\\partial_{\\alpha}n)", "taun" -> "\\tau_{n}", "Ds" -> "D_{s}",
    "u" -> "u^{r}", "al0" -> "\\alpha", "Y" -> "n", "R" -> "r", "pT" -> "(\\partial_{T}P_{h})", "pa" -> "(\\partial_{\\alpha}P_{h})", "pTT" -> "(\\partial^{2}_{T}P_{h})",
    "pTa" -> "(\\partial_{T}\\partial_{\\alpha}P_{h})", "paa" -> "(\\partial^{2}_{\\alpha}P_{h})"|>];
  rdSymTeX["uR"] = "u^{r}";
  pr = rdPrep[#, {}, rdUtCar /. {ux -> ux, uy -> uy}, "2+1D"] &;
  pc = rdPrep[#, {u -> uR}, {{Power[1 + uR^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + uR^2 :> utau^2}, {utau -> Sqrt[1 + uR^2]}}, "1+1D"] &;
  names = <|"At" -> "A_t", "Ax" -> "A_x", "Ay" -> "A_y"|>;
  wrongItems = Flatten[Function[w, With[{k = w[[1]], ij = w[[2]], lab = names[w[[1]]] <> ")_{" <> ToString[w[[2, 1]]] <> ToString[w[[2, 2]]] <> "}"},
      {{"(" <> lab <> "^{\\mathrm{legacy}}", leg[k][[ij[[1]], ij[[2]]]]}, {"(" <> lab <> "^{\\mathrm{derived}}", der[k][[ij[[1]], ij[[2]]]]},
       {"\\delta(" <> lab, dLeg[k][[ij[[1]], ij[[2]]]]}}]] /@ wrong, 1];
  items = {
    {"sec", "The matrices no production path relies on"},
    {"text", "Script: \\texttt{deriv/derive\\_upstream.wls}"},
    {"text", "Writes nothing into \\texttt{src/Matrix}."},
    {"sub", "\\texttt{two\\_d\\_viscous\\_matrix} (\\texttt{legacy = true})"},
    {"text", "Against the first section's matrix at $\\delta_{\\pi\\pi} = 0$, $\\phi = (T, u^{x}, u^{y}, \\pi^{yy}, \\pi^{\\eta}{}_{\\eta}, \\pi^{xy}, \\Pi)$. The source and 137 entries agree. The ten that differ, with $\\delta$ = legacy $-$ derived:"},
    {"entries", {#[[1]], pr[#[[2]]]} & /@ wrongItems},
    {"sub", "\\texttt{two\\_d\\_viscous\\_HQ\\_matrix}"},
    {"text", "$\\phi = (T, u^{x}, u^{y}, \\pi^{yy}, \\pi^{\\eta}{}_{\\eta}, \\pi^{xy}, \\Pi, \\alpha, \\nu^{x}, \\nu^{y})$. Rows 1--7: the derived matrix at $\\delta_{\\pi\\pi} = 0$. Rows 8--10:"},
    {"eq", {"\\nabla_{\\mu}(n u^{\\mu} + \\nu^{\\mu}) = 0,\\qquad \\tau_{n}\\Delta^{\\mu}{}_{\\beta}D\\nu^{\\beta} + \\nu^{\\mu} + \\kappa\\Delta^{\\mu\\beta}\\partial_{\\beta}\\alpha = 0,\\qquad \\kappa = D_{s}n"}},
    {"text", "The steps are those of the second section. With the medium as fields, in columns 1, 2, 3, 8, 9, 10:"},
    {"entries", Join[
      Flatten[Table[{"(A_{" <> {"t", "x", "y"}[[k]] <> "})_{" <> ToString[7 + i] <> "," <> ToString[colsH[[j]]] <> "}", pr[AsH[[k, i, j]] /. subH]}, {k, 3}, {i, 3}, {j, 6}], 2],
      Table[{"S_{" <> ToString[7 + i] <> "}", pr[srcH[[i]] /. subH]}, {i, 3}]] /. {_String, 0} -> Nothing},
    {"text", "Its wrapper is \\texttt{matrix2d\\_visc\\_HQ!} (D-HQup)."},
    {"sub", "\\texttt{one\\_d\\_ideal\\_matrix\\_ruwen7}"},
    {"eq", {rdDiag[CYL] <> ",\\qquad u^{\\mu} = (u^{\\tau}, u^{r}, 0, 0)",
      "P_{\\mathrm{tot}} = P(T) + P_{h}(T, \\alpha),\\qquad T^{\\mu\\nu} = (e + P_{\\mathrm{tot}})u^{\\mu}u^{\\nu} + P_{\\mathrm{tot}}g^{\\mu\\nu},\\qquad \\nabla_{\\mu}T^{\\mu\\nu} = 0,\\qquad \\nabla_{\\mu}(n u^{\\mu}) = 0",
      "e + P_{\\mathrm{tot}} = Ts + \\mu n = T\\partial_{T}P_{\\mathrm{tot}}\\bigl|_{\\mu} + \\mu\\partial_{\\mu}P_{\\mathrm{tot}}\\bigl|_{T} = T\\partial_{T}P_{\\mathrm{tot}}\\bigl|_{\\alpha},\\qquad \\mu = \\alpha T"}},
    {"text", "$\\phi = (T, u^{r}, \\alpha)$, the $T^{\\mu\\nu}$ as written, the same extraction. The file, raw (U-P7):"},
    {"entries", Join[rdMat[At7, "A_t", pc], rdMat[Ar7, "A_r", pc], rdVec[S7, "S", pc]]},
    {"text", "Its wrapper is \\texttt{matrix1d\\_ideal\\_HQ\\_br\\_p\\_fugacity!} (D-r7). For a Boltzmann gas of unit-charge species $P_{h} = nT$:"},
    {"eq", {"\\partial_{T}P_{h} = n + T\\partial_{T}n,\\quad \\partial_{\\alpha}P_{h} = T\\partial_{\\alpha}n,\\quad \\partial^{2}_{T}P_{h} = 2\\partial_{T}n + T\\partial^{2}_{T}n,\\quad \\partial_{T}\\partial_{\\alpha}P_{h} = \\partial_{\\alpha}n + T\\partial_{T}\\partial_{\\alpha}n,\\quad \\partial^{2}_{\\alpha}P_{h} = T\\partial^{2}_{\\alpha}n"}},
    {"sub", "\\texttt{one\\_d\\_ideal\\_matrix\\_ruwen6}"},
    {"text", "The same system in $\\phi = (T, u^{r}, n)$, with $n = \\partial_{\\alpha}P_{h}/T$ and $\\alpha = \\alpha(T, n)$:"},
    {"eq", {"\\partial_{n}\\alpha\\bigl|_{T} = \\frac{T}{\\partial^{2}_{\\alpha}P_{h}},\\qquad \\partial_{T}\\alpha\\bigl|_{n} = \\frac{\\partial_{\\alpha}P_{h}/T - \\partial_{T}\\partial_{\\alpha}P_{h}}{\\partial^{2}_{\\alpha}P_{h}},\\qquad \\partial_{k}\\alpha = \\partial_{T}\\alpha\\,\\partial_{k}T + \\partial_{n}\\alpha\\,\\partial_{k}n"}},
    {"text", "Derived:"},
    {"entries", Join[rdMat[At6, "A_t", pc], rdMat[Ar6, "A_r", pc], rdVec[S6, "S", pc]]},
    {"text", "The file agrees, raw (U-P6-cols, U-P6-row3, U-P6-c3). Its wrapper takes $P_{h} = nT$ from \\texttt{charm\\_pressure} and $\\alpha$ from \\texttt{charm\\_alpha} (D-r6)."},
    {"sub", "Checks"},
    {"checks", gateLog, {}},
    {"source", {FileNameJoin[{DERIV, "derive_upstream.wls"}]}}};
  rdWrite["upstream", "The matrices no production path relies on", items]];

(* ============================================================================
   derive_drivers.wls
   ============================================================================ *)
renderDrivers[] := Module[{items},
  rdSymTeX = rdBaseTeX;
  items = {
    {"sec", "From the matrices to the solver"},
    {"text", "Script: \\texttt{deriv/derive\\_drivers.wls}, which runs \\texttt{deriv/audit\\_matrices.jl}"},
    {"sub", "The solved form"},
    {"eq", {"\\partial_{\\tau}\\phi = -A_t^{-1}A_i\\,\\partial_{i}\\phi - A_t^{-1}S"}},
    {"text", "Each driver evaluates the EOS and the transport at the state, calls a kernel, and returns $A_t^{-1}A_i$, $A_t^{-1}S$; the solver upwinds $A_t^{-1}A_i\\partial_{i}\\phi$. $\\eta$, $\\tau_{\\pi}$, $\\zeta$, $\\tau_{\\Pi}$ are inputs, not derived. $\\delta_{\\pi\\pi} = \\tfrac{4}{3}\\tau_{\\pi}$ unless switched off."},
    {"sub", "$\\kappa$ and $\\tau_{n}$"},
    {"text", "From the consistent first moment of the section on the two-way system, at rest ($u^{r} = 0$, $\\nu^{r} = 0$ in the coefficients):"},
    {"eq", {"\\frac{D_{s}}{T}\\Delta^{r}{}_{\\lambda}\\nabla_{\\mu}T_{Q}^{\\mu\\lambda} + \\nu^{r}\\ \\to\\ \\frac{D_{s}h}{T}\\partial_{\\tau}\\nu^{r} + \\nu^{r} + D_{s}\\partial_{\\alpha}n\\,\\partial_{r}\\alpha + \\frac{D_{s}}{T}(n + T\\partial_{T}n)\\partial_{r}T",
      "\\tau_{n} = \\frac{D_{s}h}{T},\\qquad \\kappa = D_{s}\\partial_{\\alpha}n"}},
    {"text", "The Fick closure keeps these two coefficients and drops the $\\partial_{r}T$ term. In the code $\\tau_{n}$ is written with a Bessel combination:"},
    {"eq", {"K_{5} = K_{3} + \\frac{8}{z}K_{4},\\quad K_{4} = K_{2} + \\frac{6}{z}K_{3},\\quad K_{1} = K_{3} - \\frac{4}{z}K_{2}\\ \\Rightarrow\\ 2K_{1} - 3K_{3} + K_{5} = \\frac{48}{z^{2}}K_{3}",
      "\\tau_{n} = \\frac{\\sum_{i}g_{i}q_{i}^{2}\\,2\\pi D_{s}T\\,m_{i}^{5}\\,e^{q_{i}\\alpha}(2K_{1} - 3K_{3} + K_{5})/(192\\pi^{3}T^{3})}{\\sum_{i}q_{i}^{2}n_{i}} = \\frac{D_{s}}{T}\\,\\frac{\\sum_{i}q_{i}^{2}n_{i}h_{i}}{\\sum_{i}q_{i}^{2}n_{i}} = \\frac{D_{s}h}{T},\\qquad h_{i} = m_{i}\\frac{K_{3}(z_{i})}{K_{2}(z_{i})}"}},
    {"sub", "The charm closures"},
    {"eq", {"f = g\\,e^{q\\alpha - E/T},\\qquad dP = \\frac{d^{3}p}{(2\\pi)^{3}E},\\qquad I_{nq} = \\frac{1}{(2q + 1)!!}\\int dP\\,E^{n - 2q}p^{2q}f"}},
    {"text", "Angular averages, then the matching in the rest frame:"},
    {"eq", {"p_{x}^{2}\\to\\frac{p^{2}}{3},\\qquad p_{x}^{4}\\to\\frac{p^{4}}{5},\\qquad p_{x}^{2}p_{y}^{2}\\to\\frac{p^{4}}{15}",
      "\\delta f_{\\nu} = f\\,\\frac{\\nu_{a}p^{a}}{P_{0}}:\\quad \\int dP\\,p_{x}\\delta f_{\\nu} = \\frac{I_{21}}{P_{0}}\\nu^{x} = \\nu^{x},\\quad S^{\\tau\\tau x} = \\frac{I_{41}}{P_{0}}\\nu^{x} = \\mathcal{A}\\nu^{x},\\quad S^{xxx} = 3\\frac{I_{42}}{P_{0}}\\nu^{x} = 3\\mathcal{B}\\nu^{x}",
      "\\delta f_{\\pi} = f\\,\\frac{\\pi_{ab}p^{a}p^{b}}{2I_{42}}:\\quad \\int dP\\,p_{x}p_{y}\\delta f_{\\pi} = \\pi^{xy},\\quad S^{\\tau xy} = \\frac{I_{52}}{I_{42}}\\pi^{xy} = \\mathcal{C}\\pi^{xy}",
      "S^{\\tau xx}\\bigl|_{f} = I_{31},\\qquad \\delta M^{ab} = \\pi_{Q}^{ab} + \\Pi_{Q}\\Delta^{ab}\\ \\text{with the same }\\mathcal{C}\\ \\text{(closure)}"}},
    {"text", "One Boltzmann species:"},
    {"eq", {"n = \\frac{g}{2\\pi^{2}}m^{2}T\\,K_{2}(z)e^{q\\alpha},\\ \\ I_{21} = P_{0} = nT,\\ \\ I_{31} = P_{0}h,\\ \\ I_{41} = P_{0}(m^{2} + 5Th),\\ \\ I_{42} = P_{0}Th,\\ \\ I_{52} = I_{42}\\,m\\frac{K_{4}}{K_{3}}"}},
    {"text", "A list of species at one $(T, \\alpha)$, with $w_{i} = q_{i}^{2}n_{i}$, as in \\texttt{charm\\_eos}:"},
    {"eq", {"h = \\frac{\\sum w_{i}h_{i}}{\\sum w_{i}},\\quad \\mathcal{A} = \\frac{\\sum w_{i}m_{i}^{2}}{\\sum w_{i}} + 5Th,\\quad \\mathcal{B} = Th,\\quad \\mathcal{C} = \\frac{\\sum w_{i}h_{i}\\mathcal{C}_{i}}{\\sum w_{i}h_{i}},\\quad I_{31} = nTh",
      "\\partial_{T}h_{i} = -\\frac{m^{2}}{T^{2}}\\frac{-\\tfrac{1}{2}(K_{2} + K_{4})K_{2} + \\tfrac{1}{2}K_{3}(K_{1} + K_{3})}{K_{2}^{2}},\\quad \\partial_{T}\\mathcal{C}_{i} = -\\frac{m^{2}}{T^{2}}\\frac{-\\tfrac{1}{2}(K_{3} + K_{5})K_{3} + \\tfrac{1}{2}K_{4}(K_{2} + K_{4})}{K_{3}^{2}}",
      "\\partial_{T}n_{i} = \\frac{g}{2\\pi^{2}}e^{q\\alpha}\\frac{m^{2}(mK_{1} + 2TK_{2} + mK_{3})}{2T},\\quad \\partial_{\\alpha}I_{31} = \\partial_{\\alpha}n\\,Th,\\quad \\partial_{T}I_{31} = \\partial_{T}n\\,Th + nh + nT\\partial_{T}h"}},
    {"text", "A state of charge $|q| = 2$ enters with these weights; a charge-correct second moment would also carry $\\partial_{\\alpha}\\mathcal{C}$, $\\partial_{\\alpha}h$."},
    {"sub", "The vacuum regulators"},
    {"text", "Two forms. Two-way drivers, $\\epsilon$ = \\texttt{dmn\\_eps}:"},
    {"eq", {"E_{\\alpha}\\to\\nabla_{\\mu}N^{\\mu} + \\epsilon\\,u^{\\mu}\\partial_{\\mu}\\alpha:\\qquad (A_t)_{\\alpha\\alpha}\\to(A_t)_{\\alpha\\alpha} + \\epsilon u^{\\tau},\\quad (A_i)_{\\alpha\\alpha}\\to(A_i)_{\\alpha\\alpha} + \\epsilon u^{i}"}},
    {"text", "Background and 7-field drivers, $\\partial_{\\alpha}n\\to\\partial_{\\alpha}n + \\epsilon$: the same term in the number row. In the consistent rows it also reaches $\\kappa$, since $\\kappa$ is formed inside the kernel:"},
    {"eq", {"\\kappa\\to D_{s}(\\partial_{\\alpha}n + \\epsilon)\\quad\\text{(consistent, second moment)},\\qquad \\kappa\\ \\text{unchanged}\\quad\\text{(Fick)}"}},
    {"text", "It is not conservative:"},
    {"eq", {"\\partial_{\\tau}(\\tau rN^{\\tau}) + \\partial_{r}(\\tau rN^{r}) = -\\tau r\\,\\epsilon\\,u^{\\mu}\\partial_{\\mu}\\alpha"}},
    {"text", "\\texttt{matrix1d\\_visc\\_gamma\\_HQ!} accepts \\texttt{dmn\\_eps}, but its kernel has no $\\partial_{\\alpha}n$: the number row carries $n$, so the regulator does nothing there. $\\alpha$ is clamped to $|\\alpha|\\le\\alpha_{\\max}$ where the EOS is read."},
    {"sub", "The cold edge"},
    {"eq", {"x = \\frac{T - (T_{\\mathrm{edge}} - w_{T})}{w_{T}},\\qquad s(x) = x^{3}(10 - 15x + 6x^{2})",
      "s(0) = s'(0) = s''(0) = s'(1) = s''(1) = 0,\\quad s(1) = 1",
      "D_{s}\\to f(T)D_{s},\\qquad f = f_{0} + (1 - f_{0})s\\ \\text{in the window},\\quad f = 1\\ (T\\ge T_{\\mathrm{edge}}),\\quad f = f_{0}\\ (T\\le T_{\\mathrm{edge}} - w_{T})",
      "w = 1 - s = \\frac{1 - f}{1 - f_{0}}"}},
    {"text", "Only $D_{s}$ changes, and it is not differentiated, as in the matrices. Its gradient is dropped; this is a decoupling prescription, not a derived term."},
    {"sub", "The corona update"},
    {"text", "Below $T_{\\mathrm{edge}}$ the quasi-linear $\\dot\\alpha$ is blended with one that gives a flux-form number update. $V = r$ in 1+1D, $V = 1$ in 2+1D; rates of $T$, $u$, $\\nu$ are those the matrices gave:"},
    {"eq", {"N^{\\tau} = nu^{\\tau} + \\nu^{\\tau},\\qquad N^{i} = nu^{i} + \\nu^{i},\\qquad \\nu^{\\tau} = \\frac{u_{i}\\nu^{i}}{u^{\\tau}}",
      "\\frac{d}{d\\tau}(\\tau VN^{\\tau}) = VN^{\\tau} + \\tau V\\bigl(\\partial_{T}n\\,\\dot T\\,u^{\\tau} + \\partial_{\\alpha}n\\,\\dot\\alpha\\,u^{\\tau} + n\\dot u^{\\tau} + \\dot\\nu^{\\tau}\\bigr) \\equiv \\mathcal{T}",
      "\\dot\\alpha_{\\mathrm{flux}} = \\frac{(\\mathcal{T} - VN^{\\tau})/(\\tau V) - \\partial_{T}n\\,\\dot T\\,u^{\\tau} - n\\dot u^{\\tau} - \\dot\\nu^{\\tau}}{u^{\\tau}(\\partial_{\\alpha}n + \\epsilon)}",
      "\\dot u^{\\tau} = \\frac{u_{i}\\dot u^{i}}{u^{\\tau}},\\qquad \\dot\\nu^{\\tau} = \\frac{\\dot u_{i}\\nu^{i} + u_{i}\\dot\\nu^{i}}{u^{\\tau}} - \\frac{\\nu^{\\tau}\\dot u^{\\tau}}{u^{\\tau}}",
      "\\mathcal{T} = -\\frac{G_{i+1/2} - G_{i-1/2}}{\\Delta r},\\qquad G_{i+1/2} = \\tau r_{i+1/2}F_{i+1/2},\\qquad G_{1/2} = 0\\ \\text{on the axis}\\quad\\text{(1+1D)}",
      "\\mathcal{T} = -\\tau\\Bigl(\\frac{F^{x}_{i+1/2,j} - F^{x}_{i-1/2,j}}{\\Delta x} + \\frac{F^{y}_{i,j+1/2} - F^{y}_{i,j-1/2}}{\\Delta y}\\Bigr)\\quad\\text{(2+1D)}",
      "F_{i+1/2} = \\frac{N^{i}_{L} + N^{i}_{R}}{2} - a\\frac{N^{\\tau}_{R} - N^{\\tau}_{L}}{2},\\qquad a = \\min\\Bigl(1,\\ \\max\\Bigl(\\frac{|N^{i}_{L}|}{N^{\\tau}_{L}},\\ \\frac{|N^{i}_{R}|}{N^{\\tau}_{R}}\\Bigr)\\Bigr)",
      "\\dot\\alpha = (1 - w)\\dot\\alpha_{\\mathrm{ql}} + w\\,\\dot\\alpha_{\\mathrm{flux}}"}},
    {"text", "For $w = 1$ and $\\epsilon = 0$ the sum over cells telescopes to the boundary fluxes:"},
    {"eq", {"\\sum_{i}\\Delta r\\,\\frac{d}{d\\tau}(\\tau r_{i}N^{\\tau}_{i}) = -(G_{\\mathrm{end}} - G_{1/2})"}},
    {"text", "The upwind cap and the blend are numerical choices, not limits of the equations. The update is a discretisation of $E_{\\alpha}$, so the matrices are unchanged."},
    {"sub", "Checks"},
    {"text", "R-: symbolic, in the script. A-: \\texttt{audit\\_matrices.jl} rebuilds each system numerically with ForwardDiff and compares the solved form at random states, with a control that multiplies one coefficient by 1.2. D-: it calls the shipped drivers at physical states and evaluates the corona update on a grid."},
    {"checks", gateLog, {}},
    {"source", {FileNameJoin[{DERIV, "derive_drivers.wls"}]}}};
  rdWrite["drivers", "From the matrices to the solver", items]];
