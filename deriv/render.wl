(* ============================================================================
   render.wl -- the human-readable form of a derivation: a LaTeX fragment for
   deriv/rendered/derivations.tex and a notebook deriv/rendered/<name>.nb.

   Loaded by the derive_*.wls scripts only when FLUIDUM_DERIV_RENDER=1, after the
   derivation and its checks have run.  It reads the script's own objects (matrices,
   sources, field vectors, the check log) and writes nothing into src/.

   Every entry in a Result section is an expression the script computed, renamed with
   the script's own export rules, simplified for reading, and then compared with the
   unsimplified expression at random exact points.  An entry broken over several lines
   is re-assembled from its pieces and compared again.  Only the Setup equations are
   written here, transcribed from the scripts' header comments.  Problems (an unnamed
   symbol, a mismatch) are printed with the prefix [render] and listed at the end;
   they never change the script's own check counts.
   ============================================================================ *)

rdOutDir = FileNameJoin[{DERIV, "rendered"}];
If[!DirectoryQ[rdOutDir], CreateDirectory[rdOutDir, CreateIntermediateDirectories -> True]];
rdProblems = {};
rdNote[s__] := (AppendTo[rdProblems, StringJoin[ToString /@ {s}]]; Print["  [render] ", s]);
rdStats = <|"entries" -> 0, "verified" -> 0, "broken" -> 0|>;
rdW = 74;   (* line width for a displayed entry, in the units of rdBW (about one letter) *)

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
  "approx" -> "\[TildeTilde]", "in" -> "\[Element]"|>;
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
rdPrep[e_, rules_, ut_, what_ : ""] := Module[{o = e //. rules, d},
  d = TimeConstrained[Simplify[o], 10, o];
  d = d /. ut[[1]];
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
  {il, ir} = rdBreak[inner, w];
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
      StringJoin[("\\texttt{" <> rdEscId[#[[1]]] <> "} & " <> rdEsc[#[[3]]] <> " & " <>
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

(* ============================================================================
   derive_2p1d_viscous.wls
   ============================================================================ *)
renderViscous[] := Module[{ut1, ut2, p1, p2, items, prepNo},
  rdSymTeX = Join[rdBaseTeX, <|"uT" -> "T", "uR" -> "u^{r}", "uX" -> "u^{x}", "uY" -> "u^{y}",
    "pFF" -> "\\pi^{\\phi}{}_{\\phi}", "pEE" -> "\\pi^{\\eta}{}_{\\eta}", "pBB" -> "\\Pi",
    "pYY" -> "\\pi^{yy}", "pXY" -> "\\pi^{xy}", "pres" -> "p", "dtp" -> "p'(T)", "dtdtp" -> "p''(T)",
    "R" -> "r", "visc" -> "\\eta", "zeta" -> "\\zeta", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}",
    "dpp" -> "\\delta_{\\pi\\pi}"|>];
  ut1 = {{Power[1 + uR^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + uR^2 :> utau^2}, {utau -> Sqrt[1 + uR^2]}};
  ut2 = {{Power[1 + uX^2 + uY^2, pw_] :> utau^(2 pw), q_Plus /; q === 1 + uX^2 + uY^2 :> utau^2}, {utau -> Sqrt[1 + uX^2 + uY^2]}};
  p1 = rdPrep[#, subC, ut1, "1+1D"] &;
  p2 = rdPrep[#, subK, ut2, "2+1D"] &;
  items = {
    {"sec", "Viscous hydrodynamics: Israel--Stewart with $\\delta_{\\pi\\pi}$"},
    {"text", "Script: \\texttt{deriv/derive\\_2p1d\\_viscous.wls}"}, {"text", "Writes: \\texttt{src/Matrix/viscous\\_generated.jl}"},
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
    {"subsub", "Co-rotating terms (diagnostics)"},
    {"text", "Not part of the system above; exported alone so that a test can ask whether a 2+1D matrix carries them. The projected vorticity coupling, without a factor $\\tau_{\\pi}$ (\\texttt{two\\_d\\_viscous\\_vorticity}):"},
    {"eq", {"2\\,\\pi_{\\lambda}{}^{\\langle\\mu}\\omega^{\\nu\\rangle\\lambda},\\qquad \\omega^{\\mu\\nu} = \\tfrac{1}{2}\\Delta^{\\mu\\alpha}\\Delta^{\\nu\\beta}(\\nabla_{\\alpha}u_{\\beta} - \\nabla_{\\beta}u_{\\alpha})"}},
    {"text", "and the unprojected commutator (\\texttt{two\\_d\\_viscous\\_notebook\\_commutator}):"},
    {"eq", {"\\tau_{\\pi}\\bigl(\\omega^{\\mu}{}_{\\lambda}\\pi^{\\lambda\\nu} - \\pi^{\\mu\\lambda}\\omega_{\\lambda}{}^{\\nu}\\bigr),\\qquad \\omega_{\\mu\\nu} = \\tfrac{1}{2}(\\nabla_{\\mu}u_{\\nu} - \\nabla_{\\nu}u_{\\mu})"}},
    {"text", "Both carry no source (V-GV-src, V-GN-src). Nonzero entries of their matrices, $A^{(\\omega)}$ for the first and $A^{(\\mathrm{c})}$ for the second, rows and columns as above:"},
    {"sparsity", {{"A^{(\\omega)}_t", rdNZ[Map[p2, AtV, {2}]]}, {"A^{(\\omega)}_x", rdNZ[Map[p2, AxV, {2}]]},
      {"A^{(\\omega)}_y", rdNZ[Map[p2, AyV, {2}]]}}},
    {"sparsity", {{"A^{(\\mathrm{c})}_t", rdNZ[Map[p2, AtN, {2}]]}, {"A^{(\\mathrm{c})}_x", rdNZ[Map[p2, AxN, {2}]]},
      {"A^{(\\mathrm{c})}_y", rdNZ[Map[p2, AyN, {2}]]}}},
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

renderHQ[] := Module[{jr, pr, items},
  rdSymTeX = rdHQTeX[];
  jr = jrepl /. {Verbatim[dn_da] -> rdDnda, Verbatim[dn_dT] -> rdDndT};
  pr = rdPrep[#, jr, rdUtCar, "2+1D"] &;
  items = {
    {"sec", "Charm first moment on a 2+1D background ($\\nabla\\alpha$ closure)"},
    {"text", "Script: \\texttt{deriv/derive\\_hq\\_2p1d.wls}"}, {"text", "Writes: \\texttt{src/Matrix/HQ\\_2p1d\\_BG\\_generated.jl}"},
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

renderHQC[] := Module[{jr, pr, items},
  rdSymTeX = rdHQTeX[];
  jr = jrepl /. {Verbatim[dn_da] -> rdDnda, Verbatim[dn_dT] -> rdDndT};
  pr = rdPrep[#, jr, rdUtCar, "2+1D"] &;
  items = {
    {"sec", "Consistent charm first moment on a 2+1D background"},
    {"text", "Script: \\texttt{deriv/derive\\_hq\\_consistent\\_2p1d.wls}"}, {"text", "Writes: \\texttt{src/Matrix/HQ\\_2p1d\\_BG\\_consistent\\_generated.jl}"},
    {"sub", "Setup"},
    Splice[rdHQSetupCommon[]],
    {"text", "The first moment is the projected divergence of the charm stress tensor:"},
    {"eq", {"\\nabla_{\\mu}(n\\,u^{\\mu} + \\nu^{\\mu}) = 0,\\qquad \\frac{D_{s}}{T}\\,\\Delta^{\\mu}{}_{\\lambda}\\nabla_{\\nu}T_{Q}^{\\nu\\lambda} + \\nu^{\\mu} = 0",
      "T_{Q}^{\\alpha\\beta} = (n h - P_{0})\\,u^{\\alpha}u^{\\beta} + P_{0}\\,\\Delta^{\\alpha\\beta} + h\\,(u^{\\alpha}\\nu^{\\beta} + \\nu^{\\alpha}u^{\\beta}),\\qquad P_{0} = n T",
      "n = n(T, \\alpha),\\quad h = h(T),\\qquad \\tau_{n} = D_{s}h/T,\\quad \\kappa = D_{s}\\,\\partial_{\\alpha}n"}},
    {"text", "At rest the current relaxes to"},
    {"eq", {"\\nu^{i} = -\\kappa\\,\\partial_{i}\\alpha - \\frac{D_{s}}{T}\\,(n + T\\,\\partial_{T}n)\\,\\partial_{i}T"}},
    {"text", "Relative to the $\\nabla\\alpha$ row of the previous section, the current equation gains (HC-C3, HC-C5)"},
    {"eq", {"\\frac{D_{s}}{T}(n + T\\,\\partial_{T}n)\\,\\nabla^{\\mu}_{\\perp}T + \\tau_{n} n\\,a^{\\mu} + \\tau_{n}\\Delta^{\\mu}{}_{\\lambda}\\,\\nu^{\\sigma}\\nabla_{\\sigma}u^{\\lambda} + \\nu^{\\mu}\\Bigl(\\tau_{n}\\theta + \\frac{D_{s}}{T}\\,Dh\\Bigr),\\qquad a^{\\mu} = Du^{\\mu}"}},
    {"text", "Fields and rows as before:"},
    {"eq", {"\\phi = " <> rdTuple[rdTeX /@ (varsK //. jr)] <> ",\\qquad \\nu^{\\tau} = " <> rdTeX[pr[nucar[[1]]]],
      "A_t\\,\\partial_{\\tau}\\phi + A_x\\,\\partial_{x}\\phi + A_y\\,\\partial_{y}\\phi + S = 0"}},
    {"sub", "Result"},
    {"entries", rdMat[AtK, "A_t", pr]},
    {"entries", rdMat[AxK, "A_x", pr]},
    {"entries", rdMat[AyK, "A_y", pr]},
    {"entries", rdVec[srcK, "S", pr]},
    {"sub", "Checks"},
    {"checks", gateLog, {}},
    {"source", {FileNameJoin[{DERIV, "derive_hq_consistent_2p1d.wls"}]}}};
  rdWrite["hq_consistent_2p1d", "Consistent charm first moment on a 2+1D background", items]];

(* ============================================================================
   derive_two_way.wls (with _gates and _export)
   ============================================================================ *)
rdTwoWayTeX[] := Module[{fld, m},
  fld = <|"T" -> "T", "ur" -> "u^{r}", "ux" -> "u^{x}", "uy" -> "u^{y}",
    "piphi" -> "\\pi^{\\phi}{}_{\\phi}", "pieta" -> "\\pi^{\\eta}{}_{\\eta}", "piyy" -> "\\pi^{yy}",
    "pixy" -> "\\pi^{xy}", "PiB" -> "\\Pi", "al" -> "\\alpha", "nur" -> "\\nu^{r}", "nux" -> "\\nu^{x}",
    "nuy" -> "\\nu^{y}", "pl" -> "\\pi_{Q}^{\\ell\\ell}", "pphi" -> "\\pi_{Q}^{\\hat\\phi\\hat\\phi}",
    "PiQ" -> "\\Pi_{Q}", "qeta" -> "\\pi_{Q}{}^{\\eta}{}_{\\eta}", "pQxx" -> "\\pi_{Q}^{xx}",
    "pQxy" -> "\\pi_{Q}^{xy}", "pQyy" -> "\\pi_{Q}^{yy}"|>;
  m = Join[rdBaseTeX, fld, <|"pres" -> "P", "dtp" -> "P'(T)", "dtdtp" -> "P''(T)", "aSB" -> "a_{\\mathrm{SB}}",
    "zeta" -> "\\zeta", "visc" -> "\\eta", "tauS" -> "\\tau_{\\pi}", "tauB" -> "\\tau_{\\Pi}",
    "dpp" -> "\\delta_{\\pi\\pi}", "dPP" -> "\\delta_{\\Pi\\Pi}", "n" -> "n", "dndT" -> "(\\partial_{T} n)",
    "dnda" -> "(\\partial_{\\alpha} n)", "h" -> "h", "dhdT" -> "h'(T)", "I31" -> "I_{31}",
    "dI31dT" -> "(\\partial_{T} I_{31})", "dI31da" -> "(\\partial_{\\alpha} I_{31})", "aQ" -> "\\mathcal{A}",
    "bQ" -> "\\mathcal{B}", "dbQdT" -> "\\mathcal{B}'(T)", "cQ" -> "\\mathcal{C}", "dcQdT" -> "\\mathcal{C}'(T)",
    "Ds" -> "D_{s}", "c" -> "c", "wv" -> "w_{\\omega}"|>];
  (* background derivatives named by the export: dt/dr/dx/dy + field name *)
  Do[m[d[[1]] <> f] = rdDeriv[d[[2]], fld[f]],
    {d, {{"dt", "\\tau"}, {"dr", "r"}, {"dx", "x"}, {"dy", "y"}}}, {f, {"T", "ur", "ux", "uy"}}];
  m];

renderTwoWay[] := Module[{items, full, sysLab, chartLab, dirLab, fnItems, skipped, ids, stJ, fl, conv},
  rdSymTeX = rdTwoWayTeX[];
  stJ[st_] := (toJ[#, st] &);
  fl[st_, key_] := rdTuple[rdTeX /@ (fieldNames /@ st[key] /. s_String :> Symbol[s])];
  sysLab = <|"S1" -> "S1, ideal hydro", "S2" -> "S2, viscous hydro", "S3" -> "S3, charm first moment on a background",
    "S4" -> "S4, charm first and second moment on a background", "S5" -> "S5, two-way"|>;
  chartLab = <|"0p1d" -> "0+1D Bjorken", "1p1d" -> "1+1D cylindrical", "2p1d" -> "2+1D Cartesian"|>;
  dirLab = <|"0p1d" -> {"A_t"}, "1p1d" -> {"A_t", "A_r"}, "2p1d" -> {"A_t", "A_x", "A_y"}|>;
  (* typeset in full: every 0+1D system and the 1+1D conformal ones *)
  full[r_] := r["chart"] === "0p1d" || (r["chart"] === "1p1d" && r["conf"]);
  fnItems[fn_, r_] := Module[{nm = r["chart"], labs = dirLab[r["chart"]], phi, head, ents},
    phi = rdTuple[rdTeX /@ (Symbol /@ r["fields"])];
    head = {{"subsub", sysLab[r["sys"]] <> If[r["conf"], ", conformal", ""]},
      {"text", "\\texttt{" <> rdEscId[fn] <> "}, fields $\\phi = " <> phi <> "$."}};
    If[full[r],
      ents = Join[Flatten[Table[rdMat[r["A"][[ic]], labs[[ic]], Identity], {ic, Length[labs]}], 1],
        rdVec[r["src"], "S", Identity]];
      Join[head, {{"entries", ents}}],
      Join[head, {{"sparsity", Append[Table[{labs[[ic]], rdNZ[r["A"][[ic]]]}, {ic, Length[labs]}],
        {"S", List /@ (# =!= 0 & /@ r["src"])}]}}]]];
  ids = gateLog[[All, 1]];
  skipped = Select[{"R-S3-1p1d-generator", "R-S4-1p1d-generator", "R-S4c-1p1d-generator"}, !MemberQ[ids, #] &];
  skipped = {#, "These compare with reference CSVs of the 1+1D derivations of the consistent first and second moment and of the $c_M$ coupling, which are read only when \\texttt{FLUIDUM\\_EXTERNAL\\_REFPOINTS} names a directory holding them; it was " <>
      If[csvDir === "", "unset", "set, but the CSVs were not found"] <> "."} & /@ skipped;
  items = {
    {"sec", "Medium, charm and the two-way coupled system"},
    {"text", "Script: \\texttt{deriv/derive\\_two\\_way.wls}, with \\texttt{derive\\_two\\_way\\_gates.wls} and \\texttt{derive\\_two\\_way\\_export.wls}"}, {"text", "Writes: \\texttt{src/Matrix/two\\_way\\_generated.jl}"},
    {"sub", "Setup"},
    {"text", "Three charts, boost invariant: 0+1D Bjorken (fields of $\\tau$ only), 1+1D cylindrical, 2+1D Cartesian Milne:"},
    {"eq", {rdDiag[BJ] <> ",\\quad u^{\\mu} = " <> rdTuple[rdTeX /@ stateBJ["u"]], rdDiag[CYL], rdDiag[CAR]}},
    {"text", "Medium, with the EOS $P(T)$ kept symbolic:"},
    {"eq", {"T_{m}^{\\mu\\nu} = (e + P + \\Pi)\\,u^{\\mu}u^{\\nu} + (P + \\Pi)\\,g^{\\mu\\nu} + \\pi^{\\mu\\nu},\\qquad e = T P'(T) - P(T)",
      "\\tau_{\\pi}\\,\\Delta^{\\mu\\nu}_{\\alpha\\beta}D\\pi^{\\alpha\\beta} + \\pi^{\\mu\\nu} + 2\\eta\\,\\sigma^{\\mu\\nu} + \\delta_{\\pi\\pi}\\,\\pi^{\\mu\\nu}\\theta = 0,\\qquad \\tau_{\\Pi}D\\Pi + \\Pi + \\zeta\\theta + \\delta_{\\Pi\\Pi}\\Pi\\theta = 0"}},
    {"text", "($D$, $\\theta$, $\\Delta^{\\mu\\nu}$, $\\Delta^{\\mu\\nu}_{\\alpha\\beta}$, $\\sigma^{\\mu\\nu}$ as in the viscous section.) Charm, with $n = n(T, \\alpha)$, $h = h(T)$, $P_{0} = nT$ and $\\delta M^{\\mu\\nu} = \\pi_{Q}^{\\mu\\nu} + \\Pi_{Q}\\Delta^{\\mu\\nu}$:"},
    {"eq", {"\\nabla_{\\mu}(n\\,u^{\\mu} + \\nu^{\\mu}) = 0,\\qquad \\frac{D_{s}}{T}\\,\\Delta^{\\mu}{}_{\\lambda}\\nabla_{\\nu}T_{Q}^{\\nu\\lambda} + \\nu^{\\mu} = 0",
      "T_{Q}^{\\mu\\nu} = (n h - P_{0})\\,u^{\\mu}u^{\\nu} + P_{0}\\,\\Delta^{\\mu\\nu} + h\\,(u^{\\mu}\\nu^{\\nu} + \\nu^{\\mu}u^{\\nu}) + c\\,\\delta M^{\\mu\\nu}",
      "S^{\\lambda\\alpha\\beta} = I_{31}\\,(u^{\\lambda}\\Delta^{\\alpha\\beta} + 2\\ \\text{perm.}) + \\mathcal{A}\\,(u^{\\lambda}u^{\\alpha}\\nu^{\\beta} + 2\\ \\text{perm.}) + \\mathcal{B}\\,(\\Delta^{\\lambda\\alpha}\\nu^{\\beta} + 2\\ \\text{perm.}) + \\mathcal{C}\\,(u^{\\lambda}\\delta M^{\\alpha\\beta} + 2\\ \\text{perm.})",
      "\\mathcal{E}^{\\mu\\nu} = \\Delta^{\\mu}{}_{\\alpha}\\Delta^{\\nu}{}_{\\beta}\\nabla_{\\lambda}S^{\\lambda\\alpha\\beta}:\\qquad \\frac{D_{s}}{2T}\\Bigl(\\mathcal{E}^{\\mu\\nu} - \\tfrac{1}{3}\\Delta^{\\mu\\nu}\\mathcal{E}^{\\lambda}{}_{\\lambda}\\Bigr) + \\pi_{Q}^{\\mu\\nu} = 0,\\qquad \\frac{D_{s}}{2T}\\,\\tfrac{1}{3}\\mathcal{E}^{\\lambda}{}_{\\lambda} + \\Pi_{Q} = 0"}},
    {"text", "``2 perm.'' are the two other placements of the lone index. $I_{31}(T, \\alpha)$, $\\mathcal{A}(T)$, $\\mathcal{B}(T)$, $\\mathcal{C}(T)$ stay symbolic (Boltzmann: $I_{31} = nTh$, $\\mathcal{A} = m^{2} + 5Th$, $\\mathcal{B} = Th$, $\\mathcal{C} = m K_{4}/K_{3}$). $c = 1$ is $c_{M} = D_{s}/T$; $c = 0$ drops $\\delta M$ from $T_{Q}$. In 2+1D the traceless rows also carry"},
    {"eq", {"-(1 - w_{\\omega})\\,\\tau_{M}W^{\\mu\\nu},\\qquad \\tau_{M} = \\frac{D_{s}\\mathcal{C}}{2T},\\qquad W^{\\mu\\nu} = \\pi_{Q}^{\\mu}{}_{\\lambda}\\omega^{\\lambda\\nu} + \\pi_{Q}^{\\nu}{}_{\\lambda}\\omega^{\\lambda\\mu},\\qquad \\omega^{\\mu\\nu} = \\tfrac{1}{2}(\\Delta^{\\mu\\lambda}\\nabla_{\\lambda}u^{\\nu} - \\Delta^{\\nu\\lambda}\\nabla_{\\lambda}u^{\\mu})"}},
    {"text", "so $w_{\\omega} = 1$ is the covariant row and $w_{\\omega} = 0$ removes the corotation coupling. The five systems:"},
    {"eq", {"\\text{S1}:\\ \\nabla_{\\mu}T_{m}^{\\mu\\nu} = 0,\\ \\pi = \\Pi = 0\\qquad \\text{S2}:\\ \\text{S1} + \\text{the }\\pi, \\Pi\\text{ equations}",
      "\\text{S3}:\\ \\text{charm number} + \\text{first moment, background }T, u\\qquad \\text{S4}:\\ \\text{S3} + \\text{second moment, background }T, u",
      "\\text{S5}:\\ \\nabla_{\\mu}(T_{m} + T_{Q})^{\\mu\\nu} = 0,\\ \\text{the }\\pi, \\Pi\\text{ equations, the charm rows of S4; all fields evolved}"}},
    {"text", "Conformal: $P = a_{\\mathrm{SB}}T^{4}$, no bulk. Source on the left, $A_t\\,\\partial_{\\tau}\\phi + A_r\\,\\partial_{r}\\phi + S = 0$ in 1+1D and $A_t\\,\\partial_{\\tau}\\phi + A_x\\,\\partial_{x}\\phi + A_y\\,\\partial_{y}\\phi + S = 0$ in 2+1D."},
    {"subsub", "Fields"},
    {"text", "Medium, charm first moment, charm second moment, chart by chart:"},
    {"eq", {"\\text{0+1D}:\\ " <> fl[stateBJ, "med"] <> ",\\ " <> fl[stateBJ, "ch"] <> ",\\ " <> fl[stateBJ, "m2"],
      "\\text{1+1D}:\\ " <> fl[stateCYL, "med"] <> ",\\ " <> fl[stateCYL, "ch"] <> ",\\ " <> fl[stateCYL, "m2"],
      "\\text{2+1D}:\\ " <> fl[stateCAR, "med"] <> ",\\ " <> fl[stateCAR, "ch"] <> ",\\ " <> fl[stateCAR, "m2"]}},
    {"text", "Rows, in the same order: the fluid rows $\\nabla_{\\mu}T^{\\mu\\tau}$, $\\nabla_{\\mu}T^{\\mu i}$ ($T_{m} + T_{Q}$ in S5); the shear rows $E^{\\eta\\eta}$ (0+1D), $E^{\\phi\\phi}$, $E^{\\eta\\eta}$ (1+1D), $E^{yy}$, $E^{\\eta\\eta}$, $E^{xy}$ (2+1D); the bulk row; $\\nabla_{\\mu}N^{\\mu}$; the current components; the traceless second-moment rows, components $\\hat\\eta\\hat\\eta$ (0+1D), $\\ell\\ell$, $\\hat\\phi\\hat\\phi$ (1+1D), $xx$, $xy$, $yy$ (2+1D); the trace row."},
    {"text", "In 0+1D the tensors are diagonal and $\\nu = 0$:"},
    {"eq", {"\\pi^{\\mu\\nu} = \\mathrm{diag}" <> rdTuple[rdTeX /@ Diagonal[toJ[stateBJ["piM"], stateBJ]]] <>
      ",\\qquad \\pi_{Q}^{\\mu\\nu} = \\mathrm{diag}" <> rdTuple[rdTeX /@ Diagonal[toJ[stateBJ["piQ"], stateBJ]]]}},
    {"text", "In 1+1D the medium $\\pi^{\\mu\\nu}$ is as in the viscous section, and the charm stress is spanned by an orthonormal triad transverse to $u$:"},
    {"eq", {"\\pi_{Q}^{\\mu\\nu} = \\pi_{Q}^{\\ell\\ell}\\,\\ell^{\\mu}\\ell^{\\nu} + \\pi_{Q}^{\\hat\\phi\\hat\\phi}\\,\\hat\\phi^{\\mu}\\hat\\phi^{\\nu} - (\\pi_{Q}^{\\ell\\ell} + \\pi_{Q}^{\\hat\\phi\\hat\\phi})\\,\\hat\\eta^{\\mu}\\hat\\eta^{\\nu}",
      "\\ell^{\\mu} = " <> rdTuple[rdTeX /@ (toJ[stateCYL["triad"][[1]], stateCYL] /. Sqrt[1 + ur^2] -> utau)] <>
      ",\\quad \\hat\\phi^{\\mu} = " <> rdTuple[rdTeX /@ stateCYL["triad"][[2]]] <> ",\\quad \\hat\\eta^{\\mu} = " <>
      rdTuple[rdTeX /@ stateCYL["triad"][[3]]]}},
    {"text", "In 2+1D the medium $\\pi^{\\mu\\nu}$ is as in the viscous section; $\\pi_{Q}^{xx}$, $\\pi_{Q}^{xy}$, $\\pi_{Q}^{yy}$ are independent, the $\\tau$ row follows from $u_{\\mu}\\pi_{Q}^{\\mu\\nu} = 0$ and"},
    {"eq", {"\\pi_{Q}{}^{\\eta}{}_{\\eta} = " <> rdTeX[TimeConstrained[Simplify[toJ[stateCAR["qet"], stateCAR]], 20,
        toJ[stateCAR["qet"], stateCAR]] /. rdUtTW2[[1]]]}},
    {"sub", "Result"},
    {"text", "Every 0+1D system and the 1+1D conformal systems in full; the generic 1+1D systems and every 2+1D system as the pattern of nonzero entries ($\\bullet$), rows and columns in the field order given. Their entries are in \\texttt{src/Matrix/two\\_way\\_generated.jl}, one function per system."}};
  Do[Module[{fns = Select[exported[[All, 1]], renderSys[#]["chart"] === nm &]},
      AppendTo[items, {"sub", chartLab[nm]}];
      Do[items = Join[items, fnItems[fn, renderSys[fn]]], {fn, fns}]],
    {nm, Select[{"0p1d", "1p1d", "2p1d"}, KeyExistsQ[results, #] &]}];
  items = Join[items, {
    {"sub", "Checks"},
    {"checks", gateLog, skipped},
    {"text", "Not a check: the characteristic speeds of S5 are computed and reported (\\texttt{two\\_way\\_char\\_speeds.tsv} in the scratch directory)."},
    {"source", {FileNameJoin[{DERIV, "derive_two_way.wls"}], FileNameJoin[{DERIV, "derive_two_way_gates.wls"}],
      FileNameJoin[{DERIV, "derive_two_way_export.wls"}]}}}];
  rdWrite["two_way", "Medium, charm and the two-way coupled system", items]];
rdUtTW2 = {{Power[1 + ux^2 + uy^2, pw_] :> utau^(2 pw)}, {utau -> Sqrt[1 + ux^2 + uy^2]}};
