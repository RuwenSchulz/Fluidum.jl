(* ============================================================================
   julia_source_reader.wl -- read a Fluidum matrix out of its Julia source.

   A derivation is compared with the .jl file itself, not with a hand transcription
   of it, so that an edit to the code is seen by the comparison. Used by the
   1+1D reduction checks of the 2+1D and two-way generators.

     juliaMatrices[file, fname]  ->  <|"At" -> {{..}}, "Ax" -> {{..}}, "source" -> {..}, ...|>

   for every `NAME = SMatrix{a,b}(...)` / `NAME = SVector{n}(...)` inside
   `function fname(`.  Entries are WL expressions in the function's own argument names
   (`u[2]` stays `u[2]`, in the context jl`: write jl["u"][2]), with Julia `_` mapped to `U` (dP_dT -> dPUdT), `κ` to
   `kappaG`, broadcasting dots dropped, and floats rationalised (1.5 -> 3/2, `3.` -> 3).
   SMatrix entries are column-major, as in Julia.

   Supported Julia: + - * / ^ (infix and `^(a,b)`), sqrt / exp / log (also `sqrt.(`),
   indexing, literals.  Anything else fails at ToExpression with an error.
   ============================================================================ *)

juliaFunctionText[file_String, fname_String] := Module[{txt, p0, rest, pend},
  txt = Import[file, "Text"];
  p0 = StringPosition[txt, RegularExpression["(?m)^(@\\w+\\s+)*function\\s+" <> fname <> "\\("]];
  If[Length[p0] != 1, Throw[{"juliaFunctionText: not exactly one definition", file, fname, Length[p0]}]];
  rest = StringDrop[txt, p0[[1, 1]] - 1];
  pend = StringPosition[rest, RegularExpression["(?m)^end\\b"]];
  If[pend === {}, Throw[{"juliaFunctionText: no closing `end`", file, fname}]];
  (* strip comments: everything from # to end of line *)
  StringReplace[StringTake[rest, pend[[1, 1]] - 1], RegularExpression["#[^\n]*"] -> ""]];

(* the text between the "(" at position p (1-based, pointing AT the paren) and its match *)
balancedFrom[s_String, p_Integer] := Module[{depth = 0, i = p, c},
  While[True,
    c = StringTake[s, {i}];
    Which[c === "(", depth++, c === ")", depth--];
    If[depth == 0, Break[]]; i++];
  StringTake[s, {p + 1, i - 1}]];

(* split on commas at bracket depth 0 *)
topSplit[s_String] := Module[{depth = 0, cur = {}, out = {}, c},
  Do[c = StringTake[s, {i}];
    Which[MemberQ[{"(", "["}, c], depth++, MemberQ[{")", "]"}, c], depth--];
    If[c === "," && depth == 0, AppendTo[out, StringJoin[cur]]; cur = {}, AppendTo[cur, c]],
    {i, StringLength[s]}];
  AppendTo[out, StringJoin[cur]];
  Select[StringTrim /@ out, # =!= "" &]];

julia2wl[s0_String] := Module[{s = s0, out = {}, stack = {}, prev = "", c},
  s = StringReplace[s, {"\[Kappa]" -> "kappaG", "κ" -> "kappaG"}];
  s = StringReplace[s, {"sqrt.(" -> "sqrt(", "exp.(" -> "exp(", "log.(" -> "log("}];
  s = StringReplace[s, {".^" -> "^", ".*" -> "*", "./" -> "/", ".+" -> "+", ".-" -> "-"}];
  s = StringReplace[s, RegularExpression["(\\d)\\.(?!\\d)"] -> "$1"];      (* 3. -> 3 *)
  s = StringReplace[s, {"^(" -> "Power(", "sqrt(" -> "Sqrt(", "exp(" -> "Exp(", "log(" -> "Log("}];
  s = StringReplace[s, "_" -> "U"];
  Do[c = StringTake[s, {i}];
    Which[
      c === "(", If[StringMatchQ[prev, RegularExpression["[A-Za-z0-9]"]],
                    AppendTo[stack, True]; AppendTo[out, "["],
                    AppendTo[stack, False]; AppendTo[out, "("]],
      c === ")", If[Last[stack], AppendTo[out, "]"], AppendTo[out, ")"]]; stack = Most[stack],
      True, AppendTo[out, c]];
    If[!StringMatchQ[c, WhitespaceCharacter], prev = c],
    {i, StringLength[s]}];
  (* parse into the private context jl` so a Julia argument called `r` or `X` cannot collide
     with a global of the calling script. Refer to them as jl["r"]. *)
  Rationalize[Block[{$Context = "jl`", $ContextPath = {"System`"}}, ToExpression[StringJoin[out]]], 0]];

jl[name_String] := Symbol["jl`" <> name];

juliaMatrices[file_String, fname_String] := Module[{body, hits, res = <||>, name, dims, arg, ents},
  body = juliaFunctionText[file, fname];
  hits = StringPosition[body, RegularExpression["\\b(\\w+)\\s*=\\s*S(Matrix|Vector)\\{([0-9, ]+)\\}\\("], Overlaps -> False];
  Do[
    With[{m = StringTake[body, h]},
      name = First @ StringCases[m, RegularExpression["^(\\w+)"] -> "$1"];
      dims = ToExpression["{" <> First @ StringCases[m, RegularExpression["\\{([0-9, ]+)\\}"] -> "$1"] <> "}"];
      arg = StringTrim @ balancedFrom[body, h[[2]]];
      (* SMatrix{n,n}((a, b, ...)) : strip the tuple parens *)
      If[StringStartsQ[arg, "("] && StringLength[balancedFrom[arg, 1]] == StringLength[arg] - 2,
        arg = StringTake[arg, {2, -2}]];
      ents = julia2wl /@ topSplit[arg];
      If[Length[ents] != Times @@ dims, Throw[{"juliaMatrices: entry count", fname, name, Length[ents], dims}]];
      res[name] = If[Length[dims] == 1, ents, Transpose[Partition[ents, dims[[1]]]]]],
    {h, hits}];
  res];
