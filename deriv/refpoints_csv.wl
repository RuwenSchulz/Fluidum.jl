(* ============================================================================
   refpoints_csv.wl -- write reference states (the derivation evaluated at exact
   points) as CSV, so a Julia test can check the code against the derivation
   without Mathematica.

   Long format, one value per row:  system, point, key, value
     key = "in:<argument name>"        an input of the Julia function
           "<matrix>:<i>:<j>", "<vector>:<i>"   an output entry (1-based)
   Values printed to 18 significant digits from a 30-digit evaluation.
   ============================================================================ *)

fmtNum[x_] := Module[{v = N[x, 30]},
  If[!NumericQ[v], Throw[{"fmtNum: not numeric", x}]];
  ToString[NumberForm[v, 18, NumberFormat -> (If[#3 == "", #1, #1 <> "e" <> #3] &),
    ExponentFunction -> (If[-3 <= # <= 3, Null, #] &)]]];
matRows[stem_String, M_] := Flatten[Table[{stem <> ":" <> ToString[i] <> ":" <> ToString[j], M[[i, j]]},
  {i, Length[M]}, {j, Length[M[[1]]]}], 1];
vecRows[stem_String, V_] := Table[{stem <> ":" <> ToString[i], V[[i]]}, {i, Length[V]}];
(* inputs: {{julia name, WL expression standing for it, exact value}, ...} -- the WL side may be a
   pattern-looking symbol such as `dn_da` (the generators' export trick), hence Verbatim;
   an Association symbol -> value is accepted too.  outs: {{key, expression}, ...}. *)
refRows[name_String, p_Integer, inputs_List, outs_List] := With[
  {rules = (Verbatim[#[[2]]] -> #[[3]]) & /@ inputs},
  Join[({name, p, "in:" <> #[[1]], fmtNum[#[[3]]]} & /@ inputs),
       ({name, p, #[[1]], fmtNum[#[[2]] /. rules]} & /@ outs)]];
refRows[name_String, p_Integer, inputs_Association, outs_List] :=
  refRows[name, p, KeyValueMap[{SymbolName[#1], #1, #2} &, inputs], outs];
writeRefpoints[path_String, rows_List] := (
  Export[path, Prepend[rows, {"system", "point", "key", "value"}], "CSV"];
  Print["  wrote ", path, "  (", Length[rows], " rows)"]);
