(* ============================================================================
   derive_two_way_export.wls -- loaded by derive_two_way.wls last.

   Writes one Julia function per system the solver uses, into the generated block of
   that system's matrix file in src/Matrix (only when every check passed and all three
   charts ran):
     S3  one_d_HQ_BG_consistent_matrix      1d_viscous_HQ_cilindrical_BG.jl
         two_d_HQ_BG_consistent_matrix      2d_viscous_HQ_BG.jl
     S4  one_d_HQ_BG_m2_matrix              1d_viscous_HQ_cilindrical_BG_m2.jl
         two_d_HQ_BG_m2_matrix              2d_viscous_HQ_BG_m2.jl
     S5  one_d_two_way_matrix(_conformal)   1d_two_way_cilindrical.jl
         two_d_two_way_matrix(_conformal)   2d_two_way.jl
     S6  one_d_two_way_cM_matrix(_conformal)   1d_two_way_cM_cilindrical.jl
         two_d_two_way_cM_matrix(_conformal)   2d_two_way_cM.jl
   with dPP = 0 (no medium has delta_PiPi) and the 2+1D corotation coupling on (wv = 1).
   and checks
     E-same-<fn>   the simplified exported entries equal the derived ones
     E-free-<fn>   every free symbol of the generated expressions is an argument
     E-jl          the generated Julia, run on 24 random states per function (exact
                   values from the derivation), reproduces them
   ============================================================================ *)

If[doExport,
Print[""]; Print["================ export ================  (", stamp[], ")"];
toPkg = fail == 0 && SubsetQ[charts, {"bj", "cyl", "car"}];   (* a partial-chart run never writes the package *)
If[!toPkg, Print["  a check failed or not all charts ran: the export is checked, but src/Matrix is not written"]];
refDir = FileNameJoin[{outDir, "refpoints"}];
If[DirectoryQ[refDir], DeleteDirectory[refDir, DeleteContents -> True]];   (* this run's states only *)
jlFile = FileNameJoin[{outDir, "generated.jl"}];

bgArgs = <|"0p1d" -> {"T", "dtT"}, "1p1d" -> {"T", "ur", "dtT", "drT", "dtur", "drur"},
  "2p1d" -> {"T", "ux", "uy", "dtT", "dxT", "dyT", "dtux", "dxux", "dyux", "dtuy", "dxuy", "dyuy"}|>;
geoArgs = <|"0p1d" -> {"tau"}, "1p1d" -> {"tau", "rr"}, "2p1d" -> {"tau"}|>;
medGen = {"pres", "dtp", "dtdtp"}; medConf = {"aSB"};
trGen = {"zeta", "visc", "tauS", "tauB", "dpp"}; trConf = {"visc", "tauS", "dpp"};
chFM = {"n", "dndT", "dnda", "h", "dhdT"};
chM2 = {"I31", "dI31dT", "dI31da", "aQ", "bQ", "dbQdT", "cQ", "dcQdT"};
(* {system, chart, conformal} -> {function, file} *)
exportSet = {
  {"S3", "1p1d", False} -> {"one_d_HQ_BG_consistent_matrix", "1d_viscous_HQ_cilindrical_BG.jl"},
  {"S4", "1p1d", False} -> {"one_d_HQ_BG_m2_matrix", "1d_viscous_HQ_cilindrical_BG_m2.jl"},
  {"S5", "1p1d", False} -> {"one_d_two_way_matrix", "1d_two_way_cilindrical.jl"},
  {"S5", "1p1d", True} -> {"one_d_two_way_matrix_conformal", "1d_two_way_cilindrical.jl"},
  {"S3", "2p1d", False} -> {"two_d_HQ_BG_consistent_matrix", "2d_viscous_HQ_BG.jl"},
  {"S4", "2p1d", False} -> {"two_d_HQ_BG_m2_matrix", "2d_viscous_HQ_BG_m2.jl"},
  {"S5", "2p1d", False} -> {"two_d_two_way_matrix", "2d_two_way.jl"},
  {"S5", "2p1d", True} -> {"two_d_two_way_matrix_conformal", "2d_two_way.jl"},
  {"S6", "1p1d", False} -> {"one_d_two_way_cM_matrix", "1d_two_way_cM_cilindrical.jl"},
  {"S6", "1p1d", True} -> {"one_d_two_way_cM_matrix_conformal", "1d_two_way_cM_cilindrical.jl"},
  {"S6", "2p1d", False} -> {"two_d_two_way_cM_matrix", "2d_two_way_cM.jl"},
  {"S6", "2p1d", True} -> {"two_d_two_way_cM_matrix_conformal", "2d_two_way_cM.jl"}};
sysDoc = <|
  "S3" -> "Charm first moment on a prescribed background, consistent closure",
  "S4" -> "Charm first and second moment on a prescribed background",
  "S5" -> "The two-way system, Nabla_mu (T_m + T_Q)^{mu nu} = 0 with the medium and charm rows",
  "S6" -> "The fully coupled system, Nabla_mu (T_m + T_Q + dM)^{mu nu} = 0 with the medium and charm rows and c_M = D_s/T"|>;

argsFor[k_, nm_, conf_] := Module[{a},
  a = Switch[k,
    "S3", Join[bgArgs[nm], chFM, {"Ds"}],
    "S4", Join[bgArgs[nm], chFM, chM2, {"Ds"}],
    "S5" | "S6", Join[If[conf, Join[medConf, trConf], Join[medGen, trGen]], chFM, chM2, {"Ds"}]];
  Join[geoArgs[nm], a]];

utauRule[nm_] := Switch[nm, "0p1d", {}, "1p1d", {Power[1 + ur^2, p_] :> utau^(2 p)},
  "2p1d", {Power[1 + ux^2 + uy^2, p_] :> utau^(2 p)}];
jstr[e_] := StringReplace[ToString[e, InputForm, PageWidth -> Infinity],
  {"Sqrt[" -> "sqrt(", "]" -> ")", "[" -> "("}];
mat2jl[M_, name_] := Module[{n = Length[M], m = Length[First[M]]},
  "    " <> name <> " = SMatrix{" <> ToString[n] <> "," <> ToString[m] <> "}((\n" <>
  StringRiffle[Flatten@Table["        " <> jstr[M[[i, j]]], {j, m}, {i, n}], ",\n"] <> "))\n"];
vec2jl[V_, name_] := "    " <> name <> " = SVector{" <> ToString[Length[V]] <> "}((\n" <>
  StringRiffle[Table["        " <> jstr[V[[i]]], {i, Length[V]}], ",\n"] <> "))\n";

(* the derivation's field names -> the field names of Fluidum's layouts (two_way_1d, two_way_2d, ...) *)
jlFieldNames = <|"0p1d" -> {},
  "1p1d" -> {"al" -> "\[Alpha]", "piphi" -> "piphiphi", "pieta" -> "pietaeta", "PiB" -> "piB",
             "piphiphiQ" -> "piphiphi_Q", "pietaetaQ" -> "pietaeta_Q", "piBQ" -> "piB_Q"},
  "2p1d" -> {"al" -> "\[Alpha]", "pieta" -> "pizz", "PiB" -> "piB",
             "piyyQ" -> "piyy_Q", "pizzQ" -> "pizz_Q", "pixyQ" -> "pixy_Q", "piBQ" -> "piB_Q"}|>;
exported = {}; jlOut = ""; jlByFile = <||>;

Do[Module[{k, nm, conf, res, st, sys, m, vars, fnames, keep, A, src, args, fn, file, free, allowed,
    body, csvRows, hdr, ncs, piIdx},
  {k, nm, conf} = ex[[1]]; {fn, file} = ex[[2]];
  res = results[nm]; st = res["st"]; ncs = Length[st["cs"]];
    sys = res["sys"][k];
    vars = sys["vars"]; fnames = fieldNames /@ vars;
    m = toJ[{sys["A"], sys["src"]}, st] /. {dPP -> 0, wv -> 1};
    {A, src} = m;
    If[conf,
      piIdx = Flatten@Position[fnames, "PiB"];
      keep = Complement[Range[Length[vars]], piIdx];
      A = Table[A[[ic]][[keep, keep]], {ic, ncs}]; src = src[[keep]]; fnames = fnames[[keep]];
      {A, src} = {A, src} /. {pres -> aSB T^4, dtp -> 4 aSB T^3, dtdtp -> 12 aSB T^2, PiB -> 0, zeta -> 0}];
    Print["  ", fn, ": ", Length[fnames], " fields, simplifying ... (", stamp[], ")"];
    Module[{A0 = A, s0 = src, z},
      A = Map[norm, A, {3}]; src = norm /@ src;
      A = A /. utauRule[nm]; src = src /. utauRule[nm];
      (* the exported (simplified, u^tau-renamed) entries must equal the derived ones *)
      z = zeroJ[Flatten[{A, src} - {A0, s0}] /. utau -> Switch[nm, "1p1d", Sqrt[1 + ur^2], "2p1d", Sqrt[1 + ux^2 + uy^2], _, 1], 3, 5];
      gate["E-same-" <> fn, z[[1]],
        "simplified exported entries == derived entries (" <> z[[2]] <> ")"]];
    args = argsFor[k, nm, conf];
    allowed = Join[args, fnames, If[nm === "0p1d", {}, {"utau"}]];
    free = Union@Cases[{A, src}, s_Symbol /; Context[s] === "Global`" :> SymbolName[s], {0, Infinity}];
    gate["E-free-" <> fn, SubsetQ[allowed, free] &&
        FreeQ[{A, src}, Derivative] && FreeQ[{A, src}, Alternatives @@ Join[eosHeads, allHeads[st]]],
      "every free symbol is an argument or a field" <> If[SubsetQ[allowed, free], "", ": EXTRA " <> ToString[Complement[free, allowed]]]];
    (* Julia *)
    body = "\n\"\"\"\n    " <> fn <> "(u, " <> StringRiffle[args, ", "] <> ")\n" <>
      sysDoc[k] <> "; " <> Switch[nm, "0p1d", "Bjorken", "1p1d", "1+1D cylindrical", "2p1d", "2+1D Cartesian"] <>
      If[conf, ", conformal medium", ""] <> ".\nFields u = (" <>
      StringRiffle[fnames, ", "] <> "); returns " <>
      Switch[nm, "0p1d", "(A_t, source)", "1p1d", "(A_t, A_r, source)", "2p1d", "(A_t, A_x, A_y, source)"] <>
      ".\n\"\"\"\nfunction " <> fn <> "(u, " <> StringRiffle[args, ", "] <> ")\n    " <>
      StringRiffle[fnames, ", "] <> " = " <> StringRiffle[Table["u[" <> ToString[i] <> "]", {i, Length[fnames]}], ", "] <> "\n" <>
      Switch[nm, "0p1d", "", "1p1d", "    utau = sqrt(1 + ur^2)\n", "2p1d", "    utau = sqrt(1 + ux^2 + uy^2)\n"] <>
      StringJoin[Table[mat2jl[A[[ic]], {"At", "Ax", "Ay"}[[ic]]] <> "\n", {ic, ncs}]] <>
      vec2jl[src, "source"] <> "\n    return (" <> StringRiffle[Join[{"At", "Ax", "Ay"}[[1 ;; ncs]], {"source"}], ", "] <> ")\nend\n";
    body = juliaNames[body, jlFieldNames[nm]];
    jlOut = jlOut <> body;
    jlByFile[file] = Append[Lookup[jlByFile, file, {}], body];
    (* CSV *)
    SeedRandom[Hash[fn]];
    csvRows = Table[Module[{rep, vals},
        rep = Table[Symbol[s] -> Switch[s,
            "tau", RandomInteger[{5, 100}]/10, "rr", RandomInteger[{3, 90}]/10,
            "T", RandomInteger[{120, 500}]/1000,
            "ur" | "ux" | "uy", RandomInteger[{-15, 15}]/10,
            "aSB", RandomInteger[{30, 60}]/10, "pres", RandomInteger[{1, 50}]/1000,
            "dtp", RandomInteger[{10, 500}]/1000, "dtdtp", RandomInteger[{50, 2000}]/1000,
            "zeta" | "visc", RandomInteger[{1, 50}]/1000, "tauS" | "tauB" | "dpp" | "dPP", RandomInteger[{5, 50}]/10,
            "n" | "dnda" | "I31" | "dI31da", RandomInteger[{1, 500}]/10000, "dndT" | "dI31dT", RandomInteger[{0, 1000}]/1000,
            "h" | "cQ", RandomInteger[{160, 300}]/100, "dhdT" | "dbQdT", RandomInteger[{50, 300}]/100,
            "aQ", RandomInteger[{200, 500}]/100, "bQ", RandomInteger[{30, 150}]/100,
            "dcQdT", RandomInteger[{-200, 200}]/100, "Ds", RandomInteger[{5, 100}]/10,
            "al", RandomInteger[{-20, 20}]/10,
            _, If[StringStartsQ[s, "d"], RandomInteger[{-500, 500}]/1000, RandomInteger[{-100, 100}]/10000]],
          {s, Join[fnames, args]}];
        rep = DeleteDuplicatesBy[rep, First];
        rep = Join[rep, Switch[nm, "0p1d", {}, "1p1d", {utau -> Sqrt[1 + ur^2]}, "2p1d", {utau -> Sqrt[1 + ux^2 + uy^2]}] /. rep];
        vals = N[Flatten[{A, src} /. rep], 17];
        Join[N[(Symbol /@ Join[fnames, args]) /. rep, 17], vals]], {24}];
    hdr = Join[Join[fnames, args],
      Flatten@Table[{"At", "Ax", "Ay"}[[ic]] <> "_" <> ToString[i] <> "_" <> ToString[j], {ic, ncs}, {i, Length[fnames]}, {j, Length[fnames]}],
      Table["src_" <> ToString[i], {i, Length[fnames]}]];
    If[!DirectoryQ[refDir], CreateDirectory[refDir, CreateIntermediateDirectories -> True]];
    Export[FileNameJoin[{refDir, fn <> ".csv"}],
      Prepend[N[csvRows], hdr], "CSV"];
    AppendTo[exported, {fn, Length[fnames], fnames, args, ncs}];
    If[Environment["FLUIDUM_DERIV_RENDER"] === "1",   (* kept for render.wl *)
      renderSys[fn] = <|"A" -> A, "src" -> src, "fields" -> fnames, "chart" -> nm, "sys" -> k, "conf" -> conf,
        "file" -> file|>];
  ], {ex, Select[exportSet, KeyExistsQ[results, #[[1, 2]]] &]}];

Export[jlFile, jlOut, "Text"];
Print["  wrote ", jlFile, "  (", Length[exported], " functions, for the check below)"];

(* E-jl: the generated Julia reproduces the CSV it was exported with *)
Module[{chk, run, outs},
  chk = "
using StaticArrays
include(ARGS[1])
dir = ARGS[2]
worst = 0.0
for f in readdir(dir)
    endswith(f, \".csv\") || continue
    fn = Symbol(replace(f, \".csv\" => \"\"))
    lines = readlines(joinpath(dir, f))
    hdr = strip.(split(lines[1], ','), '\"')
    for ln in lines[2:end]
        v = parse.(Float64, split(ln, ','))
        # inputs: field names first, then the scalar arguments (header order); outputs after
        nf = count(h -> startswith(h, \"src_\"), hdr)
        nin = findfirst(h -> startswith(h, \"At_\"), hdr) - 1
        u = v[1:nf]; a = v[nf+1:nin]
        r = getfield(Main, fn)(u, a...)
        got = vcat([vec(permutedims(Matrix(x))) for x in r[1:end-1]]..., collect(r[end]))
        ref = v[nin+1:end]
        e = maximum(abs.(got .- ref) ./ max.(1.0, abs.(ref)))
        global worst = max(worst, e)
    end
end
println(\"EJL \", worst)
";
  Export[FileNameJoin[{jlDir, "check_generated.jl"}], chk, "Text"];
  run = RunProcess[{"julia", "--startup-file=no", "--project=" <> juliaProject,
     FileNameJoin[{jlDir, "check_generated.jl"}], jlFile, refDir}];
  outs = StringCases[run["StandardOutput"], "EJL " ~~ a : (Except[WhitespaceCharacter] ..) :> ToExpression[StringReplace[a, "e" -> "*^"]]];
  If[outs === {}, Print[run["StandardError"]]];
  gate["E-jl", outs =!= {} && outs[[1]] < 1*^-12,
    "the generated Julia (" <> ToString[Length[exported]] <> " functions) parses, runs, and reproduces its 24-state CSVs: max rel err " <>
      If[outs === {}, "julia failed", sci[outs[[1]]]]];
];
(* into the package: each function into the generated block of its system's file *)
If[toPkg && fail == 0,
  KeyValueMap[spliceJulia[FileNameJoin[{matDir, #1}], StringRiffle[#2, "\n"]] &, jlByFile],
  Print["  NOT writing src/Matrix: a check failed or not all charts ran"]];
];
