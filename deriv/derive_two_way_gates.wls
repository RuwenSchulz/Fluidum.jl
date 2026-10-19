(* ============================================================================
   derive_two_way_gates.wls -- loaded by derive_two_way.wls after the per-chart
   build (it uses `results`, `gate`, `zeroTest`, the states and the rows).

   G-reduce   the generic-EOS systems against the other generators and the
              hand-written matrices
   G-chart    2+1D -> 1+1D (axisymmetric, on the x axis) and 1+1D -> 0+1D
   G-bjorken  closed forms
   G-twoway   no momentum exchange at rest; Newton's third law for the drag
   speeds     characteristic speeds of S5 at c = 0 and c = 1 (reported, not checked)

   Three G-reduce checks (R-S3/R-S4/R-S4c) compare with reference CSVs of the 1+1D
   derivations of the consistent first moment, the consistent second moment and the
   c_M coupling.  They run only when FLUIDUM_EXTERNAL_REFPOINTS names a directory
   holding those CSVs, and are skipped otherwise.
   ============================================================================ *)

(* ------------------------------------------------ export names (also used by
   derive_two_way_export.wls) ------------------------------------------------- *)
fieldNames = <|
  TTb -> "T", peeb -> "pieta", PiBb -> "PiB", alb -> "al", qeb -> "qeta", PiQb -> "PiQ",
  TTc -> "T", urc -> "ur", pffc -> "piphi", peec -> "pieta", PiBc -> "PiB", alc -> "al",
  nurc -> "nur", plc -> "pl", pphc -> "pphi", PiQc -> "PiQ",
  TTk -> "T", uxk -> "ux", uyk -> "uy", pyyk -> "piyy", peek -> "pieta", pxyk -> "pixy",
  PiBk -> "PiB", alk -> "al", nuxk -> "nux", nuyk -> "nuy", qxxk -> "pQxx", qxyk -> "pQxy",
  qyyk -> "pQyy", PiQk -> "PiQ"|>;
dirName[nm_, k_] := Switch[nm, "0p1d", {"dt"}, "1p1d", {"dt", "dr"}, "2p1d", {"dt", "dx", "dy"}][[k]];
jRules[st_] := Module[{nm = st["name"], cs = st["cs"], nar = Length[st["cs"]], fr = {}, dr = {}},
  Do[With[{nmf = fieldNames[h]},
      AppendTo[fr, (h @@ cs) -> Symbol[nmf]];
      Do[With[{ord = UnitVector[nar, k]},
          AppendTo[dr, (Derivative @@ ord)[h] @@ cs -> Symbol[dirName[nm, k] <> nmf]]], {k, nar}]],
    {h, allHeads[st]}];
  Join[dr, fr, {
    Derivative[2][pp][_] -> dtdtp, Derivative[1][pp][_] -> dtp, pp[_] -> pres,
    Derivative[1, 0][nQ][__] -> dndT, Derivative[0, 1][nQ][__] -> dnda, nQ[__] -> n,
    Derivative[1][hQ][_] -> dhdT, hQ[_] -> h,
    Derivative[1, 0][i31Q][__] -> dI31dT, Derivative[0, 1][i31Q][__] -> dI31da, i31Q[__] -> I31,
    aQf[_] -> aQ, Derivative[1][bQf][_] -> dbQdT, bQf[_] -> bQ,
    Derivative[1][cQf][_] -> dcQdT, cQf[_] -> cQ, csw -> c}]];
toJ[e_, st_] := e /. jRules[st];

(* plain-symbol exact zero test (velocities from Pythagorean tuples) *)
zeroJ[e_, k_: 4, t_: 15] := Module[{s, syms, ok = True, rule},
  s = TimeConstrained[Simplify[e], t, $Aborted];
  If[s =!= $Aborted, Return[{AllTrue[Flatten[{s}], # === 0 &], "symbolic"}]];
  syms = Union@Cases[e, _Symbol?(! NumericQ[#] && Context[#] === "Global`" &), {0, Infinity}];
  Do[SeedRandom[7000 + i];
    rule = Table[sy -> rpos[], {sy, syms}];
    With[{p1 = RandomChoice[pyth1], p2 = RandomChoice[pyth2]},
      rule = Join[{ur -> p1[[1]]/p1[[2]], ux -> p2[[2]]/p2[[1]], uy -> p2[[3]]/p2[[1]]}, rule]];
    ok = ok && AllTrue[Flatten[{Together[e /. rule]}], # === 0 &], {i, k}];
  {ok, "exact at " <> ToString[k] <> " random rational points"}];

matsJ[sysA_, st_] := toJ[{sysA["A"], sysA["src"]}, st];
noDel = {dpp -> 0, dPP -> 0};
Get[FileNameJoin[{DERIV, "julia_source_reader.wl"}]];

(* ================================================================ G-reduce == *)
Print[""]; Print["================ G-reduce ================  (", stamp[], ")"];
jl2wl = "
function wl(x, env, mp)
    if x isa Integer
        return x < 0 ? \"($(x))\" : string(x)
    elseif x isa AbstractFloat
        isinteger(x) && return wl(Int(x), env, mp)
        r = rationalize(BigInt, x; tol = 0)
        return \"($(numerator(r))/$(denominator(r)))\"
    elseif x isa Rational
        return \"($(numerator(x))/$(denominator(x)))\"
    elseif x isa Symbol
        s = string(x)
        haskey(env, s) && return \"(\" * env[s] * \")\"
        haskey(mp, s) && return mp[s]
        return s
    elseif x isa Expr
        if x.head == :ref
            s = string(x); haskey(mp, s) && return mp[s]; error(\"unmapped ref $s\")
        elseif x.head == :call
            f = x.args[1]; a = x.args[2:end]
            f == :+ && return \"(\" * join([wl(y, env, mp) for y in a], \" + \") * \")\"
            if f == :-
                length(a) == 1 && return \"(-\" * wl(a[1], env, mp) * \")\"
                return \"(\" * wl(a[1], env, mp) * \" - \" * wl(a[2], env, mp) * \")\"
            end
            f == :* && return \"(\" * join([wl(y, env, mp) for y in a], \" * \") * \")\"
            (f == :/ || f == ://) && return \"((\" * wl(a[1], env, mp) * \")/(\" * wl(a[2], env, mp) * \"))\"
            f == :^ && return \"Power[\" * wl(a[1], env, mp) * \", \" * wl(a[2], env, mp) * \"]\"
            f == :sqrt && return \"Sqrt[\" * wl(a[1], env, mp) * \"]\"
            error(\"unknown call $f\")
        elseif x.head == :block
            return wl(filter(y -> !(y isa LineNumberNode), x.args)[end], env, mp)
        end
        error(\"unknown expr $(x.head)\")
    end
    error(\"unknown leaf $(typeof(x))\")
end
function smat(call, env, mp)
    dims = call.args[1].args[2:end]; a = call.args[2:end]
    (length(a) == 1 && a[1] isa Expr && a[1].head == :tuple) && (a = a[1].args)
    v = [wl(y, env, mp) for y in a]
    if length(dims) == 2
        n, m = dims
        return \"{\" * join([\"{\" * join([v[(j-1)*n+i] for j in 1:m], \", \") * \"}\" for i in 1:n], \",\\n \") * \"}\"
    end
    return \"{\" * join(v, \",\\n \") * \"}\"
end
function findfunc(ex, name)
    if ex isa Expr
        if ex.head == :function
            sig = ex.args[1]
            fn = sig.head == :call ? sig.args[1] : nothing
            fn == Symbol(name) && return ex
        end
        for a in ex.args
            r = findfunc(a, name); r === nothing || return r
        end
    end
    nothing
end
function convert_func(file, name, mp)
    f = findfunc(Meta.parseall(read(file, String)), name)
    f === nothing && error(\"function $name not found in $file\")
    env = Dict{String,String}(); out = Dict{String,String}()
    for st in f.args[2].args
        (st isa Expr && st.head == :(=)) || continue
        lhs, rhs = st.args
        if lhs isa Expr && lhs.head == :tuple
            for (l, v) in zip(lhs.args, rhs.args); env[string(l)] = wl(v, env, mp); end
        elseif rhs isa Expr && rhs.head == :call && rhs.args[1] isa Expr && rhs.args[1].head == :curly
            out[string(lhs)] = smat(rhs, env, mp)
        else
            env[string(lhs)] = wl(rhs, env, mp)
        end
    end
    out
end
io = IOBuffer(); i = 2
while i <= length(ARGS)
    mp = Dict{String,String}()
    for kv in split(ARGS[i+3], \";\"; keepempty = false)
        k, v = split(kv, \"=>\"); mp[strip(k)] = strip(v)
    end
    o = convert_func(ARGS[i], ARGS[i+1], mp)
    println(io, \"jlref[\\\"\", ARGS[i+2], \"\\\"] = <|\")
    println(io, join([\"  \\\"$k\\\" -> \" * o[k] for k in sort(collect(keys(o)))], \",\\n\"))
    println(io, \"|>;\")
    global i += 4
end
write(ARGS[1], String(take!(io)))
";
jlDir = FileNameJoin[{outDir, "julia_helpers"}];
If[!DirectoryQ[jlDir], CreateDirectory[jlDir]];
Export[FileNameJoin[{jlDir, "jl2wl.jl"}], jl2wl, "Text"];
refFile = FileNameJoin[{jlDir, "jlref.wl"}];
jlref = <||>;
convRun = RunProcess[{"julia", "--startup-file=no", FileNameJoin[{jlDir, "jl2wl.jl"}], refFile,
   FileNameJoin[{matDir, "1d_viscous_cilindrical.jl"}], "one_d_viscous_matrix", "S2cyl_shipped",
     "u[1]=>T;u[2]=>ur;u[3]=>piphi;u[4]=>pieta;u[5]=>PiB;R=>rr;p=>pres",
   FileNameJoin[{matDir, "viscous_generated.jl"}], "one_d_viscous_matrix_derived", "S2cyl_gen",
     "u[1]=>T;u[2]=>ur;u[3]=>piphi;u[4]=>pieta;u[5]=>PiB;R=>rr",
   FileNameJoin[{matDir, "viscous_generated.jl"}], "two_d_viscous_matrix_derived", "S2car_gen",
     "u[1]=>T;u[2]=>ux;u[3]=>uy;u[4]=>piyy;u[5]=>pieta;u[6]=>pixy;u[7]=>PiB",
   FileNameJoin[{matDir, "HQ_2p1d_BG_consistent_generated.jl"}], "hq2d_matrices_consistent", "S3car_gen",
     "nuxF=>nux;nuyF=>nuy;dn_da=>dnda;dn_dT=>dndT;hp=>dhdT"}];
gate["R0-jl2wl", convRun["ExitCode"] === 0 && FileExistsQ[refFile],
  "Julia -> WL conversion of the shipped/generated matrices (" <> ToString[convRun["ExitCode"]] <> ")"];
If[convRun["ExitCode"] =!= 0, Print[convRun["StandardError"]]];
If[FileExistsQ[refFile], Get[refFile]];

cmpMats[mine_, ref_, label_, id_] := Module[{d, z},
  d = Flatten[{mine[[1]] - ref[[1]], mine[[2]] - ref[[2]]}];
  z = zeroJ[d];
  gate[id, z[[1]], label <> ", entry by entry (" <> ToString[Length[d]] <> " entries, " <> z[[2]] <> ")"]];

If[KeyExistsQ[results, "1p1d"] && KeyExistsQ[jlref, "S2cyl_shipped"],
  Module[{st = stateCYL, m = matsJ[results["1p1d"]["sys"]["S2"], stateCYL] /. noDel, r},
    r = jlref["S2cyl_shipped"];
    cmpMats[{m[[1]], m[[2]]}, {{r["At"], r["Ax"]}, r["source"]},
      "S2 1+1D generic EOS, dpp = dPP = 0 == SHIPPED one_d_viscous_matrix (1d_viscous_cilindrical.jl)", "R-S2-1p1d-shipped"];
    r = jlref["S2cyl_gen"];
    cmpMats[{m[[1]], m[[2]]}, {{r["At"], r["Ax"]}, r["source"]},
      "S2 1+1D, dpp = dPP = 0 == one_d_viscous_matrix_derived (derive_2p1d_viscous.wls output)", "R-S2-1p1d-generator"];
  ]];
If[KeyExistsQ[results, "2p1d"] && KeyExistsQ[jlref, "S2car_gen"],
  Module[{m = matsJ[results["2p1d"]["sys"]["S2"], stateCAR] /. noDel, r},
    r = jlref["S2car_gen"];
    cmpMats[{m[[1]], m[[2]]}, {{r["At"], r["Ax"], r["Ay"]}, r["source"]},
      "S2 2+1D, dpp = dPP = 0 == two_d_viscous_matrix_derived (derive_2p1d_viscous.wls output)", "R-S2-2p1d-generator"];
    m = matsJ[results["2p1d"]["sys"]["S3"], stateCAR];
    r = jlref["S3car_gen"];
    cmpMats[{m[[1]], m[[2]]}, {{r["At"], r["Ax"], r["Ay"]}, r["src"]},
      "S3 2+1D == hq2d_matrices_consistent (derive_hq_consistent_2p1d.wls output)", "R-S3-2p1d-generator"];
  ]];

(* ---- the 7-field 1-D :MIS kernel: its medium block carries dpp = 4/3 tauS ---- *)
If[KeyExistsQ[results, "1p1d"],
  Module[{m = matsJ[results["1p1d"]["sys"]["S2"], stateCYL], sh, blk, d, ctl},
    sh = Values[juliaMatrices[FileNameJoin[{matDir, "1d_viscous_HQ_cilindrical_fugacity.jl"}], "one_d_viscous_matrix_fugacity"]] /.
      {jl["u"][1] -> T, jl["u"][2] -> ur, jl["u"][3] -> piphi, jl["u"][4] -> pieta, jl["u"][5] -> PiB,
       jl["tau"] -> tau, jl["R"] -> rr, jl["p"] -> pres, jl["dtp"] -> dtp, jl["dtdtp"] -> dtdtp, jl["zeta"] -> zeta,
       jl["visc"] -> visc, jl["tauS"] -> tauS, jl["tauB"] -> tauB, jl["dmp"] -> 0, jl["dtdmp"] -> 0, jl["dmdmp"] -> 0};
    blk = {{sh[[1, 1 ;; 5, 1 ;; 5]], sh[[2, 1 ;; 5, 1 ;; 5]]}, sh[[3, 1 ;; 5]]};
    cmpMats[m /. {dpp -> 4 tauS/3, dPP -> 0}, blk,
      "S2 1+1D with dpp = 4/3 tauS, dPP = 0 == medium block of one_d_viscous_matrix_fugacity (:MIS, read from the .jl)",
      "R-S2-1p1d-MIS"];
    d = Flatten[{(m /. noDel)[[1]] - blk[[1]], (m /. noDel)[[2]] - blk[[2]]}];
    gate["R-S2-1p1d-MIS-ctl", !zeroJ[d, 2, 10][[1]], "control: dpp = 0 does NOT reproduce the :MIS medium block (the two 1-D media differ by delta_pipi)"];
  ]];

(* ---- the reference CSVs of the 1+1D derivations (FLUIDUM_EXTERNAL_REFPOINTS) ---- *)
readCSV[f_] := Module[{d = Import[f, "CSV"]}, {d[[1]], Rest[d]}];
relErr[a_, b_] := Abs[a - b]/Max[1, Abs[b]];
(* Boltzmann charm, degeneracy gq: the closed forms of the 1+1D second-moment derivation *)
boltz[Tv_, av_, mqv_, gq_] := Module[{TT0, a0, nn, hh, z, i31, aq, bq, cq},
  z = mqv/TT0;
  nn = gq Exp[a0] mqv^2 TT0 BesselK[2, z]/(2 Pi^2);
  hh = mqv BesselK[3, z]/BesselK[2, z];
  i31 = nn TT0 hh; aq = mqv^2 + 5 TT0 hh; bq = TT0 hh; cq = mqv BesselK[4, z]/BesselK[3, z];
  N[{n -> nn, dndT -> D[nn, TT0], dnda -> D[nn, a0], h -> hh, dhdT -> D[hh, TT0],
     I31 -> i31, dI31dT -> D[i31, TT0], dI31da -> D[i31, a0], aQ -> aq, bQ -> bq,
     dbQdT -> D[bq, TT0], cQ -> cq, dcQdT -> D[cq, TT0]} /. {TT0 -> Tv, a0 -> av}, 30]];

If[KeyExistsQ[results, "1p1d"],
  Module[{st = stateCYL, S3 = matsJ[results["1p1d"]["sys"]["S3"], stateCYL],
      S4 = matsJ[results["1p1d"]["sys"]["S4"], stateCYL], hdr, rows, err, ix, vals, rep, mine, ref,
      f1 = FileNameJoin[{csvDir, "hq_consistent_refpoints.csv"}],
      f2 = FileNameJoin[{csvDir, "hq_consistent_m2_refpoints.csv"}],
      f3 = FileNameJoin[{csvDir, "hq_cm_coupling_refpoints.csv"}]},
    If[csvDir === "" || !AllTrue[{f1, f2, f3}, FileExistsQ],
      Print["  [SKIP] R-S3/R-S4/R-S4c: the 1+1D reference CSVs are not available (FLUIDUM_EXTERNAL_REFPOINTS)"],
    (* S3 1-D vs the 1+1D derivation of the consistent first moment *)
    {hdr, rows} = readCSV[f1];
    ix[s_] := First@First@Position[hdr, s];
    err = Max @ Table[
      rep = {tau -> r[[ix["tau"]]], rr -> r[[ix["rr"]]], ur -> r[[ix["ur"]]], T -> r[[ix["T"]]],
             al -> r[[ix["alpha"]]], nur -> r[[ix["nur"]]], n -> r[[ix["n"]]], dndT -> r[[ix["n10"]]],
             dnda -> r[[ix["n01"]]], h -> r[[ix["h"]]], dhdT -> r[[ix["hp"]]], Ds -> r[[ix["Ds"]]],
             dtT -> r[[ix["dtT"]]], drT -> r[[ix["drT"]]], dtur -> r[[ix["dtur"]]], drur -> r[[ix["drur"]]]};
      mine = Flatten[{S3[[1, 1]], S3[[1, 2]], S3[[2]]} /. rep];
      ref = Flatten[{Table[r[[ix["AtC" <> ToString[i] <> ToString[j]]]], {i, 2}, {j, 2}],
                     Table[r[[ix["ArC" <> ToString[i] <> ToString[j]]]], {i, 2}, {j, 2}],
                     {r[[ix["SrcC1"]]], r[[ix["SrcC2"]]]}}];
      Max[MapThread[relErr, {mine, ref}]], {r, rows}];
    gate["R-S3-1p1d-generator", err < 1*^-11,
      "S3 1+1D == the 1+1D consistent first moment (hq_consistent_refpoints.csv, " <> ToString[Length[rows]] <>
        " pts x 10 entries), max rel err " <> sci[err]];
    (* S4 1-D rows 3-5 vs the 1+1D derivation of the consistent second moment (Boltzmann, g = 1) *)
    {hdr, rows} = readCSV[f2];
    err = Max @ Table[
      rep = Join[{tau -> r[[ix["tau"]]], rr -> r[[ix["rr"]]], ur -> r[[ix["ur"]]], T -> r[[ix["T"]]],
             al -> r[[ix["alpha"]]], nur -> r[[ix["nur"]]], pl -> r[[ix["pl"]]], pphi -> r[[ix["pph"]]],
             PiQ -> r[[ix["PiQ"]]], dtT -> r[[ix["dtT"]]], drT -> r[[ix["drT"]]], dtur -> r[[ix["dtur"]]],
             drur -> r[[ix["drur"]]], Ds -> r[[ix["Ds"]]], c -> 0},
             boltz[r[[ix["T"]]], r[[ix["alpha"]]], r[[ix["mq"]]], 1]];
      mine = N[Flatten[{S4[[1, 1]][[3 ;; 5]], S4[[1, 2]][[3 ;; 5]], S4[[2]][[3 ;; 5]]} /. rep]];
      ref = Flatten[{Table[r[[ix["At3_" <> ToString[i] <> "_" <> ToString[j]]]], {i, 3}, {j, 5}],
                     Table[r[[ix["Ax3_" <> ToString[i] <> "_" <> ToString[j]]]], {i, 3}, {j, 5}],
                     Table[r[[ix["Src3_" <> ToString[i]]]], {i, 3}]}];
      Max[MapThread[relErr, {mine, ref}]], {r, rows}];
    gate["R-S4-1p1d-generator", err < 1*^-11,
      "S4 1+1D rows 3-5 (Boltzmann) == the 1+1D consistent second moment (hq_consistent_m2_refpoints.csv, " <>
        ToString[Length[rows]] <> " pts x 33 entries), max rel err " <> sci[err]];
    (* the c-coupling of row 2 vs the 1+1D derivation of the c_M coupling (per unit c_M = c Ds/T) *)
    {hdr, rows} = readCSV[f3];
    err = Max @ Table[
      rep = {tau -> r[[ix["tau"]]], rr -> r[[ix["rr"]]], ur -> r[[ix["ur"]]], pl -> r[[ix["pr"]]],
             pphi -> r[[ix["pp"]]], PiQ -> r[[ix["PQ"]]], dtur -> r[[ix["dtur"]]], drur -> r[[ix["drur"]]]};
      mine = Flatten[{D[S4[[1, 1]][[2, 3 ;; 5]], c], D[S4[[1, 2]][[2, 3 ;; 5]], c], D[S4[[2]][[2]], c]}/(Ds/T) /. rep];
      ref = {r[[ix["At1"]]], r[[ix["At2"]]], r[[ix["At3"]]], r[[ix["Ar1"]]], r[[ix["Ar2"]]], r[[ix["Ar3"]]], r[[ix["Src"]]]};
      Max[MapThread[relErr, {Simplify[mine], ref}]], {r, rows}];
    gate["R-S4c-1p1d-generator", err < 1*^-11,
      "S4/S5 1+1D c-coupling of the nu row / (c Ds/T) == the 1+1D c_M coupling (hq_cm_coupling_refpoints.csv, " <>
        ToString[Length[rows]] <> " pts x 7 entries), max rel err " <> sci[err]]];
    (* and the full 2-row first-moment block of S4 at c = 0 is S3 *)
    gate["R-S4S3-1p1d", Simplify[{S4[[1]][[All, 1 ;; 2, 1 ;; 2]], S4[[2]][[1 ;; 2]]} /. c -> 0] ===
        Simplify[{S3[[1]], S3[[2]]}],
      "S4 rows 1-2 at c = 0 == S3 (so the first-moment block of S4/S5 is the gated consistent closure)"];
  ]];

(* ---- the 2-D second moment vs Fluidum's hq2d_m2_rows (Julia, numeric) ---- *)
If[KeyExistsQ[results, "2p1d"],
  Module[{st = stateCAR, S4 = matsJ[results["2p1d"]["sys"]["S4"], stateCAR], vars, inNames, rowsOut = {},
      csvf = FileNameJoin[{jlDir, "m2_2p1d_states.csv"}], jls, run, out, gq = 6, mqv = 3/2},
    SeedRandom[20261013];
    Do[Module[{Tv = RandomReal[{0.15, 0.45}, WorkingPrecision -> 30], av = RandomReal[{-2, 2}, WorkingPrecision -> 30],
        th, rep, dphi, row, val, pick},
      th = boltz[Tv, av, mqv, gq];
      rep = Join[{tau -> RandomReal[{0.6, 10}, WorkingPrecision -> 30], T -> Tv, al -> av,
          ux -> RandomReal[{-1.2, 1.2}, WorkingPrecision -> 30], uy -> RandomReal[{-1.2, 1.2}, WorkingPrecision -> 30],
          nux -> RandomReal[{-0.02, 0.02}, WorkingPrecision -> 30], nuy -> RandomReal[{-0.02, 0.02}, WorkingPrecision -> 30],
          pQxx -> RandomReal[{-0.01, 0.01}, WorkingPrecision -> 30], pQxy -> RandomReal[{-0.01, 0.01}, WorkingPrecision -> 30],
          pQyy -> RandomReal[{-0.01, 0.01}, WorkingPrecision -> 30], PiQ -> RandomReal[{-0.01, 0.01}, WorkingPrecision -> 30],
          Ds -> RandomReal[{0.05, 1.2}, WorkingPrecision -> 30], c -> 0},
         Table[s -> RandomReal[{-0.3, 0.3}, WorkingPrecision -> 30],
           {s, {dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy}}], th];
      (* derivatives of alpha, nu^x, nu^y (columns 1-3); the m2 fields' own derivatives are set to 0,
         which removes exactly the diagonal advection tau_M u.d that hq2d_m2_rows leaves to A_t/A_x/A_y *)
      dphi = Table[Join[RandomReal[{-0.5, 0.5}, 3, WorkingPrecision -> 30], {0, 0, 0, 0}], {3}];
      Do[
        val = N[(Sum[S4[[1, k]] . dphi[[k]], {k, 3}] + S4[[2]]) /. Join[{wv -> wvv}, rep], 25];
        AppendTo[rowsOut, Join[{wvv},
          {tau, T, al, ux, uy, nux, nuy, pQxx, pQxy, pQyy, PiQ, Ds, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy,
           n, h, dhdT, cQ} /. rep,
          Flatten[dphi[[All, 1 ;; 3]]], val[[4 ;; 7]]]], {wvv, {0, 1}}]], {12}];
    Export[csvf, Prepend[N[rowsOut], Join[{"wv", "tau", "T", "al", "ux", "uy", "nux", "nuy", "pQxx", "pQxy",
        "pQyy", "PiQ", "Ds", "dtT", "dxT", "dyT", "dtux", "dxux", "dyux", "dtuy", "dxuy", "dyuy", "n", "h", "hp", "cQ",
        "dta", "dtnx", "dtny", "dxa", "dxnx", "dxny", "dya", "dynx", "dyny", "Sxx", "Sxy", "Syy", "SB"}]], "CSV"];
    jls = "
using Bessels, StaticArrays
module M2
using Bessels, StaticArrays
include(ARGS[1])
end
lines = readlines(ARGS[2]); hdr = split(lines[1], ',') .|> s -> strip(s, '\"')
col(s) = findfirst(==(s), hdr)
worst = 0.0; worst0 = 0.0; worst1 = 0.0
for ln in lines[2:end]
    v = parse.(Float64, split(ln, ','))
    g(s) = v[col(s)]
    T = g(\"T\"); Ds = g(\"Ds\"); h = g(\"h\"); hp = g(\"hp\"); n = g(\"n\")
    mq = 1.5
    taun = Ds * h / T; tauM = Ds * g(\"cQ\") / (2T); etaM = T * taun / 2
    ux = g(\"ux\"); uy = g(\"uy\"); uT = sqrt(1 + ux^2 + uy^2)
    pet = M2.hq2d_piQ_eta(ux, uy, uT, g(\"pQxx\"), g(\"pQxy\"), g(\"pQyy\"))
    X = (g(\"al\"), g(\"nux\"), g(\"nuy\"), g(\"pQxx\"), g(\"pQxy\"), g(\"pQyy\"), pet, g(\"PiQ\"))
    r = M2.hq2d_m2_rows(X, g(\"tau\"), T, ux, uy, g(\"dtT\"), g(\"dxT\"), g(\"dyT\"),
            g(\"dtux\"), g(\"dxux\"), g(\"dyux\"), g(\"dtuy\"), g(\"dxuy\"), g(\"dyuy\"),
            g(\"dxa\"), g(\"dya\"), g(\"dta\"), g(\"dtnx\"), g(\"dxnx\"), g(\"dynx\"), g(\"dtny\"), g(\"dxny\"), g(\"dyny\"),
            n, taun, Ds, h, hp, tauM, etaM, mq; projector = true, vorticity = g(\"wv\") == 1)
    SJ = -(tauM * uT) .* collect(r)
    SW = [g(\"Sxx\"), g(\"Sxy\"), g(\"Syy\"), g(\"SB\")]
    e = maximum(abs.(SJ .- SW) ./ max.(abs.(SW), 1e-3 * maximum(abs.(SW))))
    global worst = max(worst, e)
    g(\"wv\") == 0 ? (global worst0 = max(worst0, e)) : (global worst1 = max(worst1, e))
end
println(\"WORST \", worst, \" \", worst0, \" \", worst1)
";
    Export[FileNameJoin[{jlDir, "check_m2_2p1d.jl"}], jls, "Text"];
    run = RunProcess[{"julia", "--startup-file=no", "--project=" <> juliaProject,
        FileNameJoin[{jlDir, "check_m2_2p1d.jl"}], FileNameJoin[{matDir, "HQ_2p1d_BG_m2.jl"}], csvf}];
    out = StringCases[run["StandardOutput"], "WORST " ~~ a : (Except[WhitespaceCharacter] ..) ~~ " " ~~ b : (Except[WhitespaceCharacter] ..) ~~ " " ~~ cc : (Except[WhitespaceCharacter] ..) :>
        ToExpression /@ StringReplace[{a, b, cc}, "e" -> "*^"]];
    If[out === {}, Print[run["StandardError"]]];
    gate["R-S4-2p1d-fluidum", out =!= {} && out[[1, 1]] < 1*^-9,
      "S4 2+1D second-moment rows == Fluidum hq2d_m2_rows (projector = true), Julia, 24 states: " <>
        If[out === {}, "julia failed",
          "wv=0 <-> vorticity=false: " <> sci[out[[1, 2]]] <>
          ", wv=1 <-> vorticity=true: " <> sci[out[[1, 3]]]]];
  ]];

(* ================================================================= G-chart == *)
(* built with Function @@ at run time: a literal Function[{x}, Evaluate[..]] inside a Module is
   lexically renamed (x -> x$) while the evaluated body keeps the global x -- a silent constant *)
mkF3[body_] := Function @@ {{tq, xq, yq}, body};
Print[""]; Print["================ G-chart ================  (", stamp[], ")"];
If[KeyExistsQ[results, "2p1d"] && KeyExistsQ[results, "1p1d"],
  Module[{Rk = results["2p1d"]["R"], Rc = results["1p1d"]["R"], sub, rho, ax, cy, kk, pairs, z, all, Wax},
    rho = Sqrt[xq^2 + yq^2];
    With[{ur1 = urc[tq, rho], pff1 = pffc[tq, rho], pee1 = peec[tq, rho], pl1 = plc[tq, rho], pph1 = pphc[tq, rho]},
      With[{prr = -(1 + ur1^2) (pff1 + pee1), qrr = pl1 (1 + ur1^2)},
        sub = {
          TTk -> mkF3[(TTc[tq, rho])],
          uxk -> mkF3[(ur1 xq/rho)], uyk -> mkF3[(ur1 yq/rho)],
          pyyk -> mkF3[(prr yq^2/rho^2 + pff1 xq^2/rho^2)],
          pxyk -> mkF3[((prr - pff1) xq yq/rho^2)],
          peek -> mkF3[(pee1)], PiBk -> mkF3[(PiBc[tq, rho])],
          alk -> mkF3[(alc[tq, rho])],
          nuxk -> mkF3[(nurc[tq, rho] xq/rho)],
          nuyk -> mkF3[(nurc[tq, rho] yq/rho)],
          qxxk -> mkF3[(qrr xq^2/rho^2 + pph1 yq^2/rho^2)],
          qxyk -> mkF3[((qrr - pph1) xq yq/rho^2)],
          qyyk -> mkF3[(qrr yq^2/rho^2 + pph1 xq^2/rho^2)],
          PiQk -> mkF3[(PiQc[tq, rho])]}]];
    ax[e_] := (((e /. sub) /. yy -> 0) /. Power[xx^2, p_] :> xx^(2 p)) /. xx -> rr;
    cy = Join[Rc["medM"] + Rc["medQ"], Rc["ch"], Rc["m2"]];
    kk = Join[Rk["medM"] + Rk["medQ"], Rk["ch"], Rk["m2"]];
    Print["  axisymmetric substitution ... (", stamp[], ")"];
    all = ax /@ kk;
    Print["  done (", stamp[], ")"];
    (* row map 2+1D -> 1+1D (on the +x axis, x = r) *)
    pairs = {
      {1, cy[[1]], "tau-momentum"}, {2, cy[[2]], "x-momentum == r-momentum"}, {3, 0, "y-momentum == 0"},
      {4, rr^2 cy[[3]], "E^{yy} == r^2 E^{phiphi}"}, {5, cy[[4]], "E^{etaeta}"}, {6, 0, "E^{xy} == 0"},
      {7, cy[[5]], "bulk"},
      {8, cy[[6]], "charm number"}, {9, cy[[7]], "nu^x row == nu^r row"}, {10, 0, "nu^y row == 0"},
      {11, (1 + urc[tau, rr]^2) cy[[8]], "piQ^{xx} row == (u^tau)^2 x (l l) row"}, {12, 0, "piQ^{xy} row == 0"},
      {13, cy[[9]], "piQ^{yy} row == (phi phi) row"}, {14, cy[[10]], "trace row"}};
    Do[
      z = zeroTest[all[[p[[1]]]] - p[[2]], stateCYL["med"] ~Join~ stateCYL["ch"] ~Join~ stateCYL["m2"], {urc}, 30];
      gate["X-2to1-" <> ToString[p[[1]]], z[[1]], "S5 2+1D axisymmetric on the x axis: " <> p[[3]] <> " (" <> z[[2]] <> ")"],
      {p, pairs}];
    z = zeroTest[all[[1]] - cy[[2]], stateCYL["med"] ~Join~ stateCYL["ch"] ~Join~ stateCYL["m2"], {urc}, 5];
    gate["X-ctrl", !z[[1]], "control: the same test REJECTS a wrong pairing (2+1D tau row vs 1+1D r row)"];
    Wax = ax /@ Flatten[results["2p1d"]["B"]["W"]];
    gate["X-vort-axis", Simplify[Wax] === ConstantArray[0, 16],
      "the corotation coupling W vanishes identically on an axisymmetric state (so G-chart is blind to wv, as documented)"];
  ]];
If[KeyExistsQ[results, "1p1d"] && KeyExistsQ[results, "0p1d"],
  Module[{Rc = results["1p1d"]["R"], Rb = results["0p1d"]["R"], sub, cy, bj, all, pairs, z},
    sub = {TTc -> Function[{tq, rq}, TTb[tq]], urc -> Function[{tq, rq}, 0],
      pffc -> Function[{tq, rq}, -peeb[tq]/2], peec -> Function[{tq, rq}, peeb[tq]],
      PiBc -> Function[{tq, rq}, PiBb[tq]], alc -> Function[{tq, rq}, alb[tq]], nurc -> Function[{tq, rq}, 0],
      plc -> Function[{tq, rq}, -qeb[tq]/2], pphc -> Function[{tq, rq}, -qeb[tq]/2],
      PiQc -> Function[{tq, rq}, PiQb[tq]]};
    cy = Join[Rc["medM"] + Rc["medQ"], Rc["ch"], Rc["m2"]];
    bj = Join[Rb["medM"] + Rb["medQ"], Rb["ch"], Rb["m2"]];
    all = cy /. sub;
    pairs = {{1, bj[[1]], "energy"}, {2, 0, "r-momentum == 0"}, {3, -(tau^2/(2 rr^2)) bj[[2]], "E^{phiphi} == -(tau^2/2r^2) E^{etaeta}"},
      {4, bj[[2]], "E^{etaeta}"}, {5, bj[[3]], "bulk"}, {6, bj[[4]], "charm number"}, {7, 0, "nu^r row == 0"},
      {8, -bj[[5]]/2, "(l l) row == -1/2 (eta eta) row"}, {9, -bj[[5]]/2, "(phi phi) row == -1/2 (eta eta) row"},
      {10, bj[[6]], "trace row"}};
    Do[
      z = zeroTest[all[[p[[1]]]] - p[[2]], stateBJ["med"] ~Join~ stateBJ["ch"] ~Join~ stateBJ["m2"], {}, 30];
      gate["X-1to0-" <> ToString[p[[1]]], z[[1]], "S5 1+1D at u^r = 0, no r-dependence: " <> p[[3]] <> " (" <> z[[2]] <> ")"],
      {p, pairs}];
  ]];

(* =============================================================== G-bjorken == *)
Print[""]; Print["================ G-bjorken ================  (", stamp[], ")"];
If[KeyExistsQ[results, "0p1d"],
  Module[{sy = results["0p1d"]["sys"], st = stateBJ, conf, sol, Tf = TTb[tau], pe = peeb[tau], Pb = PiBb[tau], ok,
      s5, rate},
    conf = {pp -> Function[x, aSB x^4]};
    sol = -sy["S1"]["src"][[1]]/sy["S1"]["A"][[1, 1, 1]];
    gate["B-ideal-conf", Simplify[(sol /. conf) + Tf/(3 tau)] === 0, "conformal ideal Bjorken: dT/dtau == -T/(3 tau)"];
    gate["B-ideal-gen", Simplify[sol + pp'[Tf]/(pp''[Tf] tau)] === 0,
      "generic-EOS ideal Bjorken: dT/dtau == -P'/(P'' tau) == -c_s^2 T/tau, c_s^2 = P'/(T P'')"];
    sol = -Inverse[sy["S2"]["A"][[1]]] . sy["S2"]["src"];
    gate["B-IS-gen", Simplify[sol - {-(Tf pp'[Tf] + Pb + pe)/(tau Tf pp''[Tf]), -(pe + 4 visc/(3 tau) + dpp pe/tau)/tauS,
         -(Pb + zeta/tau + dPP Pb/tau)/tauB}] === {0, 0, 0},
      "generic IS Bjorken: T' = -(T P' + Pi + pi)/(tau T P''), tauS pi' = -pi - 4 visc/(3 tau) - dpp pi/tau, tauB Pi' = -Pi - zeta/tau - dPP Pi/tau"];
    sol = Simplify[sol /. conf /. PiBb -> Function[x, 0]];
    gate["B-IS-conf", Simplify[sol[[1 ;; 2]] - {-Tf/(3 tau) - pe/(12 aSB Tf^3 tau), -(pe + 4 visc/(3 tau) + dpp pe/tau)/tauS}] === {0, 0},
      "conformal IS Bjorken: dT/dtau = -T/(3tau) - pi/(12 a T^3 tau), tauS dpi/dtau = -pi - 4 visc/(3 tau) - dpp pi/tau  [dpp = 4/3 tauS: the DNMR delta_pipi form; no tau_pipi, lambda_1]"];
    gate["B-NS-conf", Simplify[(sol[[1]] /. peeb -> Function[x, -4 visc/(3 x)]) - (-Tf/(3 tau) + (4/9) (visc/(4 aSB Tf^3))/tau^2)] === 0,
      "first-order (NS) limit pi = -4 visc/(3 tau): dT/dtau = -T/(3tau) + (4/9)(eta/s)/tau^2, s = 4 a T^3"];
    (* the two-way energy row in Bjorken *)
    s5 = results["0p1d"]["R"]["medM"][[1]] + results["0p1d"]["R"]["medQ"][[1]];
    gate["B-twoway-energy", Simplify[s5 - (D[Tf pp'[Tf] - pp[Tf] + nQ[Tf, alb[tau]] hQ[Tf] - nQ[Tf, alb[tau]] Tf, tau] +
         (Tf pp'[Tf] + Pb + pe + nQ[Tf, alb[tau]] hQ[Tf] + csw (qeb[tau] + PiQb[tau]))/tau)] === 0,
      "S5 Bjorken energy row: d_tau(e + e_Q) + (e+P+Pi+pi + n h + c(piQ^eta_eta + Pi_Q))/tau"];
    (* charm number in Bjorken: d_tau(n tau) = 0 *)
    gate["B-charm-number", Simplify[results["0p1d"]["R"]["ch"][[1]] - D[tau nQ[Tf, alb[tau]], tau]/tau] === 0,
      "Bjorken charm number: (1/tau) d_tau (n tau) = 0"];
  ]];

(* ================================================================ G-twoway == *)
Print[""]; Print["================ G-twoway ================  (", stamp[], ")"];
twowayLimits[st_, res_] := Module[{nm = st["name"], B = res["B"], R = res["R"], rest, drag, nuh, tr, ch,
    mom, solv, Tf, dnu, ok},
  Tf = TT0[tau];
  Switch[nm,
   "1p1d",
     rest = {TTc -> Function[{tq, rq}, TT0[tq]], urc -> Function[{tq, rq}, 0], alc -> Function[{tq, rq}, al0[tq]],
             nurc -> Function[{tq, rq}, 0], plc -> Function[{tq, rq}, 0], pphc -> Function[{tq, rq}, 0],
             PiQc -> Function[{tq, rq}, 0]};
     drag = {TTc -> Function[{tq, rq}, TT0[tq]], urc -> Function[{tq, rq}, 0], alc -> Function[{tq, rq}, al0[tq]],
             nurc -> Function[{tq, rq}, nu0[tq]], plc -> Function[{tq, rq}, 0], pphc -> Function[{tq, rq}, 0],
             PiQc -> Function[{tq, rq}, 0]};
     mom = {2}; ch = {2},
   "2p1d",
     rest = {TTk -> Function[{tq, xq, yq}, TT0[tq]], uxk -> Function[{tq, xq, yq}, 0], uyk -> Function[{tq, xq, yq}, 0],
             alk -> Function[{tq, xq, yq}, al0[tq]], nuxk -> Function[{tq, xq, yq}, 0], nuyk -> Function[{tq, xq, yq}, 0],
             qxxk -> Function[{tq, xq, yq}, 0], qxyk -> Function[{tq, xq, yq}, 0], qyyk -> Function[{tq, xq, yq}, 0],
             PiQk -> Function[{tq, xq, yq}, 0]};
     (* drag: uniform current in a static uniform medium, with a homogeneous second moment kept on
        (it must not spoil the balance at c = 1) *)
     drag = {TTk -> Function[{tq, xq, yq}, TT0[tq]], uxk -> Function[{tq, xq, yq}, 0], uyk -> Function[{tq, xq, yq}, 0],
             alk -> Function[{tq, xq, yq}, al0[tq]], nuxk -> Function[{tq, xq, yq}, nu0[tq]],
             nuyk -> Function[{tq, xq, yq}, nu1[tq]],
             qxxk -> Function[{tq, xq, yq}, q0[tq]], qxyk -> Function[{tq, xq, yq}, q1[tq]],
             qyyk -> Function[{tq, xq, yq}, q2[tq]], PiQk -> Function[{tq, xq, yq}, q3[tq]]};
     mom = {2, 3}; ch = {2, 3}];
  (* (i) charm at rest and homogeneous, static medium: no transverse momentum exchange *)
  ok = Simplify[(B["consQ"][[mom]] /. rest)] === ConstantArray[0, Length[mom]] &&
       Simplify[(R["ch"][[ch]] /. rest)] === ConstantArray[0, Length[ch]];
  gate["T-rest-" <> nm, ok,
    "nu = 0, dM = 0, grad alpha = 0, medium at rest: transverse momentum exchange Delta.Nabla.T_Q == 0 and the nu row is satisfied (c symbolic)"];
  (* (ii) stationary drag: on the charm row's shell the charm loses (T/Ds) nu, the medium gains it *)
  dnu = If[nm === "1p1d", {nu0'[tau]}, {nu0'[tau], nu1'[tau]}];
  solv = Solve[Thread[(R["ch"][[ch]] /. drag) == 0], dnu][[1]];
  nuh = If[nm === "1p1d", {nu0[tau]}, {nu0[tau], nu1[tau]}];
  ok = Simplify[((B["consQ"][[mom]] /. drag) /. solv) + (Tf/Ds) nuh] === ConstantArray[0, Length[mom]];
  gate["T-drag-" <> nm, ok,
    "drag: on shell, Nabla.T_Q^i = -(T/Ds) nu^i (drag on charm) and the fluid row then gives Nabla.T_m^i = +(T/Ds) nu^i (Newton III), exact, c symbolic"];
  (* and the fluid row splits exactly into medium + charm *)
  ok = Simplify[((R["medM"][[mom]] + R["medQ"][[mom]]) - B["consM"][[mom]] - B["consQ"][[mom]])] === ConstantArray[0, Length[mom]];
  gate["T-split-" <> nm, ok, "fluid momentum rows == Nabla.T_m + Nabla.T_Q identically"];
];
If[KeyExistsQ[results, "1p1d"], twowayLimits[stateCYL, results["1p1d"]]];
If[KeyExistsQ[results, "2p1d"], twowayLimits[stateCAR, results["2p1d"]]];

(* ============================================================ char. speeds == *)
(* Principal part of S5: eigenvalues of A_t^{-1} A(n), A(n) = n_x A_x + n_y A_y, exact.
   Physical states: conformal medium (P = a T^4, a = (16 + 21 Nf/2) pi^2/90, Nf = 3, no bulk),
   eta/s = 0.16, tauS = 5 (eta/s)/T; Boltzmann charm g = 6, m = 1.5 GeV (n, h, I31, A, B, C and their
   derivatives from the closed forms); Ds = DsT/T; natural units.  Every number is rationalised
   (25 digits) so that the characteristic polynomial, its factorisation, and the rank test for
   repeated roots are exact:
     hyperbolic  <=> all roots real AND every repeated root has geometric = algebraic multiplicity.
   u is drawn from Pythagorean tuples so that u^tau is rational.  The results are reported, not checked. *)
Print[""]; Print["================ characteristic speeds of S5 (exact) ================  (", stamp[], ")"];
speedRows = {};
ratz[x_] := Rationalize[SetPrecision[x, 30], 10^-25];
charAnalyse[M_] := Module[{x, cp, fac, rts, imax, rmax, defect = 0},
  cp = CharacteristicPolynomial[M, x];
  fac = FactorList[cp];
  rts = x /. NSolve[cp == 0, x, WorkingPrecision -> 40];
  imax = Max[Abs[Im[rts]]]; rmax = Max[Abs[rts]];
  Do[If[f[[2]] > 1 && Exponent[f[[1]], x] >= 1,
      Do[defect += f[[2]] - (Length[M] - MatrixRank[M - r0 IdentityMatrix[Length[M]]]),
        {r0, x /. Solve[f[[1]] == 0, x]}]], {f, fac}];
  {rmax, imax, defect, Sort[Re[rts]]}];
speedScan[nm_] := Module[{st = If[nm === "1p1d", stateCYL, stateCAR], m, aSBv, states, Mt, Mx, My, dirs, uset, an},
  m = matsJ[results[nm]["sys"]["S5"], st];
  aSBv = (16 + 21 3/2) Pi^2/90;
  m = m /. {pres -> aSB T^4, dtp -> 4 aSB T^3, dtdtp -> 12 aSB T^2, zeta -> 0, tauB -> 1, PiB -> 0, dPP -> 0};
  dirs = If[nm === "1p1d", {{1, 0}}, {{1, 0}, {0, 1}, {3/5, 4/5}}];
  uset = If[nm === "1p1d", {{0, 0}, {3/4, 0}, {12/5, 0}}, {{0, 0}, {4/7, 4/7}, {2, 2}}];
  states = Flatten[Table[{Tv, uv, dsT, load}, {Tv, If[nm === "1p1d", {0.16, 0.25, 0.40}, {0.16, 0.40}]}, {uv, uset},
      {dsT, {0.116, 0.40}}, {load, {"phys", "heavy"}}], 3];
  Do[Module[{Tv, uv, dsT, load, av, th, rep, s0, P4},
      {Tv, uv, dsT, load} = sv;
      (* alpha: 'phys' -> e^alpha = 1 (n_c ~ 1e-4 GeV^3); 'heavy' -> e_Q/e_m = 0.3 *)
      av = If[load === "phys", 0., Log[0.3 (3 aSBv Tv^4)/((n (h - Tv)) /. boltz[Tv, 0, 1.5, 6])]];
      th = boltz[Tv, av, 1.5, 6];
      s0 = n /. th; P4 = aSBv Tv^4;
      Do[Do[Do[
        rep = Join[{T -> Tv, tau -> 2/0.1973, rr -> 3/0.1973, aSB -> aSBv,
            ux -> uv[[1]], uy -> uv[[2]], ur -> uv[[1]],
            piyy -> 0.02 P4, pieta -> -0.05 P4, pixy -> 0.01 P4, piphi -> 0.025 P4,
            al -> av, nux -> 0.05 s0, nuy -> -0.03 s0, nur -> 0.05 s0,
            pQxx -> 0.02 s0 Tv, pQxy -> 0.01 s0 Tv, pQyy -> -0.01 s0 Tv, pl -> 0.02 s0 Tv, pphi -> -0.01 s0 Tv,
            PiQ -> 0.01 s0 Tv, visc -> 0.16 4 aSBv Tv^3, tauS -> 5 0.16/Tv, dpp -> dpf 5 0.16/Tv, Ds -> dsT/Tv}, th];
        rep = Join[{ux -> uv[[1]], uy -> uv[[2]], ur -> uv[[1]], c -> cv, wv -> wq}, (#[[1]] -> ratz[N[#[[2]], 30]]) & /@ rep];
        Mt = m[[1, 1]] /. rep; Mx = m[[1, 2]] /. rep; My = If[nm === "2p1d", m[[1, 3]] /. rep, 0 Mx];
        Do[
          an = charAnalyse[Inverse[Mt] . (d[[1]] Mx + d[[2]] My)];
          AppendTo[speedRows, {nm, Tv, N[Norm[uv]], dsT, load, cv, wq, d, an[[1]], an[[2]], an[[3]], an[[4]], dpf}], {d, dirs}],
        {cv, {0, 1}}], {wq, If[nm === "2p1d", {0, 1}, {1}]}], {dpf, {0, 4/3}}]],
    {sv, states}];
];
If[KeyExistsQ[results, "1p1d"], speedScan["1p1d"]];
Print["  1p1d scanned (", stamp[], ")"];
If[KeyExistsQ[results, "2p1d"], speedScan["2p1d"]];
If[speedRows =!= {},
  Export[FileNameJoin[{outDir, "two_way_char_speeds.tsv"}],
    Prepend[Map[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], #[[6]], #[[7]], ToString[#[[8]], InputForm],
        N[#[[9]], 8], N[#[[10]], 4], #[[11]], StringRiffle[ToString[NumberForm[N[#], 6]] & /@ #[[12]], " "], ToString[#[[13]], InputForm]} &, speedRows],
      {"chart", "T", "|u|", "DsT", "load", "c", "wv", "dir", "max|lambda|", "max|Im lambda|", "Jordan defect", "Re lambda (sorted)", "dpp/tauS"}],
    "TSV"];
  Do[Module[{sel = Select[speedRows, #[[1]] === nm && #[[6]] === cv &], hyp, nbad = 0},
      If[sel =!= {},
        hyp = Select[sel, #[[10]] < 10^-20 && #[[11]] == 0 &];
        info[nm, "  c=", cv, ":  ", Length[sel], " (state, direction) pairs;  max|Re lambda| = ",
          NumberForm[N[Max[sel[[All, 9]]]], 5], ";  max|Im lambda| = ", sci[N[Max[sel[[All, 10]]]]],
          ";  Jordan defects ", Total[sel[[All, 11]]],
          ";  hyperbolic in ", Length[hyp], "/", Length[sel], ";  subluminal in ", Count[sel, r_ /; r[[9]] <= 1], "/", Length[sel]];
        Do[If[!(r[[10]] < 10^-20 && r[[11]] == 0) && (nbad = nbad + 1) <= 3,
            info["      NOT hyperbolic: T=", r[[2]], " |u|=", r[[3]], " DsT=", r[[4]], " ", r[[5]], " wv=", r[[7]], " dir=", r[[8]],
              " dpp/tauS=", r[[13]], "  max|Im| = ", sci[N[r[[10]]]], "  defect ", r[[11]]]], {r, sel}]]],
    {nm, {"1p1d", "2p1d"}}, {cv, {0, 1}}];
  info["full table: ", FileNameJoin[{outDir, "two_way_char_speeds.tsv"}]]];
