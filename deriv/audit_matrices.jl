# audit_matrices.jl -- an independent audit of every matrix in src/Matrix, and of the drivers around them.
# Shares no code with the Wolfram scripts: each system is rebuilt here from its tensor equations, numerically,
# with ForwardDiff for the partial derivatives and explicit Christoffels. The residual R(q, dq) is affine in the
# first derivatives, so A_t, A_x, A_y, S are read off by unit vectors. Comparison is on the solved form
# (A_t^-1 A_x, A_t^-1 A_y, A_t^-1 S), which does not care about row order or row normalisation.
#
# Prints one line per check, "CHK <id> <PASS|FAIL> <worst> <what>"; exit 1 if a check fails. Run by
# derive_drivers.wls, or by itself:  julia --project=<env with Fluidum, StaticArrays, ForwardDiff> deriv/audit_matrices.jl
using ForwardDiff, LinearAlgebra, Random, Printf, StaticArrays
using Fluidum
const F = Fluidum

# ---------------------------------------------------------------- geometry
ginv_diag(chart, x) = chart === :cyl ? [-1one(x[1]), 1one(x[1]), 1/x[2]^2, 1/x[1]^2] :
                                       [-1one(x[1]), 1one(x[1]), 1one(x[1]), 1/x[1]^2]
g_diag(chart, x)    = chart === :cyl ? [-1one(x[1]), 1one(x[1]), x[2]^2, x[1]^2] :
                                       [-1one(x[1]), 1one(x[1]), 1one(x[1]), x[1]^2]
function christoffel(chart, x)
    dg = ForwardDiff.jacobian(y -> g_diag(chart, y), x)    # dg[a, b] = d_b g_aa
    g = g_diag(chart, x)
    Γ = zeros(4, 4, 4)                                     # Γ[a, b, c] = Γ^a_bc
    for a in 1:4, b in 1:4, c in 1:4
        s = 0.0
        a == c && (s += dg[a, b])
        a == b && (s += dg[a, c])
        b == c && (s -= dg[b, a])
        Γ[a, b, c] = s / (2 * g[a])
    end
    Γ
end
# ∇_d of an all-upper tensor valued function; index d is LAST
function cov(f, x0, chart)
    F0 = f(x0)
    r = ndims(F0)
    J = ForwardDiff.jacobian(y -> vec(f(y)), x0)
    D = reshape(J, ntuple(_ -> 4, r)..., 4)
    Γ = christoffel(chart, x0)
    out = Array{Float64}(D)
    for I in CartesianIndices(out)
        idx = Tuple(I); d = idx[end]
        for slot in 1:r, c in 1:4
            jdx = Base.setindex(idx[1:r], c, slot)
            out[I] += Γ[idx[slot], d, c] * F0[jdx...]
        end
    end
    out
end

# ---------------------------------------------------------------- tensors from the field layouts
uvec(chart, a, b = zero(a)) = [sqrt(1 + a^2 + b^2), a, b * one(a), zero(a)]   # (u^τ, u^1, u^2, 0)
function pi_cyl(ur, pφ, pη, x)
    uτ = sqrt(1 + ur^2); prr = -uτ^2 * (pη + pφ); ptr = ur * prr / uτ; ptt = ur * ptr / uτ
    M = [ptt ptr zero(ptt) zero(ptt); ptr prr zero(ptt) zero(ptt);
         zero(ptt) zero(ptt) pφ/x[2]^2 zero(ptt); zero(ptt) zero(ptt) zero(ptt) pη/x[1]^2]
    M
end
function pi_car(ux, uy, pyy, pe, pxy, x)
    uτ2 = 1 + ux^2 + uy^2; uτ = sqrt(uτ2)
    # trace: -π^ττ + π^xx + π^yy + π^η_η = 0, π^ττ = (ux²π^xx + 2 ux uy π^xy + uy²π^yy)/uτ²
    pxx = (-(pyy + pe) * uτ2 + 2ux * uy * pxy + uy^2 * pyy) / (1 + uy^2)
    ptx = (ux * pxx + uy * pxy) / uτ; pty = (ux * pxy + uy * pyy) / uτ; ptt = (ux * ptx + uy * pty) / uτ
    z = zero(ptt)
    [ptt ptx pty z; ptx pxx pxy z; pty pxy pyy z; z z z pe/x[1]^2]
end
nu_up(u, n1, n2 = zero(n1)) = [(u[2] * n1 + u[3] * n2) / u[1], n1, n2 * one(n1), zero(n1)]
Delta(chart, u, x) = Diagonal(ginv_diag(chart, x)) + u * u'

# projectors, evaluated at x0 only
function proj_tt(A, u, chart, x)             # Δ^{μν}_{αβ} A^{αβ}
    g = g_diag(chart, x); ul = g .* u
    Dmix = Matrix{Float64}(I, 4, 4) + u * ul'        # Δ^μ_α
    Dup = Delta(chart, u, x)
    Dlow = Diagonal(g) + ul * ul'
    S = Dmix * ((A + A') / 2) * Dmix'
    S - Dup * sum(Dlow .* A) / 3
end

# ---------------------------------------------------------------- medium
struct MedEOS{A}; pf::A; T0::Float64; end
pres(m, T) = m.pf(T)
dpres(m, T) = ForwardDiff.derivative(m.pf, T)
taylor_p(T0, p0, p1, p2) = T -> p0 + p1 * (T - T0) + p2 * (T - T0)^2 / 2
Tm_up(m, T, u, pi, Pi, chart, x) = begin
    p = pres(m, T); e = T * dpres(m, T) - p
    (e + p + Pi) .* (u * u') .+ (p + Pi) .* Diagonal(ginv_diag(chart, x)) .+ pi
end

function shear_row(P, x0, chart, η, τπ, dpp)
    s = P(x0); u = s.u
    du = cov(y -> P(y).u, x0, chart)                 # du[a, c] = ∇_c u^a
    θ = tr(du)
    gi = ginv_diag(chart, x0)
    gradu = [gi[a] * du[b, a] for a in 1:4, b in 1:4]   # ∇^a u^b
    σ = proj_tt(gradu, u, chart, x0)
    dπ = cov(y -> P(y).pi, x0, chart)
    Dπ = [sum(u[c] * dπ[a, b, c] for c in 1:4) for a in 1:4, b in 1:4]
    τπ .* proj_tt(Dπ, u, chart, x0) .+ s.pi .+ 2η .* σ .+ dpp .* s.pi .* θ
end
function bulk_row(P, x0, chart, ζ, τΠ)
    s = P(x0)
    du = cov(y -> P(y).u, x0, chart); θ = tr(du)
    dΠ = ForwardDiff.gradient(y -> P(y).Pi, x0)
    τΠ * dot(s.u, dΠ) + s.Pi + ζ * θ
end

# ---------------------------------------------------------------- charm
struct CharmEOS{A,B,C,D,E}; nf::A; hf::B; If::C; aQ::Float64; bf::D; cf::E; end
TQ_up(c, T, α, u, ν, chart, x) = begin
    n = c.nf(T, α); h = c.hf(T); P0 = n * T
    (n * h - P0) .* (u * u') .+ P0 .* Delta(chart, u, x) .+ h .* (u * ν' .+ ν * u')
end
function S3_up(c, T, α, u, ν, dM, chart, x)
    Dl = Delta(chart, u, x)
    I31 = c.If(T, α); A = c.aQ; B = c.bf(T); C = c.cf(T)
    [I31 * (u[l] * Dl[a, b] + u[a] * Dl[l, b] + u[b] * Dl[l, a]) +
     A * (u[l] * u[a] * ν[b] + u[l] * u[b] * ν[a] + u[a] * u[b] * ν[l]) +
     B * (Dl[l, a] * ν[b] + Dl[l, b] * ν[a] + Dl[a, b] * ν[l]) +
     C * (u[l] * dM[a, b] + u[a] * dM[l, b] + u[b] * dM[l, a]) for l in 1:4, a in 1:4, b in 1:4]
end
function div_N(P, x0, chart, c)
    d = cov(y -> (s = P(y); c.nf(s.T, s.α) .* s.u .+ s.ν), x0, chart)
    tr(d)
end
function nu_row_consistent(P, x0, chart, c, Ds; cm = false)    # (D_s/T) Δ^μ_λ ∇_ν (T_Q + c δM)^{νλ} + ν^μ
    s = P(x0)
    dTQ = cov(y -> (q = P(y); TQ_up(c, q.T, q.α, q.u, q.ν, chart, y) .+ (cm ? q.dM : zero(q.u * q.u'))), x0, chart)
    div = [sum(dTQ[n, l, n] for n in 1:4) for l in 1:4]
    g = g_diag(chart, x0); ul = g .* s.u
    Dmix = Matrix{Float64}(I, 4, 4) + s.u * ul'
    (Ds / s.T) .* (Dmix * div) .+ s.ν
end
function nu_row_fick(P, x0, chart, τn, κ, dαf = y -> P(y).α)
    s = P(x0); u = s.u
    dν = cov(y -> P(y).ν, x0, chart)
    Dν = [sum(u[c] * dν[a, c] for c in 1:4) for a in 1:4]
    g = g_diag(chart, x0); ul = g .* u
    Dmix = Matrix{Float64}(I, 4, 4) + u * ul'
    dα = ForwardDiff.gradient(dαf, x0)
    τn .* (Dmix * Dν) .+ s.ν .+ κ .* (Delta(chart, u, x0) * dα)
end
function m2_rows(P, x0, chart, c, Ds)              # returns (traceless E^{μν}, trace row)
    s = P(x0); u = s.u
    dS = cov(y -> (q = P(y); S3_up(c, q.T, q.α, q.u, q.ν, q.dM, chart, y)), x0, chart)
    divS = [sum(dS[l, a, b, l] for l in 1:4) for a in 1:4, b in 1:4]
    g = g_diag(chart, x0); ul = g .* u
    Dmix = Matrix{Float64}(I, 4, 4) + u * ul'
    Dlow = Diagonal(g) + ul * ul'
    Dup = Delta(chart, u, x0)
    E = Dmix * divS * Dmix'
    trE = sum(Dlow .* E)        # 𝓔^λ_λ
    pre = Ds / (2 * s.T)
    TL = pre .* (E .- Dup .* trE ./ 3) .+ s.piQ
    TR = pre * trE / 3 + s.PiQ
    TL, TR
end

# ---------------------------------------------------------------- reading matrices off a residual
function extract(R, nf, ndim)    # R(dq) with dq = ntuple(ndim) of derivative vectors
    zv = zeros(nf)
    S = R(ntuple(_ -> zv, ndim))
    A = ntuple(a -> begin
        M = zeros(nf, nf)
        for j in 1:nf
            e = zeros(nf); e[j] = 1
            M[:, j] = R(ntuple(b -> b == a ? e : zv, ndim)) .- S
        end
        M
    end, ndim)
    (A..., S)
end
line(q0, dq, x, x0) = begin    # fields linear in x around x0
    v = q0 .+ zero(x[1])
    for a in eachindex(dq); v = v .+ dq[a] .* (x[a] - x0[a]); end
    v
end
function compare(mine, kern; label = "")
    nd = length(mine) - 2
    worst = 0.0
    Aim = mine[1]; Aik = Matrix(kern[1])
    for a in 2:nd+1
        X = Aim \ mine[a]; Y = Aik \ Matrix(kern[a])
        sc = max(maximum(abs, X), 1e-300)
        worst = max(worst, maximum(abs, X .- Y) / sc)
    end
    X = Aim \ mine[end]; Y = Aik \ Vector(kern[end])
    worst = max(worst, maximum(abs, X .- Y) / max(maximum(abs, X), 1e-300))
    worst
end

# ---------------------------------------------------------------- random inputs
rs(a, b) = a + (b - a) * rand()
const OVERRIDE = Ref{Any}((;))      # driver checks: physical inputs replace the random ones
const LAST = Ref{Any}(nothing)      # the state a system was built on
randthermo() = merge(randthermo0(), OVERRIDE[])
function randthermo0()
    (T = rs(.6, 1.6), α = rs(-1, 1), p0 = rs(.5, 2), p1 = rs(1, 3), p2 = rs(1, 3), η = rs(.1, 1), τπ = rs(.5, 2),
     ζ = rs(.05, .5), τΠ = rs(.5, 2), dpp = rs(.2, 2), n = rs(.5, 2), nT = rs(.2, 2), na = rs(.5, 2), h = rs(1, 4), hT = rs(.2, 1),
     I = rs(.5, 3), IT = rs(.2, 2), Ia = rs(.2, 2), aQ = rs(.5, 3), b = rs(.5, 3), bT = rs(.2, 2), c = rs(.5, 3), cT = rs(.2, 2),
     Ds = rs(.2, 1.5), τn = rs(.5, 2), κ = rs(.2, 1.5), a = rs(.5, 3))
end
charm_of(p) = CharmEOS((T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α), T -> p.h + p.hT * (T - p.T),
                       (T, α) -> p.I + p.IT * (T - p.T) + p.Ia * (α - p.α), p.aQ,
                       T -> p.b + p.bT * (T - p.T), T -> p.c + p.cT * (T - p.T))
med_of(p; conformal = false) = conformal ? MedEOS(T -> p.a * T^4, p.T) : MedEOS(taylor_p(p.T, p.p0, p.p1, p.p2), p.T)

# ---------------------------------------------------------------- the systems
# Each returns (mine, kernel) on a random state; `bad` multiplies one coefficient in the KERNEL's arguments
# (the control: bad = 1.2 must fail).

function sys_visc1d(; bad = 1.0)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]
    q0 = [p.T, rs(-1, 1), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3)]; m = med_of(p)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:cyl, w[2]); (T = w[1], u = u, pi = pi_cyl(w[2], w[3], w[4], y), Pi = w[5]))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :cyl, y)), x0, :cyl)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :cyl, p.η, p.τπ, p.dpp)
        [dv[1], dv[2], E[3, 3], E[4, 4], bulk_row(P, x0, :cyl, p.ζ, p.τΠ)]
    end
    LAST[] = Base.@locals()
    extract(R, 5, 2), F.one_d_viscous_matrix(SVector{5}(q0), τ, r, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ, p.dpp)
end
function sys_visc1d_MIS_block(; bad = 1.0)   # the medium of the 7-field :MIS kernel, dpp = 4/3 τπ
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]
    q0 = [p.T, rs(-1, 1), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), p.α, rs(-.3, .3)]; m = med_of(p)
    nf = (T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:cyl, w[2]);
                  (T = w[1], u = u, pi = pi_cyl(w[2], w[3], w[4], y), Pi = w[5], α = w[6], ν = nu_up(u, w[7])))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :cyl, y)), x0, :cyl)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :cyl, p.η, p.τπ, 4/3 * p.τπ)
        dN = tr(cov(y -> (s = P(y); nf(s.T, s.α) .* s.u .+ s.ν), x0, :cyl))
        nr = nu_row_fick(P, x0, :cyl, p.τn, p.κ)
        [dv[1], dv[2], E[3, 3], E[4, 4], bulk_row(P, x0, :cyl, p.ζ, p.τΠ), dN, nr[2]]
    end
    LAST[] = Base.@locals()
    extract(R, 7, 2), F.one_d_viscous_matrix_fugacity(SVector{7}(q0), τ, r, p.p1, p.p1, p.p2, 0, 0, 0, p.ζ, p.η, p.τπ, p.τΠ, p.n, p.nT, p.na, p.τn, p.κ * bad)
end
function sys_visc1d_gamma(; bad = 1.0)       # φ6 = γ, n = γ n0(T), κ multiplying ∂γ
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]
    q0 = [p.T, rs(-1, 1), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(.3, 2), rs(-.3, .3)]; m = med_of(p)
    n0f = T -> p.n + p.nT * (T - p.T)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:cyl, w[2]);
                  (T = w[1], u = u, pi = pi_cyl(w[2], w[3], w[4], y), Pi = w[5], α = w[6], ν = nu_up(u, w[7])))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :cyl, y)), x0, :cyl)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :cyl, p.η, p.τπ, 4/3 * p.τπ)
        dN = tr(cov(y -> (s = P(y); s.α * n0f(s.T) .* s.u .+ s.ν), x0, :cyl))
        nr = nu_row_fick(P, x0, :cyl, p.τn, p.κ)
        [dv[1], dv[2], E[3, 3], E[4, 4], bulk_row(P, x0, :cyl, p.ζ, p.τΠ), dN, nr[2]]
    end
    LAST[] = Base.@locals()
    extract(R, 7, 2), F.one_d_viscous_matrix_gamma(SVector{7}(q0), τ, r, p.p1, p.p1, p.p2, 0, 0, 0, p.ζ, p.η, p.τπ, p.τΠ, p.n, p.nT, p.na, p.τn, p.κ * bad)
end
function ideal_base(p, q0, τ, r, x0, nofq)
    m = med_of(p)
    dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:cyl, w[2]); (T = w[1], u = u, q3 = w[3]))
        z = zeros(4, 4)
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, z, 0.0, :cyl, y)), x0, :cyl)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        dN = tr(cov(y -> (s = P(y); nofq(s.T, s.q3) .* s.u), x0, :cyl))
        [dv[1], dv[2], dN]
    end
end
function sys_ideal_density(; bad = 1.0)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; q0 = [p.T, rs(-1, 1), rs(.3, 2)]
    LAST[] = Base.@locals()
    extract(ideal_base(p, q0, τ, r, x0, (T, n) -> n), 3, 2), F.one_d_ideal_matrix_ruwen4(SVector{3}(q0), τ, r, p.p1 * bad, p.p2)
end
function sys_ideal_fugacity(; bad = 1.0)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; q0 = [p.T, rs(-1, 1), p.α]
    nf = (T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α)
    LAST[] = Base.@locals()
    extract(ideal_base(p, q0, τ, r, x0, nf), 3, 2), F.one_d_ideal_matrix_ruwen5(SVector{3}(q0), τ, r, p.p1 * bad, p.p2, p.n, p.nT, p.na)
end
# the back-reaction files: P_tot = P(T) + P_h(T, α). enth = :T  e + P = T ∂_T P_tot|α   ; :Tμ  e + P = T ∂_T P_tot + α ∂_α P_tot
function pbr_model(p, q0, τ, r, x0, enth, nmode)
    Ph1, Ph2, Hh11, Hh12, Hh22 = p.I, p.IT, p.aQ, p.b, p.c   # dPh/dT, dPh/dα, hessian
    Ph(T, α) = p.h + Ph1 * (T - p.T) + Ph2 * (α - p.α) + (Hh11 * (T - p.T)^2 + 2Hh12 * (T - p.T) * (α - p.α) + Hh22 * (α - p.α)^2) / 2
    dT(T, α) = Ph1 + Hh11 * (T - p.T) + Hh12 * (α - p.α)
    da(T, α) = Ph2 + Hh12 * (T - p.T) + Hh22 * (α - p.α)
    m = med_of(p)
    dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (T = w[1], u = uvec(:cyl, w[2]), q3 = w[3]))
        # q = n = ∂_α P_h / T, the charge density at μ = α T
        αof(T, q) = nmode === :density ? p.α + (q * T - Ph2 - Hh12 * (T - p.T)) / Hh22 : q
        d = cov(y -> (s = P(y); α = αof(s.T, s.q3);
                 Ptot = pres(m, s.T) + Ph(s.T, α);
                 W = s.T * (dpres(m, s.T) + dT(s.T, α)) + (enth === :Tμ ? α * da(s.T, α) : zero(α));
                 W .* (s.u * s.u') .+ Ptot .* Diagonal(ginv_diag(:cyl, y))), x0, :cyl)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        dN = tr(cov(y -> (s = P(y); nn = nmode === :density ? s.q3 : p.n + p.nT * (s.T - p.T) + p.na * (s.q3 - p.α); nn .* s.u), x0, :cyl))
        [dv[1], dv[2], dN]
    end
end
function sys_pbr_density(; enth = :T, swap = false)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; q0 = [p.T, rs(-1, 1), p.IT / p.T]
    At, Ax, src = F.one_d_ideal_matrix_ruwen6(SVector{3}(q0), τ, r, p.p1, p.p2, p.I, p.IT, p.aQ, p.b, p.c)
    At, Ax = Matrix(At), Matrix(Ax)
    swap && ((At[1:2, 3], Ax[1:2, 3]) = (Ax[1:2, 3], At[1:2, 3]))
    LAST[] = Base.@locals()
    extract(pbr_model(p, q0, τ, r, x0, enth, :density), 3, 2), (At, Ax, src)
end
function sys_pbr_fugacity(; enth = :T)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; q0 = [p.T, rs(-1, 1), p.α]
    LAST[] = Base.@locals()
    extract(pbr_model(p, q0, τ, r, x0, enth, :fug), 3, 2),
        F.one_d_ideal_matrix_ruwen7(SVector{3}(q0), τ, r, p.p1, p.p2, p.I, p.IT, p.aQ, p.b, p.c, p.n, p.na, p.nT)
end

# background systems ---------------------------------------------------------------
function bg1d()
    (T = nothing,)
end
function bgP1d(p, x0, ur, g)   # returns (T(y), ur(y))
    τ, r = x0[1], x0[2]
    y -> (p.T + g.dtT * (y[1] - τ) + g.drT * (y[2] - r), ur + g.dtur * (y[1] - τ) + g.drur * (y[2] - r))
end
randg1d() = (dtT = rs(-1, 1), drT = rs(-1, 1), dtur = rs(-1, 1), drur = rs(-1, 1))
function sys_BG_fick1d(; bad = 1.0)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; ur = rs(-1, 1); g = randg1d(); bgf = bgP1d(p, x0, ur, g)
    q0 = [p.α, rs(-.3, .3)]
    nf = (T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (Tb, ub) = bgf(y); u = uvec(:cyl, ub); (T = Tb, u = u, α = w[1], ν = nu_up(u, w[2])))
        dN = tr(cov(y -> (s = P(y); nf(s.T, s.α) .* s.u .+ s.ν), x0, :cyl))
        nr = nu_row_fick(P, x0, :cyl, p.τn, p.κ)
        [dN, nr[2]]
    end
    LAST[] = Base.@locals()
    extract(R, 2, 2), F.one_d_viscous_matrix_fugacity_BG_only(SVector{2}(q0), τ, r, ur, g.dtT, g.drT, g.drur, g.dtur, p.n, p.na, p.nT, p.τn, p.κ * bad)
end
function sys_BG_cons1d(; bad = 1.0)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; ur = rs(-1, 1); g = randg1d(); bgf = bgP1d(p, x0, ur, g)
    q0 = [p.α, rs(-.3, .3)]; c = charm_of(p)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (Tb, ub) = bgf(y); u = uvec(:cyl, ub); (T = Tb, u = u, α = w[1], ν = nu_up(u, w[2])))
        [div_N(P, x0, :cyl, c), nu_row_consistent(P, x0, :cyl, c, p.Ds)[2]]
    end
    LAST[] = Base.@locals()
    extract(R, 2, 2), F.one_d_HQ_BG_consistent_matrix(SVector{2}(q0), τ, r, p.T, ur, g.dtT, g.drT, g.dtur, g.drur, p.n, p.nT, p.na, p.h, p.hT, p.Ds * bad)
end
function sys_BG_m21d(; bad = 1.0)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; ur = rs(-1, 1); g = randg1d(); bgf = bgP1d(p, x0, ur, g)
    q0 = [p.α, rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3)]; c = charm_of(p)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (Tb, ub) = bgf(y); u = uvec(:cyl, ub); pQ = pi_cyl(ub, w[3], w[4], y);
                  (T = Tb, u = u, α = w[1], ν = nu_up(u, w[2]), piQ = pQ, PiQ = w[5], dM = pQ .+ w[5] .* Delta(:cyl, u, y)))
        TL, TR = m2_rows(P, x0, :cyl, c, p.Ds)
        [div_N(P, x0, :cyl, c), nu_row_consistent(P, x0, :cyl, c, p.Ds)[2], TL[3, 3], TL[4, 4], TR]
    end
    LAST[] = Base.@locals()
    extract(R, 5, 2), F.one_d_HQ_BG_m2_matrix(SVector{5}(q0), τ, r, p.T, ur, g.dtT, g.drT, g.dtur, g.drur, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds * bad)
end

randg2d() = (dtT = rs(-1, 1), dxT = rs(-1, 1), dyT = rs(-1, 1), dtux = rs(-1, 1), dxux = rs(-1, 1), dyux = rs(-1, 1), dtuy = rs(-1, 1), dxuy = rs(-1, 1), dyuy = rs(-1, 1))
bgP2d(p, x0, ux, uy, g) = y -> begin
    dt, dx, dy = y[1] - x0[1], y[2] - x0[2], y[3] - x0[3]
    (p.T + g.dtT * dt + g.dxT * dx + g.dyT * dy, ux + g.dtux * dt + g.dxux * dx + g.dyux * dy, uy + g.dtuy * dt + g.dxuy * dx + g.dyuy * dy)
end
function sys_BG_fick2d(; bad = 1.0)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]; ux, uy = rs(-1, 1), rs(-1, 1); g = randg2d(); bgf = bgP2d(p, x0, ux, uy, g)
    q0 = [p.α, rs(-.3, .3), rs(-.3, .3)]
    nf = (T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (Tb, a, b) = bgf(y); u = uvec(:car, a, b); (T = Tb, u = u, α = w[1], ν = nu_up(u, w[2], w[3])))
        dN = tr(cov(y -> (s = P(y); nf(s.T, s.α) .* s.u .+ s.ν), x0, :car))
        nr = nu_row_fick(P, x0, :car, p.τn, p.κ)
        [dN, nr[2], nr[3]]
    end
    LAST[] = Base.@locals()
    extract(R, 3, 3), F.two_d_HQ_BG_matrix(q0[1], q0[2], q0[3], τ, p.T, ux, uy, g.dtT, g.dxT, g.dyT, g.dtux, g.dxux, g.dyux, g.dtuy, g.dxuy, g.dyuy, p.n, p.nT, p.na, p.τn, p.κ * bad)
end
function sys_live_rows(; bad = 1.0)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]
    q0 = [p.T, rs(-1, 1), rs(-1, 1), p.α, rs(-.3, .3), rs(-.3, .3)]
    nf = (T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:car, w[2], w[3]); (T = w[1], u = u, α = w[4], ν = nu_up(u, w[5], w[6])))
        dN = tr(cov(y -> (s = P(y); nf(s.T, s.α) .* s.u .+ s.ν), x0, :car))
        nr = nu_row_fick(P, x0, :car, p.τn, p.κ)
        [dN, nr[2], nr[3]]
    end
    # 3x6: extract by hand, compare after multiplying by the inverse of the charm block
    zv = zeros(6); S = R((zv, zv, zv))
    A = [begin M = zeros(3, 6); for j in 1:6; e = zeros(6); e[j] = 1; M[:, j] = R(ntuple(b -> b == a ? e : zv, 3)) .- S; end; M end for a in 1:3]
    AtL, AxL, AyL, srcL = F.two_d_HQ_live_rows(q0[4], q0[5], q0[6], τ, q0[1], q0[2], q0[3], p.n, p.nT, p.na, p.τn, p.κ * bad)
    Bm = A[1][:, 4:6]; Bk = Matrix(AtL)[:, 4:6]
    worst = 0.0
    for (X, Y) in ((A[1], AtL), (A[2], AxL), (A[3], AyL))
        a = Bm \ X; b = Bk \ Matrix(Y); worst = max(worst, maximum(abs, a .- b) / maximum(abs, a))
    end
    a = Bm \ S; b = Bk \ Vector(srcL); worst = max(worst, maximum(abs, a .- b) / maximum(abs, a))
    worst
end
function sys_BG_cons2d(; bad = 1.0)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]; ux, uy = rs(-1, 1), rs(-1, 1); g = randg2d(); bgf = bgP2d(p, x0, ux, uy, g)
    q0 = [p.α, rs(-.3, .3), rs(-.3, .3)]; c = charm_of(p)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (Tb, a, b) = bgf(y); u = uvec(:car, a, b); (T = Tb, u = u, α = w[1], ν = nu_up(u, w[2], w[3])))
        nr = nu_row_consistent(P, x0, :car, c, p.Ds)
        [div_N(P, x0, :car, c), nr[2], nr[3]]
    end
    LAST[] = Base.@locals()
    extract(R, 3, 3), F.two_d_HQ_BG_consistent_matrix(SVector{3}(q0), τ, p.T, ux, uy, g.dtT, g.dxT, g.dyT, g.dtux, g.dxux, g.dyux, g.dtuy, g.dxuy, g.dyuy, p.n, p.nT, p.na, p.h, p.hT, p.Ds * bad)
end
function sys_BG_m22d(; bad = 1.0)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]; ux, uy = rs(-1, 1), rs(-1, 1); g = randg2d(); bgf = bgP2d(p, x0, ux, uy, g)
    q0 = [p.α, rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3)]; c = charm_of(p)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); (Tb, a, b) = bgf(y); u = uvec(:car, a, b); pQ = pi_car(a, b, w[4], w[5], w[6], y);
                  (T = Tb, u = u, α = w[1], ν = nu_up(u, w[2], w[3]), piQ = pQ, PiQ = w[7], dM = pQ .+ w[7] .* Delta(:car, u, y)))
        nr = nu_row_consistent(P, x0, :car, c, p.Ds)
        TL, TR = m2_rows(P, x0, :car, c, p.Ds)
        [div_N(P, x0, :car, c), nr[2], nr[3], TL[3, 3], TL[4, 4], TL[2, 3], TR]
    end
    LAST[] = Base.@locals()
    extract(R, 7, 3), F.two_d_HQ_BG_m2_matrix(SVector{7}(q0), τ, p.T, ux, uy, g.dtT, g.dxT, g.dyT, g.dtux, g.dxux, g.dyux, g.dtuy, g.dxuy, g.dyuy, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds * bad)
end

# 2+1D medium -----------------------------------------------------------------------
function sys_visc2d(kern; dppv = nothing, bad = 1.0)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]
    q0 = [p.T, rs(-1, 1), rs(-1, 1), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3)]; m = med_of(p)
    dpp = dppv === nothing ? p.dpp : dppv
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:car, w[2], w[3]); (T = w[1], u = u, pi = pi_car(w[2], w[3], w[4], w[5], w[6], y), Pi = w[7]))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :car, y)), x0, :car)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :car, p.η, p.τπ, dpp)
        [dv[1], dv[2], dv[3], E[3, 3], E[4, 4], E[2, 3], bulk_row(P, x0, :car, p.ζ, p.τΠ)]
    end
    LAST[] = Base.@locals()
    k = kern === :derived ? F.two_d_viscous_matrix_derived(SVector{7}(q0), τ, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ, dpp) :
                            F.two_d_viscous_matrix(SVector{7}(q0), τ, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ)
    extract(R, 7, 3), k
end
function sys_visc_HQ_upstream(; dpp = 0.0, medium = :legacy)   # two_d_viscous_HQ_matrix, 10 fields (upstream, not derived)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]
    q0 = [p.T, rs(-1, 1), rs(-1, 1), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), p.α, rs(-.3, .3), rs(-.3, .3)]; m = med_of(p)
    nf = (T, α) -> p.n + p.nT * (T - p.T) + p.na * (α - p.α)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:car, w[2], w[3]);
                  (T = w[1], u = u, pi = pi_car(w[2], w[3], w[4], w[5], w[6], y), Pi = w[7], α = w[8], ν = nu_up(u, w[9], w[10])))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :car, y)), x0, :car)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :car, p.η, p.τπ, dpp === :phys ? p.dpp : dpp)
        dN = tr(cov(y -> (s = P(y); nf(s.T, s.α) .* s.u .+ s.ν), x0, :car))
        nr = nu_row_fick(P, x0, :car, p.τn, p.κ)
        [dv[1], dv[2], dv[3], E[3, 3], E[4, 4], E[2, 3], bulk_row(P, x0, :car, p.ζ, p.τΠ), dN, nr[2], nr[3]]
    end
    LAST[] = Base.@locals()
    k = F.two_d_viscous_HQ_matrix(SVector{10}(q0), τ, p.p1, p.p1, p.p2, p.ζ, p.η, p.τπ, p.τΠ, p.n, p.nT, p.na, p.τn, p.κ / p.n)   # its argument is κ/n: the rows carry D_s n
    if medium === :derived    # the kernel's charm rows, with the derived medium block in place of its legacy one
        d = F.two_d_viscous_matrix_derived(SVector{7}(q0[1:7]), τ, p.p0, p.p1, p.p2, p.ζ, p.η, p.τπ, p.τΠ, dpp === :phys ? p.dpp : dpp)
        k = ntuple(i -> begin
            M = copy(i <= 3 ? Matrix(k[i]) : Vector(k[i]))
            i <= 3 ? (M[1:7, 1:7] = Matrix(d[i])) : (M[1:7] = Vector(d[i]))
            M end, 4)
    end
    extract(R, 10, 3), k
end

# two-way ---------------------------------------------------------------------------
function sys_tw1d(; conformal = false, bad = 1.0, cm = false)
    p = randthermo(); τ, r = rs(.5, 2), rs(.5, 2); x0 = [τ, r, 1.0, 0.0]; m = med_of(p; conformal); c = charm_of(p)
    if conformal
        q0 = [p.T, rs(-1, 1), rs(-.3, .3), rs(-.3, .3), p.α, rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3)]
        ix = (T = 1, u = 2, pf = 3, pe = 4, Π = 0, α = 5, ν = 6, Qf = 7, Qe = 8, QΠ = 9)
    else
        q0 = [p.T, rs(-1, 1), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), p.α, rs(-.3, .3), rs(-.3, .3), rs(-.3, .3), rs(-.3, .3)]
        ix = (T = 1, u = 2, pf = 3, pe = 4, Π = 5, α = 6, ν = 7, Qf = 8, Qe = 9, QΠ = 10)
    end
    nf = length(q0)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:cyl, w[ix.u]); pQ = pi_cyl(w[ix.u], w[ix.Qf], w[ix.Qe], y);
                  (T = w[ix.T], u = u, pi = pi_cyl(w[ix.u], w[ix.pf], w[ix.pe], y), Pi = ix.Π == 0 ? zero(w[1]) : w[ix.Π],
                   α = w[ix.α], ν = nu_up(u, w[ix.ν]), piQ = pQ, PiQ = w[ix.QΠ], dM = pQ .+ w[ix.QΠ] .* Delta(:cyl, u, y)))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :cyl, y) .+ TQ_up(c, s.T, s.α, s.u, s.ν, :cyl, y) .+ (cm ? s.dM : zero(s.u * s.u'))), x0, :cyl)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :cyl, p.η, p.τπ, p.dpp)
        TL, TR = m2_rows(P, x0, :cyl, c, p.Ds)
        rows = Any[dv[1], dv[2], E[3, 3], E[4, 4]]
        conformal || push!(rows, bulk_row(P, x0, :cyl, p.ζ, p.τΠ))
        append!(rows, [div_N(P, x0, :cyl, c), nu_row_consistent(P, x0, :cyl, c, p.Ds; cm)[2], TL[3, 3], TL[4, 4], TR])
        Float64.(rows)
    end
    LAST[] = Base.@locals()
    k = cm ? (conformal ?
        F.one_d_two_way_cM_matrix_conformal(SVector{9}(q0), τ, r, p.a, p.η * bad, p.τπ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds) :
        F.one_d_two_way_cM_matrix(SVector{10}(q0), τ, r, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds)) :
        conformal ?
        F.one_d_two_way_matrix_conformal(SVector{9}(q0), τ, r, p.a, p.η * bad, p.τπ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds) :
        F.one_d_two_way_matrix(SVector{10}(q0), τ, r, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds)
    extract(R, nf, 2), k
end
function sys_tw2d(; conformal = false, bad = 1.0, cm = false)
    p = randthermo(); τ = rs(.5, 2); x0 = [τ, rs(-1, 1), rs(-1, 1), 0.0]; m = med_of(p; conformal); c = charm_of(p)
    rp() = rs(-.3, .3)
    if conformal
        q0 = [p.T, rs(-1, 1), rs(-1, 1), rp(), rp(), rp(), p.α, rp(), rp(), rp(), rp(), rp(), rp()]
        ix = (yy = 4, e = 5, xy = 6, Π = 0, α = 7, n1 = 8, n2 = 9, Qyy = 10, Qe = 11, Qxy = 12, QΠ = 13)
    else
        q0 = [p.T, rs(-1, 1), rs(-1, 1), rp(), rp(), rp(), rp(), p.α, rp(), rp(), rp(), rp(), rp(), rp()]
        ix = (yy = 4, e = 5, xy = 6, Π = 7, α = 8, n1 = 9, n2 = 10, Qyy = 11, Qe = 12, Qxy = 13, QΠ = 14)
    end
    nf = length(q0)
    R = dq -> begin
        P = y -> (w = line(q0, dq, y, x0); u = uvec(:car, w[2], w[3]); pQ = pi_car(w[2], w[3], w[ix.Qyy], w[ix.Qe], w[ix.Qxy], y);
                  (T = w[1], u = u, pi = pi_car(w[2], w[3], w[ix.yy], w[ix.e], w[ix.xy], y), Pi = ix.Π == 0 ? zero(w[1]) : w[ix.Π],
                   α = w[ix.α], ν = nu_up(u, w[ix.n1], w[ix.n2]), piQ = pQ, PiQ = w[ix.QΠ], dM = pQ .+ w[ix.QΠ] .* Delta(:car, u, y)))
        d = cov(y -> (s = P(y); Tm_up(m, s.T, s.u, s.pi, s.Pi, :car, y) .+ TQ_up(c, s.T, s.α, s.u, s.ν, :car, y) .+ (cm ? s.dM : zero(s.u * s.u'))), x0, :car)
        dv = [sum(d[n, l, n] for n in 1:4) for l in 1:4]
        E = shear_row(P, x0, :car, p.η, p.τπ, p.dpp)
        TL, TR = m2_rows(P, x0, :car, c, p.Ds)
        nr = nu_row_consistent(P, x0, :car, c, p.Ds; cm)
        rows = Any[dv[1], dv[2], dv[3], E[3, 3], E[4, 4], E[2, 3]]
        conformal || push!(rows, bulk_row(P, x0, :car, p.ζ, p.τΠ))
        append!(rows, [div_N(P, x0, :car, c), nr[2], nr[3], TL[3, 3], TL[4, 4], TL[2, 3], TR])
        Float64.(rows)
    end
    LAST[] = Base.@locals()
    k = cm ? (conformal ?
        F.two_d_two_way_cM_matrix_conformal(SVector{13}(q0), τ, p.a, p.η * bad, p.τπ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds) :
        F.two_d_two_way_cM_matrix(SVector{14}(q0), τ, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds)) :
        conformal ?
        F.two_d_two_way_matrix_conformal(SVector{13}(q0), τ, p.a, p.η * bad, p.τπ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds) :
        F.two_d_two_way_matrix(SVector{14}(q0), τ, p.p0, p.p1, p.p2, p.ζ, p.η * bad, p.τπ, p.τΠ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds)
    extract(R, nf, 3), k
end

# ---------------------------------------------------------------- run
Random.seed!(20261030)
const N = 12
scan(f; kw...) = maximum(begin
    r = f(; kw...)
    r isa Tuple ? compare(r[1], r[2]) : r
end for _ in 1:N)
const NFAIL = Ref(0); const NCHK = Ref(0)
function chk(id, ok, w, what)
    NCHK[] += 1; ok || (NFAIL[] += 1)
    @printf("CHK %s %s %.2e %s\n", id, ok ? "PASS" : "FAIL", w, what)
end
# agrees to rounding, and the same comparison rejects a coefficient off by 20 %
function report(id, name, f; control = true, kw...)
    w = scan(f; kw...)
    ctl = control ? scan(f; bad = 1.2, kw...) : Inf
    chk(id, w < 1e-8 && ctl > 1e-4, w, name)
end
# known to disagree: the disagreement is the statement
report_differs(id, name, f; kw...) = (w = scan(f; kw...); chk(id, w > 1e-3, w, name))

report("A-visc2d", "two_d_viscous_matrix_derived", (; kw...) -> sys_visc2d(:derived; kw...))
report("A-BGf2d", "two_d_HQ_BG_matrix", sys_BG_fick2d)
report("A-live", "two_d_HQ_live_rows", sys_live_rows)
report("A-BGc1d", "one_d_HQ_BG_consistent_matrix", sys_BG_cons1d)
report("A-BGm1d", "one_d_HQ_BG_m2_matrix", sys_BG_m21d)
report("A-BGc2d", "two_d_HQ_BG_consistent_matrix", sys_BG_cons2d)
report("A-BGm2d", "two_d_HQ_BG_m2_matrix", sys_BG_m22d)
report("A-tw1d", "one_d_two_way_matrix", sys_tw1d)
report("A-tw1dc", "one_d_two_way_matrix_conformal", sys_tw1d; conformal = true)
report("A-tw2d", "two_d_two_way_matrix", sys_tw2d)
report("A-tw2dc", "two_d_two_way_matrix_conformal", sys_tw2d; conformal = true)
report("A-cm1d", "one_d_two_way_cM_matrix", sys_tw1d; cm = true)
report("A-cm1dc", "one_d_two_way_cM_matrix_conformal", sys_tw1d; conformal = true, cm = true)
report("A-cm2d", "two_d_two_way_cM_matrix", sys_tw2d; cm = true)
report("A-cm2dc", "two_d_two_way_cM_matrix_conformal", sys_tw2d; conformal = true, cm = true)
report_differs("A-cm0", "control: the S5 kernel differs from the system with dM in T_Q", (; kw...) -> begin
    r = sys_tw1d(; cm = true); r1 = LAST[]; p = r1[:p]; q0 = r1[:q0]
    (r[1], F.one_d_two_way_matrix(SVector{10}(q0), r1[:τ], r1[:r], p.p0, p.p1, p.p2, p.ζ, p.η, p.τπ, p.τΠ, p.dpp, p.n, p.nT, p.na, p.h, p.hT, p.I, p.IT, p.Ia, p.aQ, p.b, p.bT, p.c, p.cT, p.Ds)) end)
report("A-visc1d", "one_d_viscous_matrix", sys_visc1d)
report("A-MIS", "one_d_viscous_matrix_fugacity", sys_visc1d_MIS_block)
report("A-gamma", "one_d_viscous_matrix_gamma", sys_visc1d_gamma)
report("A-BGf1d", "one_d_viscous_matrix_fugacity_BG_only", sys_BG_fick1d)
report("A-r4", "one_d_ideal_matrix_ruwen4", sys_ideal_density)
report("A-r5", "one_d_ideal_matrix_ruwen5", sys_ideal_fugacity)
report_differs("A-legacy", "two_d_viscous_matrix differs from the equations (its ten wrong entries)", (; kw...) -> sys_visc2d(:legacy; dppv = 0.0, kw...))
report("A-HQ", "two_d_viscous_HQ_matrix: derived medium at dpp = 0, Fick charm rows with kappa = Ds n", sys_visc_HQ_upstream; control = false)
report("A-r7", "one_d_ideal_matrix_ruwen7, e + P = T d_T P_tot at fixed alpha", sys_pbr_fugacity; control = false)
report_differs("A-r7mu", "ruwen7 differs from e + P = T d_T P_tot + alpha d_alpha P_tot (the wrong enthalpy)", sys_pbr_fugacity; enth = :Tμ)
report("A-r6", "one_d_ideal_matrix_ruwen6 == n = d_alpha P_h / T", sys_pbr_density; control = false)
report_differs("A-r6s", "control: ruwen6 with column 3, rows 1-2 of A_t and A_x exchanged differs", sys_pbr_density; swap = true)

# ---------------------------------------------------------------- the drivers, at physical states
# The driver is called with the shipped EOS and transport; the same system is rebuilt by AD with those inputs.
const CH = F.charm_quark_list()
const SHEAR = F.QGPViscosity(0.16, 0.2); const BULK = F.SimpleBulkViscosity(0.2, 0.3)
const DIF = F.ConstDiffusion(0.116, 1.5)
const HQE = F.Heavy_Quark(hadron_list = CH)
const CEOS = F.ConformalEOS(47.5)
const PH = F.FluidProperties(HQE, SHEAR, BULK, DIF)
const PG = F.FluidProperties(F.FluiduMEoS(), SHEAR, BULK, DIF)
const PC = F.FluidProperties(CEOS, SHEAR, F.ZeroBulkViscosity(), DIF)

function physical(; T = rs(.18, .5), α = rs(-1, 1), eosm = F.FluiduMEoS(), f = 1.0, conformal = false)
    if conformal
        a = F.conformal_aSB(CEOS); p0, p1, p2 = a * T^4, 4a * T^3, 12a * T^2
    else
        th = F.thermodynamic(T, eosm); p0, p1, p2 = th.pressure, th.pressure_derivative[1], th.pressure_hessian[1]
    end
    c = F.charm_eos(T, α, CH)
    τπ = F.τ_shear(T, p1, SHEAR)
    (T = T, α = α, p0 = p0, p1 = p1, p2 = p2, η = F.viscosity(T, p1, SHEAR), τπ = τπ, dpp = 4/3 * τπ,
     ζ = conformal ? 0.0 : F.bulk_viscosity(T, p1, BULK), τΠ = conformal ? 1.0 : F.τ_bulk(T, p1, p2, BULK),
     n = c.n, nT = c.dndT, na = c.dnda, h = c.h, hT = c.dhdT, I = c.I31, IT = c.dI31dT, Ia = c.dI31da,
     aQ = c.aQ, b = c.bQ, bT = c.dbQdT, c = c.cQ, cT = c.dcQdT,
     Ds = F.DsT(DIF, T) / T / F.fmGeV * f, τn = F.τ_diffusion_hadron(T, α, HQE, DIF), κ = F.diffusion_hadron(T, α, HQE, DIF),
     a = conformal ? F.conformal_aSB(CEOS) : 1.0)
end
function run_driver(m, ϕ, X, params, τ)
    nf = length(ϕ)
    A = ntuple(_ -> zeros(MMatrix{nf,nf,Float64}), length(X)); S = zeros(MVector{nf,Float64})
    m(A, S, SVector{nf}(ϕ), τ, X, params)
    map(Matrix, A), Vector(S)
end
# worst relative difference of the driver's solved form from the AD rebuild
function cmp_driver(mine, A, S)
    w = 0.0
    for a in eachindex(A)
        X = mine[1] \ mine[a+1]; w = max(w, maximum(abs, X .- A[a]) / max(maximum(abs, X), 1e-300))
    end
    X = mine[1] \ mine[end]
    max(w, maximum(abs, X .- S) / max(maximum(abs, X), 1e-300))
end
function check_driver(id, name, sysf, build; ctl = true, N = 8, physkw = (;), syskw = (;), over = p -> (;))
    w = 0.0
    for _ in 1:N
        p = physical(; physkw...)
        OVERRIDE[] = merge(p, over(p))
        mine, _ = sysf(; syskw...)
        L = LAST[]
        A, S = build(L, p)
        w = max(w, cmp_driver(mine, A, S))
    end
    OVERRIDE[] = (;)
    chk(id, w < 1e-7, w, name)
end
bg1(L) = (τ, r) -> (L[:p].T, L[:ur], L[:g].dtT, L[:g].drT, L[:g].drur, L[:g].dtur)
bg2(L) = (τ, x, y) -> (L[:p].T, L[:ux], L[:uy], L[:g].dtT, L[:g].dxT, L[:g].dyT, L[:g].dtux, L[:g].dxux, L[:g].dyux, L[:g].dtuy, L[:g].dxuy, L[:g].dyuy)

check_driver("D-visc1", "matrix1d_visc!: EOS, eta, tau_pi, zeta, tau_Pi, dpp = 4/3 tau_pi into the kernel", sys_visc1d,
    (L, p) -> run_driver(F.matrix1d_visc!(), L[:q0], (L[:r],), PG, L[:τ]))
check_driver("D-visc2", "matrix2d_visc!: the same, derived matrix", (; kw...) -> sys_visc2d(:derived; kw...),
    (L, p) -> run_driver(F.matrix2d_visc!(), L[:q0], (L[:x0][2], L[:x0][3]), PG, L[:τ]))
check_driver("D-MIS", "matrix1d_visc_HQ!, dmn_eps = 0: kappa, tau_n from the list; the p slot holds the entropy and is unused", sys_visc1d_MIS_block,
    (L, p) -> run_driver(F.matrix1d_visc_HQ!(; dmn_eps = 0.0), L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos))
check_driver("D-MISe", "matrix1d_visc_HQ!, dmn_eps = 1e-3: only dn/dalpha of the number row, not kappa", sys_visc1d_MIS_block,
    (L, p) -> run_driver(F.matrix1d_visc_HQ!(; dmn_eps = 1e-3), L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-gamma", "matrix1d_visc_gamma_HQ!: n, kappa, tau_n at alpha = 0", sys_visc1d_gamma,
    (L, p) -> run_driver(F.matrix1d_visc_gamma_HQ!(; dmn_eps = 1e-3), L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos, α = 0.0))
check_driver("D-BGf1", "matrix1d_visc_HQ_BG!, dmn_eps = 1e-3: dn/dalpha only, not kappa", sys_BG_fick1d,
    (L, p) -> run_driver(F.matrix1d_visc_HQ_BG!(; dmn_eps = 1e-3, background_fields = bg1(L)), L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-BGc1", "matrix1d_visc_HQ_BG_consistent!, dmn_eps = 1e-3: enters dn/dalpha in both rows", sys_BG_cons1d,
    (L, p) -> run_driver(F.matrix1d_visc_HQ_BG_consistent!(; dmn_eps = 1e-3, background_fields = bg1(L)), L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-BGm1", "matrix1d_visc_HQ_BG_consistent_m2!, dmn_eps = 1e-3", sys_BG_m21d,
    (L, p) -> run_driver(F.matrix1d_visc_HQ_BG_consistent_m2!(; dmn_eps = 1e-3, background_fields = bg1(L)), L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-BGf2", "matrix2d_HQ_BG!, dmn_eps = 1e-3", sys_BG_fick2d,
    (L, p) -> run_driver(F.matrix2d_HQ_BG!(; dmn_eps = 1e-3, background_fields = bg2(L)), L[:q0], (L[:x0][2], L[:x0][3]), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-BGc2", "matrix2d_HQ_BG_consistent!, dmn_eps = 1e-3", sys_BG_cons2d,
    (L, p) -> run_driver(F.matrix2d_HQ_BG_consistent!(; dmn_eps = 1e-3, background_fields = bg2(L)), L[:q0], (L[:x0][2], L[:x0][3]), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-BGm2", "matrix2d_HQ_BG_m2!, dmn_eps = 1e-3", sys_BG_m22d,
    (L, p) -> run_driver(F.matrix2d_HQ_BG_m2!(; dmn_eps = 1e-3, background_fields = bg2(L)), L[:q0], (L[:x0][2], L[:x0][3]), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-3))
check_driver("D-live", "matrix2d_visc_HQ_BG!, dmn_eps = 0: derived medium with dpp, passive medium, live charm rows", sys_visc_HQ_upstream,
    (L, p) -> run_driver(F.matrix2d_visc_HQ_BG!(; dmn_eps = 0.0), L[:q0], (L[:x0][2], L[:x0][3]), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), syskw = (; dpp = :phys, medium = :derived))
check_driver("D-tw1", "matrix1d_two_way!, dmn_eps = 0, no edge", sys_tw1d,
    (L, p) -> run_driver(F.matrix1d_two_way!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:r],), PG, L[:τ]))
check_driver("D-tw1c", "matrix1d_two_way!, conformal", sys_tw1d,
    (L, p) -> run_driver(F.matrix1d_two_way!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:r],), PC, L[:τ]);
    physkw = (; conformal = true), syskw = (; conformal = true))
check_driver("D-tw2", "matrix2d_two_way!, dmn_eps = 0, no edge", sys_tw2d,
    (L, p) -> run_driver(F.matrix2d_two_way!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:x0][2], L[:x0][3]), PG, L[:τ]))
check_driver("D-tw2c", "matrix2d_two_way!, conformal", sys_tw2d,
    (L, p) -> run_driver(F.matrix2d_two_way!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:x0][2], L[:x0][3]), PC, L[:τ]);
    physkw = (; conformal = true), syskw = (; conformal = true))
check_driver("D-r4", "matrix1d_ideal_HQ_density!: the medium EOS into ruwen4", sys_ideal_density,
    (L, p) -> run_driver(F.matrix1d_ideal_HQ_density!, L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos))
check_driver("D-r5", "matrix1d_ideal_HQ_fugacity!: the medium EOS and n from the list into ruwen5", sys_ideal_fugacity,
    (L, p) -> run_driver(F.matrix1d_ideal_HQ_fugacity!, L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos))
# the p_br wrappers: P_h = n T of the list; rebuilt with the same P_h, its derivatives and n
function pbr_phys(p)
    c = F.charm_pressure(p.T, p.α, CH); th = F.thermodynamic(p.T, p.α, CH)
    (; I = c.pressure_derivative[1], IT = c.pressure_derivative[2], aQ = c.pressure_hessian[1], b = c.pressure_hessian[2], c = c.pressure_hessian[3],
       h = c.pressure, n = th.pressure, nT = th.pressure_derivative[1], na = th.pressure_derivative[2])
end
check_driver("D-r7", "matrix1d_ideal_HQ_br_p_fugacity!: P_h = n T, n from the list", sys_pbr_fugacity,
    (L, p) -> run_driver(F.matrix1d_ideal_HQ_br_p_fugacity!, L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = pbr_phys)
check_driver("D-r6", "matrix1d_ideal_HQ_br_p_density!: alpha from n, P_h = n T", sys_pbr_density,
    (L, p) -> run_driver(F.matrix1d_ideal_HQ_br_p_density!, L[:q0], (L[:r],), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = pbr_phys)
check_driver("D-HQup", "matrix2d_visc_HQ!: phi_8 is the fugacity, P in the pressure slot, D_s n = kappa", sys_visc_HQ_upstream,
    (L, p) -> run_driver(F.matrix2d_visc_HQ!, L[:q0], (L[:x0][2], L[:x0][3]), PH, L[:τ]);
    physkw = (; eosm = HQE.eos), over = p -> (; na = p.na + 1e-4, κ = p.n * F.DsT(DIF, p.T) / p.T / F.fmGeV, τn = F.τ_diffusion(p.T, DIF)))
# the Hessian of the list against differences of its gradient, and the two small repairs
let w = 0.0, δ = 1e-5
    for (T, α) in ((0.16, 0.5), (0.3, -2.0), (0.45, 1.0))
        th = F.thermodynamic(T, α, CH); up = F.thermodynamic(T + δ, α, CH); dn = F.thermodynamic(T - δ, α, CH); ua = F.thermodynamic(T, α + δ, CH); da = F.thermodynamic(T, α - δ, CH)
        h = th.pressure_hessian
        w = max(w, abs(h[1] / ((up.pressure_derivative[1] - dn.pressure_derivative[1]) / 2δ) - 1),
                abs(h[2] / ((up.pressure_derivative[2] - dn.pressure_derivative[2]) / 2δ) - 1),
                abs(h[3] / ((ua.pressure_derivative[2] - da.pressure_derivative[2]) / 2δ) - 1))
    end
    chk("D-hess", w < 1e-6, w, "the Hessian of a charm list == differences of its gradient")
    y = 0.0
    for (T, α) in ((0.16, 0.5), (0.3, -2.0))
        n = F.thermodynamic(T, α, CH).pressure
        y = max(y, abs(F.charm_pressure(T, α, CH).pressure_derivative[2] / T / n - 1), abs(F.charm_alpha(T, n, CH) - α))
    end
    chk("D-Ph", y < 1e-12, y, "d_alpha P_h / T == n, and charm_alpha inverts n")
    v = F.hq_pressure(0.3, 0.2; m = 1.5).pressure / F.thermodynamic(0.3, 0.2, CH).pressure - 1
    chk("D-hqp", abs(v) < 1e-14, abs(v), "hq_pressure runs: the density of one charm quark")
    ϕ5 = SVector(0.3, 0.2, 1e-3, -1e-3, 1e-3)
    m1 = F.matrix1d_visc(ϕ5, 1.3, (2.7,), PG)
    chk("D-mv", m1 !== nothing, 0.0, "matrix1d_visc runs")
end
check_driver("D-cm1", "matrix1d_two_way_cM!, dmn_eps = 0, no edge", sys_tw1d,
    (L, p) -> run_driver(F.matrix1d_two_way_cM!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:r],), PG, L[:τ]); syskw = (; cm = true))
check_driver("D-cm1c", "matrix1d_two_way_cM!, conformal", sys_tw1d,
    (L, p) -> run_driver(F.matrix1d_two_way_cM!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:r],), PC, L[:τ]);
    physkw = (; conformal = true), syskw = (; conformal = true, cm = true))
check_driver("D-cm2", "matrix2d_two_way_cM!, dmn_eps = 0, no edge", sys_tw2d,
    (L, p) -> run_driver(F.matrix2d_two_way_cM!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:x0][2], L[:x0][3]), PG, L[:τ]); syskw = (; cm = true))
check_driver("D-cm2c", "matrix2d_two_way_cM!, conformal", sys_tw2d,
    (L, p) -> run_driver(F.matrix2d_two_way_cM!(; species = CH, edge = nothing, dmn_eps = 0.0), L[:q0], (L[:x0][2], L[:x0][3]), PC, L[:τ]);
    physkw = (; conformal = true), syskw = (; conformal = true, cm = true))
# the cold edge: D_s -> f(T) D_s and nothing else
let e = F.ColdEdge()
    check_driver("D-edge", "matrix1d_two_way! at T below T_edge: D_s -> f(T) D_s and nothing else", sys_tw1d,
        (L, p) -> run_driver(F.matrix1d_two_way!(; species = CH, edge = e, dmn_eps = 0.0), L[:q0], (L[:r],), PG, L[:τ]);
        physkw = (; T = 0.145, f = F.edge_factor(0.145, e)))
    check_driver("D-edge2", "matrix2d_two_way! at T below T_edge: D_s -> f(T) D_s and nothing else", sys_tw2d,
        (L, p) -> run_driver(F.matrix2d_two_way!(; species = CH, edge = e, dmn_eps = 0.0), L[:q0], (L[:x0][2], L[:x0][3]), PG, L[:τ]);
        physkw = (; T = 0.145, f = F.edge_factor(0.145, e)))
    w = 0.0; w2 = 0.0
    for T in (0.137, 0.142, 0.15)
        f = F.edge_factor(T, e); c = F.corona_weight(T, e)
        w = max(w, abs((1 - f) / (1 - e.floor) - c))
    end
    chk("D-edgew", w < 1e-13, w, "corona_weight = (1 - f)/(1 - floor) inside the window")
    chk("D-edgeb", F.edge_factor(e.T_edge, e) == 1 && F.edge_factor(e.T_edge - e.width, e) == e.floor &&
        F.corona_weight(e.T_edge, e) == 0 && F.corona_weight(e.T_edge - e.width, e) == 1, 0.0, "f = 1 above T_edge, f = floor below T_edge - width; the weight the other way round")
end
# the vacuum regulator of the two-way drivers: u.d alpha on the number row's alpha column only
let ε = 1e-3, w = 0.0
    for _ in 1:8
        p = physical(); OVERRIDE[] = p; mine, _ = sys_tw1d(); L = LAST[]
        M = (copy(mine[1]), copy(mine[2]), mine[3])
        uτ = sqrt(1 + L[:q0][2]^2)
        M[1][6, 6] += ε * uτ; M[2][6, 6] += ε * L[:q0][2]
        A, S = run_driver(F.matrix1d_two_way!(; species = CH, edge = nothing, dmn_eps = ε), L[:q0], (L[:r],), PG, L[:τ])
        w = max(w, cmp_driver(M, A, S))
    end
    OVERRIDE[] = (;)
    chk("D-twe", w < 1e-7, w, "matrix1d_two_way!, dmn_eps = 1e-3: + eps u.d alpha on the number row's alpha column, nothing else")
end

# ---------------------------------------------------------------- kappa, tau_n, the closures, for a species list
let
    L3 = F.HadronResonaceGas_ccbar_grandcanonical(F.StructArray([
        F.Particle("D", 1.8646, 0.0, 4, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 421),
        F.Particle("Ds", 1.9683, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 431),
        F.Particle("Lc", 2.2865, 0.0, 4, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4122)]), 1.0)
    H3 = F.Heavy_Quark(hadron_list = L3)
    w1 = w2 = w3 = w4 = 0.0
    for (T, α) in ((0.16, 0.5), (0.3, -2.0), (0.45, 1.0))
        for (l, h) in ((CH, HQE), (L3, H3))
            c = F.charm_eos(T, α, l); Ds = F.DsT(DIF, T) / T / F.fmGeV
            w1 = max(w1, abs(F.diffusion_hadron(T, α, h, DIF) / (Ds * c.dnda) - 1))
            w2 = max(w2, abs(F.τ_diffusion_hadron(T, α, h, DIF) / (Ds * c.h / T) - 1))
            w3 = max(w3, abs(c.I31 / (c.n * T * c.h) - 1), abs(c.bQ / (T * c.h) - 1))
        end
        # the species sums against a momentum quadrature (composite Simpson)
        for (m, g) in ((1.5, 6),)
            N = 40000; pmax = 60T + 6
            mom(fn) = (d = pmax / N; s = 0.0;
                for i in 0:N; p_ = i * d; E = sqrt(p_^2 + m^2);
                    s += (i == 0 || i == N ? 1 : (isodd(i) ? 4 : 2)) * fn(p_, E) * exp(-E / T) end; s * d / 3)
            n0 = mom((p_, E) -> p_^2)
            c = F.charm_eos(T, α, CH)
            w4 = max(w4, abs(mom((p_, E) -> p_^4 / E^0 / 3 / E * E) / n0 / (T * c.h) - 1))        # I31 / n
            w4 = max(w4, abs(mom((p_, E) -> E^2 * p_^4 / E / 3) / n0 / T / (c.aQ) - 1))           # I41 / P0
            w4 = max(w4, abs(mom((p_, E) -> p_^6 / E / 15) / n0 / (T * c.bQ) - 1))                # I42 / P0
            w4 = max(w4, abs(mom((p_, E) -> p_^6 / 15) / mom((p_, E) -> p_^6 / E / 15) / c.cQ - 1))  # I52 / I42
        end
    end
    chk("D-kappa", w1 < 1e-11, w1, "diffusion_hadron == D_s dn/dalpha, one species and a list")
    chk("D-taun", w2 < 1e-11, w2, "tau_diffusion_hadron == D_s h/T, one species and a list")
    chk("D-I31", w3 < 1e-13, w3, "I31 = n T h, B = T h")
    chk("D-mom", w4 < 1e-7, w4, "I31/n, I41/P0, I42/P0, I52/I42 of a Boltzmann species, momentum quadrature == T h, A, B, C")
end

# ---------------------------------------------------------------- the corona update, on a grid
let
    e = F.ColdEdge(); ε = 0.0; αmax = 200.0
    rr = collect(0.25:0.5:10.25); dr = 0.5; nr = length(rr); τ = 1.7
    u = zeros(10, nr); du = randn(10, nr) .* 1e-3
    for i in 1:nr
        u[1, i] = 0.12 + 0.002 * sin(i); u[2, i] = 0.3 * sin(0.7i); u[6, i] = -1 + 0.2 * cos(i); u[7, i] = 1e-4 * sin(1.3i)
    end
    current(T, ur, α, ν) = (uτ = sqrt(1 + ur^2); n = F.thermodynamic(T, α, CH).pressure; (n * uτ + ur * ν / uτ, n * ur + ν))
    llf(UL, UR, FL, FR) = (a = min(1.0, max(UL > 0 ? abs(FL) / UL : 1.0, UR > 0 ? abs(FR) / UR : 1.0)); (FL + FR) / 2 - a * (UR - UL) / 2)
    du0 = copy(du)
    F.corona_update_1d!(du, u, τ, rr, dr, 6, CH, e, ε, αmax)
    w = 0.0; tele = 0.0
    G = zeros(nr + 1)
    for i in 2:nr-1
        (Ui, Fi) = current(u[1, i], u[2, i], u[6, i], u[7, i])
        rp = rr[i] + dr / 2
        (Uj, Fj) = current(u[1, i+1], u[2, i+1], u[6, i+1], u[7, i+1])
        G[i+1] = τ * rp * llf(Ui, Uj, Fi, Fj)
    end
    (U1, F1) = current(u[1, 1], u[2, 1], u[6, 1], u[7, 1]); (U2, F2) = current(u[1, 2], u[2, 2], u[6, 2], u[7, 2])
    G[2] = τ * (rr[2] - dr / 2) * llf(U1, U2, F1, F2)
    total = 0.0
    for i in 2:nr-1
        # d(τ r N^τ)/dτ along the updated rates
        th = F.thermodynamic(u[1, i], u[6, i], CH); nn = th.pressure; nT, na = th.pressure_derivative
        uτf(s) = sqrt(1 + (u[2, i] + s * du[2, i])^2)
        ντf(s) = (u[2, i] + s * du[2, i]) * (u[7, i] + s * du[7, i]) / uτf(s)
        Ṅ = (nT * du[1, i] + na * du[6, i]) * uτf(0.0) + nn * ForwardDiff.derivative(uτf, 0.0) + ForwardDiff.derivative(ντf, 0.0)
        rate = rr[i] * (current(u[1, i], u[2, i], u[6, i], u[7, i])[1] + τ * Ṅ)
        target = -(G[i+1] - G[i]) / dr
        w = max(w, abs(rate - target) / max(abs(target), 1e-30))
        total += rate * dr
    end
    tele = abs(total + (G[nr] - G[2])) / max(abs(G[nr] - G[2]), 1e-30)
    chk("D-cor", w < 1e-9, w, "corona_update_1d!, T below the window: d(tau r N^tau)/dtau == -(G_{i+1/2} - G_{i-1/2})/dr, LLF, with the updated alpha rate")
    chk("D-tele", tele < 1e-9, tele, "the flux form telescopes: the sum of r dr dN/dtau over the cells is the boundary flux difference")
    # blend: inside the window the rate is (1 - w) quasi-linear + w flux form; above the window nothing changes
    u2 = copy(u); u2[1, :] .= 0.2; du2 = copy(du0); F.corona_update_1d!(du2, u2, τ, rr, dr, 6, CH, e, ε, αmax)
    chk("D-cor0", du2 == du0, 0.0, "above T_edge the update changes nothing")
    u3 = copy(u); u3[1, :] .= 0.146; du3 = copy(du0); F.corona_update_1d!(du3, u3, τ, rr, dr, 6, CH, e, ε, αmax)
    u4 = copy(u3); du4 = copy(du0); F.corona_update_1d!(du4, u4, τ, rr, dr, 6, CH, F.ColdEdge(; T_edge = 0.2, width = 0.01), ε, αmax)
    wgt = F.corona_weight(0.146, e); bl = maximum(abs, (du3[6, 2:nr-1] .- ((1 - wgt) .* du0[6, 2:nr-1] .+ wgt .* du4[6, 2:nr-1])))
    chk("D-corb", bl < 1e-14 * maximum(abs, du4[6, :]) + 1e-16 && maximum(abs, du4[6, 2:nr-1] .- du0[6, 2:nr-1]) > 0, bl, "inside the window: du_alpha = (1 - w) quasi-linear + w flux form")
    # 2+1D, with the 14-field layout: alpha = 8, nu = 9, 10
    nx, ny = 9, 9; dx, dy = 0.5, 0.5
    v = zeros(14, nx, ny); dv = randn(14, nx, ny) .* 1e-3
    for j in 1:ny, i in 1:nx
        v[1, i, j] = 0.12 + 0.002 * sin(i + j); v[2, i, j] = 0.2 * sin(0.7i); v[3, i, j] = 0.2 * cos(0.6j)
        v[8, i, j] = -1 + 0.2 * cos(i - j); v[9, i, j] = 1e-4 * sin(1.3j); v[10, i, j] = 1e-4 * cos(1.1i)
    end
    cur2(i, j) = (ux = v[2, i, j]; uy = v[3, i, j]; uτ = sqrt(1 + ux^2 + uy^2); n = F.thermodynamic(v[1, i, j], v[8, i, j], CH).pressure;
                  (n * uτ + (ux * v[9, i, j] + uy * v[10, i, j]) / uτ, n * ux + v[9, i, j], n * uy + v[10, i, j]))
    F.corona_update_2d!(dv, v, τ, dx, dy, 8, CH, e, ε, αmax)
    w2 = 0.0
    for j in 2:ny-1, i in 2:nx-1
        U, X, Y = cur2(i, j); UL, XL, _ = cur2(i - 1, j); UR, XR, _ = cur2(i + 1, j); UD, _, YD = cur2(i, j - 1); UU, _, YU = cur2(i, j + 1)
        target = -τ * ((llf(U, UR, X, XR) - llf(UL, U, XL, X)) / dx + (llf(U, UU, Y, YU) - llf(UD, U, YD, Y)) / dy)
        th = F.thermodynamic(v[1, i, j], v[8, i, j], CH); nn = th.pressure; nT, na = th.pressure_derivative
        uτf(s) = sqrt(1 + (v[2, i, j] + s * dv[2, i, j])^2 + (v[3, i, j] + s * dv[3, i, j])^2)
        ντf(s) = ((v[2, i, j] + s * dv[2, i, j]) * (v[9, i, j] + s * dv[9, i, j]) + (v[3, i, j] + s * dv[3, i, j]) * (v[10, i, j] + s * dv[10, i, j])) / uτf(s)
        rate = U + τ * ((nT * dv[1, i, j] + na * dv[8, i, j]) * uτf(0.0) + nn * ForwardDiff.derivative(uτf, 0.0) + ForwardDiff.derivative(ντf, 0.0))
        w2 = max(w2, abs(rate - target) / max(abs(target), 1e-30))
    end
    chk("D-cor2", w2 < 1e-9, w2, "corona_update_2d!, T below the window: d(tau N^tau)/dtau == the LLF flux differences")
end
println("evaluated ", NCHK[], ", failed ", NFAIL[])
exit(NFAIL[] == 0 ? 0 : 1)
