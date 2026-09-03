#!/usr/bin/env julia
#=
02 — how fine a transverse grid do I need?  The honest answer is set by NUMERICAL viscosity.

The spatial scheme is first-order upwind.  A first-order upwind scheme does not merely converge
slowly, it adds a diffusivity of its own, and in a VISCOUS hydrodynamics code that diffusivity is
indistinguishable from the physics you are trying to measure.  This file measures it, on the one
configuration where the exact answer is available in closed form.

THE MEASUREMENT.  Linearising Israel-Stewart about a static uniform state gives, for a plane sound
wave along x,

    omega^2 = k^2 c_s^2 - i omega k^2 Z(omega)/w,   Z = (4/3) eta/(1 - i omega tau_pi) + zeta/(1 - i omega tau_Pi)

with w = e + P.  Seeding the exact right-moving eigenvector and fitting the decay of the amplitude
gives Im omega.  Whatever exceeds the exact root is the scheme's own dissipation, and this file
shows it is FIRST ORDER in Delta x, which is what identifies it as numerical rather than as a
missing term in the equations.

THE ANSWER, for eta/s = 0.1 at T = 0.3 GeV:  D_num ~ 0.75 Delta x, against a physical
Gamma_s = ((4/3) eta + zeta)/(e+P) = 0.090 fm.  The two are EQUAL at Delta x ~ 0.12 fm.  Above that
the transport in the run is mostly the grid.

⚠ THIS IS WHY A CODE-TO-CODE COMPARISON MUST BE RUN AT SEVERAL RESOLUTIONS.  The Fluidum-vs-FiVo
2+1D comparison (Julia/Projects/FiVoFluidumComparison/COMPARISON_2P1D.md) took its finest grid at
Delta x = 0.15 fm, i.e. just inside the regime where this solver's numerical viscosity exceeds its
physical one.  Its conclusion — first-order convergence of the two codes to each other — is exactly
the right form of statement to make there, and a single-resolution number would not have been.

⚠ TAU_0 = 200 fm/c IS DELIBERATE.  The Milne source terms go as 1/tau; at tau = 200 they are ~0.5 %
of omega, so the measurement is of the flat-space dispersion relation.  The amplitude is taken from
the DIFFERENCE of a perturbed and an unperturbed run, which removes the residual background drift.
=#
ENV["GKSwstype"] = "100"
using Fluidum, Printf, Plots
gr(); default(; fontfamily = "sans-serif", framestyle = :box, grid = false, dpi = 150, lw = 2)
const F   = Fluidum
const FIG = joinpath(@__DIR__, "figures"); isdir(FIG) || mkpath(FIG)

const EOS = FluiduMEoS()
const PAR = F.FluidProperties(EOS, QGPViscosity(0.1, 0.2), SimpleBulkViscosity(0.083, 15.0), ZeroDiffusion())
ent(T) = F.pressure_derivative(T, Val(1), EOS); dst(T) = F.pressure_derivative(T, Val(2), EOS)
enth(T) = T*ent(T); cs2(T) = ent(T)/(T*dst(T))
function coeffs(T)
    th = thermodynamic(T, EOS)
    (F.viscosity(T, th, PAR.shear), F.τ_shear(T, th, PAR.shear),
     F.bulk_viscosity(T, th, PAR.bulk), F.τ_bulk(T, th, PAR.bulk))
end
"the exact MIS sound root, by Newton on the complex omega"
function mis_omega(k, T)
    eta, taupi, zet, tauPi = coeffs(T); w = enth(T); c2 = cs2(T)
    Z(om) = (4/3)*eta/(1 - im*om*taupi) + zet/(1 - im*om*tauPi)
    f(om) = om^2 - k^2*c2 + im*om*k^2*Z(om)/w
    om = sqrt(c2)*k - 0.5im*k^2*((4/3)*eta + zet)/w
    for _ in 1:200
        h = 1e-7*max(abs(om), 1e-3); om -= f(om)/((f(om+h) - f(om-h))/(2h))
    end
    om
end
"the right-moving eigenvector, normalised to delta e = 1"
function eigvec(k, T)
    om = mis_omega(k, T); eta, taupi, zet, tauPi = coeffs(T); w = enth(T)
    v = om/(k*w); pxx = -(4/3)*eta*(im*k*v)/(1 - im*om*taupi)
    (om, 1/(T*dst(T)), v, -pxx/2, -pxx/2, -zet*(im*k*v)/(1 - im*om*tauPi))
end
function measure(T0, k, L, N, tau0, tau1; xmax = 12.0, nout = 17, amp = 1e-4)
    om, dT, v, pyy, pzz, PiB = eigvec(k, T0); s = amp*T0/abs(dT)
    d = DiscreteFields(F.viscous_2d(),
        CartesianDiscretization(F.SymmetricInterval(N, L), F.SymmetricInterval(4, 2.0)), Float64)
    taus = collect(range(tau0, tau1; length = nout)); res = Dict()
    for pert in (false, true)
        q = pert ? s : 0.0
        p = set_array((x,y) -> T0 + q*real(dT*exp(im*k*x)), :temperature, d)
        set_array!(p, (x,y) -> q*real(v  *exp(im*k*x)), :ux,   d)
        set_array!(p, (x,y) -> q*real(pyy*exp(im*k*x)), :piyy, d)
        set_array!(p, (x,y) -> q*real(pzz*exp(im*k*x)), :pizz, d)
        set_array!(p, (x,y) -> q*real(PiB*exp(im*k*x)), :piB,  d)
        set_array!(p, (x,y) -> 0.0, :uy, d); set_array!(p, (x,y) -> 0.0, :pixy, d)
        sl = F.oneshoot(d, F.matrix2d_visc!, PAR, p, (tau0, tau1); reltol = 1e-11, abstol = 1e-13)
        for t in taus; res[(pert, t)] = sl(t); end
    end
    xs = [-L + (i-0.5)*2L/N for i in 1:N]; sel = findall(x -> abs(x) <= xmax, xs)
    a = ComplexF64[]
    for t in taus
        dlt = [res[(true,t)][1,i+1,3] - res[(false,t)][1,i+1,3] for i in 1:N]
        push!(a, 2*sum(dlt[i]*exp(-im*k*xs[i]) for i in sel)/length(sel))
    end
    lm = log.(abs.(a)); ph = angle.(a)
    for i in 2:length(ph)
        while ph[i]-ph[i-1] >  pi; ph[i] -= 2pi; end
        while ph[i]-ph[i-1] < -pi; ph[i] += 2pi; end
    end
    X = hcat(ones(nout), taus .- tau0)
    (-(X\ph)[2], -(X\lm)[2], om, taus, a)
end

const T0 = 0.30; const LAM = 6.0; const K = 2pi/LAM
eta, taupi, zet, tauPi = coeffs(T0); Gs = ((4/3)*eta + zet)/enth(T0)
om = mis_omega(K, T0)
@printf("\n  T = %.2f GeV, lambda = %.1f fm (k = %.4f 1/fm), eta/s = 0.1\n", T0, LAM, K)
@printf("  exact MIS root      omega = %.6f - %.6f i\n", real(om), -imag(om))
@printf("  Navier-Stokes limit omega = %.6f - %.6f i   (omega tau_pi = %.2f, so NS is already %.0f %% off)\n",
        sqrt(cs2(T0))*K, 0.5*K^2*Gs, abs(real(om))*taupi, 100*abs(0.5*K^2*Gs + imag(om))/abs(imag(om)))
@printf("  physical Gamma_s = ((4/3)eta + zeta)/(e+P) = %.4f fm\n\n", Gs)

@printf("  %6s %8s %12s %10s %12s %12s %8s\n", "N", "dx [fm]", "Re omega", "rel", "Im omega", "excess", "order")
Ns = (100, 200, 400, 800); ei = Float64[]; er = Float64[]; series = []
for N in Ns
    wr, wi, _, taus, a = measure(T0, K, 20.0, N, 200.0, 208.0)
    push!(er, abs(wr - real(om))); push!(ei, wi + imag(om)); push!(series, (N, taus, a))
    @printf("  %6d %8.4f %12.6f %10.2e %12.6f %12.6f %8s\n", N, 40.0/N, wr, er[end]/abs(real(om)), wi, ei[end],
            length(ei) < 2 ? "-" : @sprintf("%.2f", ei[end-1]/ei[end]))
end
Dn = 2*ei[end]/K^2; dx = 40.0/Ns[end]
@printf("\n  numerical viscosity  D_num = %.4f fm at dx = %.3f fm, i.e. %.2f * dx\n", Dn, dx, Dn/dx)
@printf("  it equals the physical Gamma_s = %.4f fm at dx = %.3f fm and DOMINATES above that.\n", Gs, Gs/(Dn/dx))
@printf("  Richardson from the two finest grids: Im omega -> %.6f against the exact %.6f (%.1f %%)\n",
        -imag(om) + 2ei[end] - ei[end-1], -imag(om), 100*abs(2ei[end] - ei[end-1])/abs(imag(om)))

plt = plot(layout = (1,2), size = (1040, 400))
for (N, taus, a) in series
    plot!(plt[1], taus .- 200.0, abs.(a)./abs(a[1]); label = "N = $N", yscale = :log10,
          xlabel = "tau - tau_0  [fm/c]", ylabel = "amplitude (normalised)", title = "sound decay")
end
plot!(plt[1], series[1][2] .- 200.0, exp.(imag(om) .* (series[1][2] .- 200.0));
      label = "exact MIS", ls = :dash, lc = :black, lw = 3)
dxs = [40.0/N for N in Ns]
plot!(plt[2], dxs, ei; marker = :circle, xscale = :log10, yscale = :log10, label = "measured excess",
      xlabel = "dx  [fm]", ylabel = "excess damping  [1/fm]", title = "it is first order => numerical")
plot!(plt[2], dxs, (0.5*K^2*(Dn/dx)) .* dxs;              # the first-order law, one fitted slope
      ls = :dash, lc = :black, label = @sprintf("0.5 k^2 * %.2f dx", Dn/dx))
savefig(plt, joinpath(FIG, "ex02_numerical_viscosity.png"))
println("\n  -> ", joinpath(FIG, "ex02_numerical_viscosity.png"))
println("""
  READING IT
    The excess halves when dx halves. A MISSING PHYSICAL TERM would not do that -- it would sit at a
    fixed size as the grid refines. That is the whole diagnostic: refine, and see whether the
    discrepancy follows the grid or stays put.
    Practical consequence: at dx >= 0.12 fm this solver's effective eta/s at T = 0.3 GeV is roughly
    DOUBLE the one you asked for. Either refine, or quote the resolution alongside every transport
    coefficient, or -- best -- report the convergence, as the code-to-code comparison did.""")
