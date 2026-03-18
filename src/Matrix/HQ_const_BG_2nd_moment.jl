# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, mu, nur

function matrix1d_visc_HQ_BG_second_moment!(;dmn_eps=1e-6, background_fields)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_second_moment!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
function matrix1d_visc_HQ_BG_second_moment!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-10,background_fields= nothing)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])
    
    thermo = thermodynamic(T,ϕ[1],params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu+= dmn_eps

    kappa = diffusion_hadron(T,ϕ[1],params.eos,params.diffusion) #diffusion coefficient for hadrons
    tauDiffusion=τ_diffusion_hadron(T,ϕ[1],params.eos,params.diffusion) #tau diffusion for hadrons
    
    #tau_phi = tauDiffusion/2

    tau_phi = tauDiffusion/2
    etaPhi =  T * tau_phi

    #etaPhi  = T * tau_phi
    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )
    #(At,Ax, source)=one_d_viscous_HQ_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix_fugacity_BG_only_second_moment(ϕ,tau,X[1],ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,kappa,tau_phi,etaPhi)

    Ainv= inv(At)

    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

    end 



function one_d_viscous_matrix_fugacity_BG_only_second_moment(X,tau,r,ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,kappa,tauPhi,etaPhi)
    #these are the matrices from 1d viscous hydro
    At=SMatrix{3,3}(
((1  .+ ur .^2) .^2 .*r .*tau .*dn_dmu,

ur .*(1  .+ ur .^2) .^2 .*kappa .*r .*tau,

0,

ur .*(1  .+ ur .^2) .*r .*tau,

(1  .+ ur .^2) .^2 .*tauDiffusion .*r .*tau,

(4 .*etaPhi .*ur .*(1  .+ ur .^2) .^2 .*r .*tau) ./3.,

0,

0,

(tauPhi .*(1  .+ ur .^2) .^3 .*(3  .+ 2 .*ur .^2) .*r .*tau) ./3.)
    )
        #########################################################################################################################################################################
        
    Ax=SMatrix{3,3}(
(ur .*(1  .+ ur .^2) .^1.5 .*r .*tau .*dn_dmu,

(1  .+ ur .^2) .^2.5 .*kappa .*r .*tau,

0,

(1  .+ ur .^2) .^1.5 .*r .*tau,

ur .*(1  .+ ur .^2) .^1.5 .*tauDiffusion .*r .*tau,

(4 .*etaPhi .*(1  .+ ur .^2) .^2.5 .*r .*tau) ./3.,

0,

(1  .+ ur .^2) .^2.5 .*r .*tau,

(tauPhi .*ur .*(1  .+ ur .^2) .^2.5 .*(3  .+ 2 .*ur .^2) .*r .*tau) ./3.)

    )
        #########################################################################################################################################################################
    
    source=SVector{3}(
(X[2] .*(sqrt.(1  .+ ur .^2) .*tau  .+ ur .^2 .*sqrt.(1  .+ ur .^2) .*tau  .+ r .*(ur  .+ ur .^3  .+ dtur .*tau))  .+ (1  .+ ur .^2) .*n .*(ur .*sqrt.(1  .+ ur .^2) .*tau  .+ r .*(1  .+ ur .^2  .+ dtur .*ur .*tau  .+ drur .*sqrt.(1  .+ ur .^2) .*tau))  .+ (1  .+ ur .^2) .*(dtT  .+ dtT .*ur .^2  .+ drT .*ur .*sqrt.(1  .+ ur .^2)) .*r .*tau .*dn_dT,

r .*(2 .*ur .*(1  .+ ur .^2) .^2 .*X[3]  .+ (sqrt.(1  .+ ur .^2)  .- ur .*(dtur .*tauDiffusion  .+ dtur .*ur .^2 .*tauDiffusion  .+ ur .*sqrt.(1  .+ ur .^2) .*( .-1  .+ drur .*tauDiffusion))) .*X[2] .*tau),

 .-0.3333333333333333 .*((1  .+ ur .^2) .*(4 .*etaPhi .*ur .*(dtur .*ur  .+ drur .*sqrt.(1  .+ ur .^2)) .*r .*X[2] .*tau  .- 3 .*sqrt.(1  .+ ur .^2) .*r .*X[3] .*tau  .+ 2 .*tauPhi .*(1  .+ ur .^2) .*X[3] .*(2 .*(1  .+ ur .^2) .*r  .- ur .*sqrt.(1  .+ ur .^2) .*tau))))
    )
    return (At,Ax, source)
end


"""
    local_qnm_spectrum_BG_second_moment(ϕ, tau, r, params; background_fields, k=0.0, dmn_eps=1e-10)

Compute *local* quasi-normal modes (QNM) of the 1D HQ background-only second-moment system
defined by `one_d_viscous_matrix_fugacity_BG_only_second_moment`.

We interpret the PDE in linearized form

    A_t ∂_τ δX + A_r ∂_r δX = J δX,

where `J = ∂(source)/∂X` evaluated at `X=ϕ` and the chosen background.
For Fourier perturbations δX ∝ exp(s τ + i k r), the modes satisfy

    (J - i k A_r) v = s A_t v.

Returns a `Dict` with fields:
- `:modes`: vector of NamedTuples `(s, omega, gamma)` sorted by increasing `gamma` (slowest first)
- `:At`, `:Ar`, `:J`, `:source0`, `:k`, `:background`
"""
function local_qnm_spectrum_BG_second_moment(
    ϕ::AbstractVector,
    tau::Real,
    r::Real,
    params;
    background_fields,
    k::Real=0.0,
    dmn_eps::Real=1e-10,
)
    T, ur, dtT, drT, drur, dtur = background_fields(tau, r)

    thermo = thermodynamic(T, ϕ[1], params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu += dmn_eps

    kappa = diffusion_hadron(T, ϕ[1], params.eos, params.diffusion)
    tauDiffusion = τ_diffusion_hadron(T, ϕ[1], params.eos, params.diffusion)

    tau_phi = tauDiffusion / 2
    etaPhi = T * tau_phi

    x0 = collect(float.(ϕ))

    # Build matrices at x0 (At, Ar) and source(x0)
    At, Ar, source0 = one_d_viscous_matrix_fugacity_BG_only_second_moment(
        x0, tau, r, ur, dtT, drT, drur, dtur,
        n, dn_dmu, dn_dT, tauDiffusion, kappa, tau_phi, etaPhi
    )

    # Jacobian of the source term w.r.t. the state X
    fsource(x) = begin
        _, _, src = one_d_viscous_matrix_fugacity_BG_only_second_moment(
            x, tau, r, ur, dtT, drT, drur, dtur,
            n, dn_dmu, dn_dT, tauDiffusion, kappa, tau_phi, etaPhi
        )
        collect(src)
    end
    J = ForwardDiff.jacobian(fsource, x0)

    # Local Fourier-mode generator L(k) = At^{-1} (J - i k Ar)
    L = Matrix(At) \ (ComplexF64.(J) .- (im * k) .* ComplexF64.(Matrix(Ar)))
    s = eigvals(L)

    modes = [(s = sk, omega = imag(sk), gamma = -real(sk)) for sk in s]
    # Sort by slowest decay (smallest positive gamma) first; push growing modes to the end.
    sort!(modes, by=m -> (m.gamma < 0 ? 1e9 : m.gamma, -abs(m.omega)))

    return Dict(
        :modes => modes,
        :At => At,
        :Ar => Ar,
        :J => J,
        :source0 => source0,
        :k => k,
        :background => (T=T, ur=ur, dtT=dtT, drT=drT, drur=drur, dtur=dtur),
    )
end
