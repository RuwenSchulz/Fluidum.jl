# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, mu, nur

function matrix1d_visc_HQ_BG_second_moment!(;dmn_eps=1e-6, background_fields)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_second_moment!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
function matrix1d_visc_HQ_BG_second_moment!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-6,background_fields= nothing)

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


using LinearAlgebra
using ForwardDiff
using StaticArrays

"""
    local_qnm_spectrum_BG_second_moment_robust(
        ϕ, tau, r, params;
        background_fields,
        k=0.0,
        dmn_eps=1e-10,
        imag_tol=1e-12,
        svd_reg=0.0,
        verbose=false,
    )

Robust local frozen-coefficient QNM spectrum for the BG second-moment system.

Conventions:
- PDE linearization is treated locally as
      A_t ∂_τ δX + A_r ∂_r δX = J δX
- Fourier ansatz:
      δX ~ exp(s τ + i k r)
  gives
      (J - i k A_r) v = s A_t v
- We therefore define
      L(k) = A_t^{-1} (J - i k A_r)
  and compute eigvals(L).

Returned fields:
- `:modes`             vector of NamedTuples `(s, gamma, omega, stable, residual)`
- `:eigvals`
- `:eigvecs`
- `:At`, `:Ar`, `:J`
- `:source0`
- `:cond_At`
- `:sigma_min_At`
- `:det_At`
- `:rank_J`
- `:background`
"""
function local_qnm_spectrum_BG_second_moment_robust(
    ϕ::AbstractVector,
    tau::Real,
    r::Real,
    params;
    background_fields,
    k::Real = 0.0,
    dmn_eps::Real = 1e-10,
    imag_tol::Real = 1e-12,
    svd_reg::Real = 0.0,
    verbose::Bool = false,
)
    x0 = collect(float.(ϕ))
    @assert length(x0) == 3 "Expected 3-component state X = [X1,X2,X3]"

    T, ur, dtT, drT, drur, dtur = background_fields(tau, r)

    # ------------------------------------------------------------
    # coefficient builder with full dependence on x[1]
    # ------------------------------------------------------------
    function build_mats_and_source(x::AbstractVector)
        thermo = thermodynamic(T, x[1], params.eos.hadron_list)
        n = thermo.pressure
        dn_dT, dn_dmu = thermo.pressure_derivative
        dn_dmu += dmn_eps

        kappa = diffusion_hadron(T, x[1], params.eos, params.diffusion)
        tauDiffusion = τ_diffusion_hadron(T, x[1], params.eos, params.diffusion)

        tau_phi = tauDiffusion / 2
        etaPhi = T * tau_phi

        At, Ar, src = one_d_viscous_matrix_fugacity_BG_only_second_moment(
            x, tau, r, ur, dtT, drT, drur, dtur,
            n, dn_dmu, dn_dT, tauDiffusion, kappa, tau_phi, etaPhi
        )

        return At, Ar, src, (
            n = n,
            dn_dT = dn_dT,
            dn_dmu = dn_dmu,
            kappa = kappa,
            tauDiffusion = tauDiffusion,
            tau_phi = tau_phi,
            etaPhi = etaPhi,
        )
    end

    At0, Ar0, source0, coeffs0 = build_mats_and_source(x0)

    # ------------------------------------------------------------
    # Jacobian of the source term with full x[1]-dependence
    # ------------------------------------------------------------
    fsource(x) = begin
        _, _, src, _ = build_mats_and_source(x)
        collect(src)
    end
    J = ForwardDiff.jacobian(fsource, x0)

    # Optional: also inspect how At/Ar vary with x if needed later
    # fAt(x) = vec(Matrix(first(build_mats_and_source(x))))
    # fAr(x) = vec(Matrix(build_mats_and_source(x)[2]))

    Atm = Matrix(At0)
    Arm = Matrix(Ar0)

    # ------------------------------------------------------------
    # diagnostics
    # ------------------------------------------------------------
    S = svd(Atm)
    sigma_min_At = minimum(S.S)
    cond_At = sigma_min_At > 0 ? maximum(S.S) / sigma_min_At : Inf
    det_At = det(Atm)
    rank_J = rank(J)

    # ------------------------------------------------------------
    # robust build of L = At^{-1}(J - i k Ar)
    # optional tiny regularization if At is close to singular
    # ------------------------------------------------------------
    RHS = ComplexF64.(J) .- (im * k) .* ComplexF64.(Arm)

    L = if svd_reg > 0
        (ComplexF64.(Atm) .+ svd_reg * I) \ RHS
    else
        ComplexF64.(Atm) \ RHS
    end

    F = eigen(L)
    eigvals_raw = collect(F.values)
    eigvecs = Matrix(F.vectors)

    # ------------------------------------------------------------
    # clean tiny imaginary parts
    # ------------------------------------------------------------
    eigvals = ComplexF64[
        abs(imag(λ)) < imag_tol ? complex(real(λ), 0.0) : λ
        for λ in eigvals_raw
    ]

    # residual diagnostic: ||L v - λ v|| / ||v||
    function eig_residual(L, λ, v)
        nv = norm(v)
        nv == 0 && return Inf
        return norm(L * v .- λ .* v) / nv
    end

    modes = NamedTuple[]
    for i in eachindex(eigvals)
        λ = eigvals[i]
        v = eigvecs[:, i]
        γ = -real(λ)
        ω = imag(λ)
        push!(modes, (
            s = λ,
            gamma = γ,
            omega = ω,
            stable = (γ >= 0),
            residual = eig_residual(L, λ, v),
        ))
    end

    # ------------------------------------------------------------
    # robust sorting:
    # 1) stable modes first
    # 2) then by increasing gamma (slowest decay first)
    # 3) then by |omega|
    # 4) then by residual
    # ------------------------------------------------------------
    sort!(modes, by = m -> (
        m.stable ? 0 : 1,
        m.gamma,
        abs(m.omega),
        m.residual
    ))

    if verbose
        println("--------------------------------------------------")
        println("local_qnm_spectrum_BG_second_moment_robust")
        println("tau = ", tau, ", r = ", r, ", k = ", k)
        println("background: T=$T ur=$ur dtT=$dtT drT=$drT drur=$drur dtur=$dtur")
        println("coeffs: n=$(coeffs0.n), dn_dmu=$(coeffs0.dn_dmu), kappa=$(coeffs0.kappa), tauDiff=$(coeffs0.tauDiffusion)")
        println("cond(At) = ", cond_At)
        println("sigma_min(At) = ", sigma_min_At)
        println("det(At) = ", det_At)
        println("rank(J) = ", rank_J)
        println("eigvals(L) = ", eigvals)
        println("modes = ", modes)
        println("--------------------------------------------------")
    end

    return Dict(
        :modes => modes,
        :eigvals => eigvals,
        :eigvecs => eigvecs,
        :At => At0,
        :Ar => Ar0,
        :J => J,
        :source0 => source0,
        :cond_At => cond_At,
        :sigma_min_At => sigma_min_At,
        :det_At => det_At,
        :rank_J => rank_J,
        :k => k,
        :background => (
            T = T, ur = ur, dtT = dtT, drT = drT, drur = drur, dtur = dtur
        ),
        :coeffs => coeffs0,
    )
end