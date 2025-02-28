using LinearAlgebra
using Statistics
using Distributions
using Plots


struct ModelParams
    β::Float64          
    γ
    ::Float64          
    α::Float64          
    δ::Float64          
    A::Float64        
    ρ::Float64          
    σz::Float64         
    τ::Float64          
    λ::Float64          
    amin::Float64      
    amax::Float64   
    nA::Int            
    nZ::Int             
end

function calibrate_baseline()
   

    β  = 0.96      
    γ  = 2.0        
    α  = 0.36
    δ  = 0.08
    A  = 1.12       
    ρ  = 0.9
    σz = 0.4
    τ  = 0.05        
    λ  = 0.0         

    amin = 0.0
    amax = 50.0
    nA   = 200

    nZ   = 5

    return ModelParams(β, γ, α, δ, A, ρ, σz, τ, λ, amin, amax, nA, nZ)
end


function tauchen_discretize(ρ, σ, nZ; m=3)
    σz = σ / sqrt(1 - ρ^2)
    zmax =  m*σz
    zmin = -m*σz
    step = (zmax - zmin)/(nZ-1)
    zgrid_ln = [zmin + step*(i-1) for i in 1:nZ]


    P = zeros(nZ,nZ)
    for j in 1:nZ
        for k in 1:nZ
            mid_up   = (k < nZ) ? (zgrid_ln[k] + zgrid_ln[k+1]) / 2 : 1e10
            mid_down = (k > 1)  ? (zgrid_ln[k] + zgrid_ln[k-1]) / 2 : -1e10
            mean_ln  = ρ*zgrid_ln[j]
            P[j,k]   = cdf(Normal(mean_ln, σ), mid_up) - cdf(Normal(mean_ln, σ), mid_down)
        end
    end

    zgrid = exp.(zgrid_ln)

    πz = stationary_dist_of_z(P)
    meanz = dot(zgrid, πz)
    zgrid .*= (1 / meanz)

    return zgrid, P
end

function stationary_dist_of_z(P; tol=1e-14, maxiter=10_000)
    nZ = size(P,1)
    dist = fill(1/nZ, nZ)
    dist_new = similar(dist)
    for _ in 1:maxiter
        dist_new = dist * P
        if maximum(abs.(dist_new .- dist)) < tol
            return dist_new
        end
        dist = dist_new
    end
    return dist
end


function asset_grid(params::ModelParams)
    return collect(range(params.amin, params.amax, length=params.nA))
end


function post_tax_income(y, ybar, τ, λ)
    return (1 - τ)^(1 - λ) * (y^λ) * (ybar)^(1 - λ)
end

Tax paid T(y) = y - post_tax_income(y).

function tax_paid(y, ybar, τ, λ)
    return y - post_tax_income(y, ybar, τ, λ)
end

function CRRA(c, γ)
    if c <= 0
        return -1.0e15
    elseif abs(γ - 1.0) < 1e-8
        return log(c)
    else
        return c^(1 - γ) / (1 - γ)
    end
end

function solve_household(params::ModelParams, r, w, zgrid, P)
    Agrid = asset_grid(params)
    nA, nZ = params.nA, params.nZ

    ybar = w
    netlab = [post_tax_income(w*z, ybar, params.τ, params.λ) for z in zgrid]

    V  = zeros(nA, nZ)
    Vn = similar(V)
    policy_idx = fill(1, nA, nZ)

    β, γ = params.β, params.γ
    R = 1 + r

    maxiter = 1000
    tol     = 1e-6

    for iter in 1:maxiter
        for iz in 1:nZ
            for ia in 1:nA
                income = netlab[iz] + R*Agrid[ia]
                best_val = -1e15
                best_aidx = 1
                for iap in 1:nA
                    a_next = Agrid[iap]
                    c = income - a_next
                    val = CRRA(c, γ)
                    if val < -1e14
                        # c <= 0 => skip
                        continue
                    end
                    # continuation
                    EV = 0.0
                    for izp in 1:nZ
                        EV += P[iz, izp]*V[iap, izp]
                    end
                    val += β*EV
                    if val > best_val
                        best_val = val
                        best_aidx = iap
                    end
                end
                Vn[ia, iz] = best_val
                policy_idx[ia, iz] = best_aidx
            end
        end
        diff = maximum(abs.(Vn .- V))
        V .= Vn
        if diff < tol
            break
        end
    end

    return V, policy_idx
end


function stationary_distribution(params::ModelParams, policy_idx, P)
    Agrid = asset_grid(params)
    nA, nZ = params.nA, params.nZ

    dist = fill(1.0/(nA*nZ), nA, nZ)
    dist_new = similar(dist)

    tol = 1e-12
    maxiter = 10_000

    for _ in 1:maxiter
        dist_new .= 0.0
        for iz in 1:nZ
            for ia in 1:nA
                iap = policy_idx[ia, iz]
                mass = dist[ia, iz]
                for izp in 1:nZ
                    dist_new[iap, izp] += mass * P[iz, izp]
                end
            end
        end
        if maximum(abs.(dist_new .- dist)) < tol
            return dist_new
        end
        dist .= dist_new
    end
    return dist
end

function equilibrium_error(Kguess::Float64, params::ModelParams, zgrid, P)
    r = params.α * params.A * Kguess^(params.α - 1) - params.δ
    w = (1 - params.α) * params.A * Kguess^(params.α)

    _, policy_idx = solve_household(params, r, w, zgrid, P)
    dist = stationary_distribution(params, policy_idx, P)

    Agrid = asset_grid(params)
    impliedK = 0.0
    for iz in 1:params.nZ
        for ia in 1:params.nA
            impliedK += Agrid[ia] * dist[ia, iz]
        end
    end

    return impliedK - Kguess
end

function find_equilibrium_K(params::ModelParams, zgrid, P;
                            Kmin=0.01, Kmax=200.0,
                            tol=1e-4, maxiter=100)
    for i in 1:maxiter
        Kmid = 0.5*(Kmin + Kmax)
        err  = equilibrium_error(Kmid, params, zgrid, P)
        if abs(err) < tol
            return Kmid
        elseif err > 0
            Kmin = Kmid
        else
            Kmax = Kmid
        end
    end
    return 0.5*(Kmin + Kmax)
end


function revenue_error(tau_guess::Float64, params::ModelParams, zgrid, P, Gtarget::Float64)
    old_tau = params.τ
    params.τ = tau_guess

    Kstar = find_equilibrium_K(params, zgrid, P)

    r = params.α * params.A * Kstar^(params.α - 1) - params.δ
    w = (1 - params.α) * params.A * Kstar^(params.α)

    _, policy_idx = solve_household(params, r, w, zgrid, P)
    dist = stationary_distribution(params, policy_idx, P)

    Agrid = asset_grid(params)
    total_revenue = 0.0
    for iz in 1:params.nZ
        for ia in 1:params.nA
            mass = dist[ia, iz]
            y = w*zgrid[iz]
            total_revenue += tax_paid(y, w, tau_guess, params.λ) * mass
        end
    end

    params.τ = old_tau  # restore
    return total_revenue - Gtarget
end

function find_tau_for_revenue(params::ModelParams, zgrid, P, Gtarget;
                              tmin=0.0, tmax=0.9, tol=1e-4)
    for _ in 1:100
        tmid = 0.5*(tmin + tmax)
        err  = revenue_error(tmid, params, zgrid, P, Gtarget)
        if abs(err) < tol
            return tmid
        elseif err > 0
            # revenue too high => reduce τ
            tmax = tmid
        else
            # revenue too low => increase τ
            tmin = tmid
        end
    end
    return 0.5*(tmin + tmax)
end

function gini_coefficient(x::Vector{Float64}, ω::Vector{Float64})
    idx = sortperm(x)
    x_sorted = x[idx]
    ω_sorted = ω[idx]

    total_x = dot(x_sorted, ω_sorted)
    if total_x <= 0
        return 0.0
    end

    cum_x = cumsum(x_sorted .* ω_sorted) ./ total_x
    # Weighted sum of cumulative shares
    B = 0.0
    for i in 1:length(x)
        B += ω_sorted[i]*cum_x[i]
    end
    return 1.0 - 2.0*B
end

function main()
  
    base_params = calibrate_baseline()

    zgrid, P = tauchen_discretize(base_params.ρ, base_params.σz, base_params.nZ)

    Kstar_base = find_equilibrium_K(base_params, zgrid, P)

    r_base = base_params.α * base_params.A * Kstar_base^(base_params.α - 1) - base_params.δ
    w_base = (1 - base_params.α) * base_params.A * Kstar_base^base_params.α

    Vb, pb_idx = solve_household(base_params, r_base, w_base, zgrid, P)
    distb = stationary_distribution(base_params, pb_idx, P)

    Y_base = base_params.A * Kstar_base^base_params.α * 1.0^(1 - base_params.α) # L=1
    rev_base = 0.0
    Agrid = asset_grid(base_params)
    for iz in 1:base_params.nZ
        for ia in 1:base_params.nA
            mass = distb[ia, iz]
            y = w_base*zgrid[iz]
            rev_base += tax_paid(y, w_base, base_params.τ, base_params.λ)*mass
        end
    end
    KoverY_base = Kstar_base / Y_base

    N = base_params.nA * base_params.nZ
    after_tax_vec = zeros(N)
    weight_vec    = zeros(N)
    idx = 1
    for iz in 1:base_params.nZ
        for ia in 1:base_params.nA
            mass = distb[ia, iz]
            inc  = post_tax_income(w_base*zgrid[iz], w_base, base_params.τ, base_params.λ)
            after_tax_vec[idx] = inc
            weight_vec[idx]    = mass
            idx += 1
        end
    end
    gini_aftertax_base = gini_coefficient(after_tax_vec, weight_vec)

    asset_vec = zeros(N)
    idx = 1
    for iz in 1:base_params.nZ
        for ia in 1:base_params.nA
            asset_vec[idx] = Agrid[ia]
            idx += 1
        end
    end
    gini_assets_base = gini_coefficient(asset_vec, weight_vec)

    println("=== Baseline (λ=0) ===")
    println("r        = ", r_base)
    println("w        = ", w_base)
    println("K        = ", Kstar_base)
    println("Y        = ", Y_base)
    println("K/Y      = ", KoverY_base)
    println("τ (base) = ", base_params.τ)
    println("Rev/Y    = ", rev_base / Y_base)
    println("Gini(After-tax Labor) = ", gini_aftertax_base)
    println("Gini(Assets)          = ", gini_assets_base)
    println("--------------------------------------------------")

    
    reform_params = deepcopy(base_params)
    reform_params.λ = 0.15

    τ_star = find_tau_for_revenue(reform_params, zgrid, P, rev_base; tmin=0.0, tmax=0.9, tol=1e-4)
    reform_params.τ = τ_star

    Kstar_ref = find_equilibrium_K(reform_params, zgrid, P)
    r_ref = reform_params.α * reform_params.A * Kstar_ref^(reform_params.α - 1) - reform_params.δ
    w_ref = (1 - reform_params.α) * reform_params.A * Kstar_ref^reform_params.α

    Vref, pref_idx = solve_household(reform_params, r_ref, w_ref, zgrid, P)
    dist_ref = stationary_distribution(reform_params, pref_idx, P)

    Y_ref = reform_params.A * Kstar_ref^reform_params.α
    rev_ref = 0.0
    for iz in 1:reform_params.nZ
        for ia in 1:reform_params.nA
            mass = dist_ref[ia, iz]
            y = w_ref*zgrid[iz]
            rev_ref += tax_paid(y, w_ref, reform_params.τ, reform_params.λ)*mass
        end
    end
    KoverY_ref = Kstar_ref / Y_ref

    after_tax_vec_ref = similar(after_tax_vec)
    weight_vec_ref    = similar(weight_vec)
    idx = 1
    for iz in 1:reform_params.nZ
        for ia in 1:reform_params.nA
            mass = dist_ref[ia, iz]
            inc  = post_tax_income(w_ref*zgrid[iz], w_ref, reform_params.τ, reform_params.λ)
            after_tax_vec_ref[idx] = inc
            weight_vec_ref[idx]    = mass
            idx += 1
        end
    end
    gini_aftertax_ref = gini_coefficient(after_tax_vec_ref, weight_vec_ref)

    asset_vec_ref = similar(asset_vec)
    idx = 1
    for iz in 1:reform_params.nZ
        for ia in 1:reform_params.nA
            asset_vec_ref[idx] = Agrid[ia]
            idx += 1
        end
    end
    gini_assets_ref = gini_coefficient(asset_vec_ref, weight_vec_ref)

    println("=== Reform (λ=0.15) ===")
    println("Found τ = ", τ_star)
    println("r        = ", r_ref)
    println("w        = ", w_ref)
    println("K        = ", Kstar_ref)
    println("Y        = ", Y_ref)
    println("K/Y      = ", KoverY_ref)
    println("Rev/Y    = ", rev_ref / Y_ref)
    println("Gini(After-tax Labor) = ", gini_aftertax_ref)
    println("Gini(Assets)          = ", gini_assets_ref)

end

function lorenz_curve(x::Vector{Float64}, ω::Vector{Float64})
    idx = sortperm(x)
    x_sorted = x[idx]
    ω_sorted = ω[idx]
    cum_pop = cumsum(ω_sorted)
    total_x = dot(x_sorted, ω_sorted)
    cum_x = cumsum(x_sorted .* ω_sorted) ./ total_x
    return cum_pop, cum_x
end
