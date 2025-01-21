TASK 1
function solve_basil_problem(N::Int, X::Float64, c::Float64, q::Float64,
                             p_min::Float64, p_max::Float64; p_step=1.0)
 
    price_grid = collect(p_min:p_step:p_max)
    M = length(price_grid)

   
    v = zeros(Float64, N+1)

    σ_approach = zeros(Int, N+1)

    σ_buy = zeros(Int, N+1, M)

    
    v[N+1] = -N*c  
    σ_approach[N+1] = 0  
  
    for n in reverse(0:N-1)
        vT = -n*c

       
        vA = -c  
        buy_reject_sum = 0.0
        for (pidx, p) in enumerate(price_grid)
            buy_payoff = X - p - (n+1)*c
            if buy_payoff >= v[n+2]
                buy_reject_sum += buy_payoff
                σ_buy[n+1, pidx] = 1
            else
                buy_reject_sum += v[n+2]
                σ_buy[n+1, pidx] = 0
            end
        end
        buy_reject_avg = buy_reject_sum / M

        vA += q * buy_reject_avg + (1-q)*v[n+2]

        if vA >= vT
            v[n+1] = vA
            σ_approach[n+1] = 1
        else
            v[n+1] = vT
            σ_approach[n+1] = 0
            σ_buy[n+1, :] .= 0
        end
    end

    return v, σ_approach, σ_buy, price_grid
end

TASK 2
function solve_job_search_with_separations()
    wages::Vector{Float64}, 
    pi_w::Vector{Float64},   
    beta::Float64,          
    p::Float64,             
    c::Float64;           
    tol::Float64 = 1e-8,    
    maxiter::Int = 10_000   
)
    
    M = length(wages)
    @assert length(pi_w) == M
    @assert abs(sum(pi_w) - 1.0) < 1e-12  


    VE = zeros(Float64, M)   
    VU = zeros(Float64, M)  
    policy = zeros(Int, M)   

    for iter in 1:maxiter
     
        VE_old = copy(VE)
        VU_old = copy(VU)

        EVU = dot(VU, pi_w)  

    
        denom = 1.0 - beta*(1-p)
        for i in 1:M
            VE[i] = ( wages[i] + beta*p*EVU ) / denom
        end

       
        for i in 1:M
            outside = c + beta*EVU
            if VE[i] >= outside
                VU[i] = VE[i]
                policy[i] = 1  
            else
                VU[i] = outside
                policy[i] = 0
            end
        end

    
        diffVE = maximum(abs.(VE .- VE_old))
        diffVU = maximum(abs.(VU .- VU_old))
        if max(diffVE, diffVU) < tol
            println("Converged in iteration $iter")
            break
        end
        if iter == maxiter
            println("Warning: Reached max iterations without full convergence.")
        end
    end

    return VU, VE, policy
end

using Printf

wages = [10.0, 15.0, 20.0, 25.0]     
pi_w  = [0.2,   0.3,   0.3,   0.2]   
beta  = 0.95
p     = 0.05  
c     = 1.0   

VU, VE, policy = solve_job_search_with_separations(wages, pi_w, beta, p, c)

println("Unemployed values: ", VU)
println("Employed values:   ", VE)
println("Accept/reject policy (1=accept, 0=reject): ", policy)


res_index = findfirst(x->x==1, policy)
if res_index == nothing
    println("No wage is high enough to accept!")
else
    println("Reservation wage = ", wages[res_index])
end

TASK 4
function solve_model_with_two_states()

    z_states = [:z1, :z2, :z3]
    P = [
        0.5  0.3  0.2
        0.2  0.7  0.1
        0.3  0.3  0.4
    ]


    x_states = [0,1,2,3,4,5]


    sigma(x,z) = begin
        if z == :z1
            return 0
        elseif z == :z2
            return x
        else
            return x < 5 ? x+1 : 3
        end
    end

    M, state_to_ind = build_joint_transition_matrix(P, sigma, x_states, z_states)


    pi_joint = stationary_distribution(M)

    pX = marginal_X(pi_joint, x_states, z_states, state_to_ind)

 
    E_X = expectation_X(x_states, pX)

    return (M, pi_joint, pX, E_X)
end

(M, pi_joint, pX, E_X) = solve_model_with_two_states()

println("Joint transition matrix M is size ", size(M))
println("Stationary distribution over (X,Z): ", pi_joint)
println("Marginal distribution over X:       ", pX)
println("Expected value of X:               ", E_X)
