using LinearAlgebra, DataFrame, Roots, Optim 

problem 1
function iterative_solver(f, x0; α=0.0, ϵ=1e-6, maxiter=1000) 
    x = x0 
    x_values = [x] 
    residuals = [abs(f(x))] 
    for iter in 1:maxiter 
        g_x = f(x) + x 
        x_new = (1 - α) * g_x + α * x 
        push!(x_values, x_new) 
        push!(residuals, abs(f(x_new))) 
        if abs(x_new - x) / (1 + abs(x)) < ϵ 
            return (0, x_new, f(x_new), abs(x_new - x), x_values, residuals) 
        end 
        
        x = x_new 
    end 
    return (1, NaN, NaN, abs(x_new - x), x_values, residuals) 
end 


problem 2
α = 2.0
β = 1.0
function exact_solutions(α,β)
    x5=1
    x4=x5
    x3=x4
    x2=x3
    x1=x2 + β - (α - β)
    return [x1,x2,x3,x4,x5]
end

function solve(α,β)

A= [1 -1 0 α-β β;
    0 -1 1  0  0;
    0  0 1 -1  0;
    0  0 0  1 -1;
    0  0 0  0  1;]
b= [α, 0, 0, 0, 1]

x_num = A/b

x_exact = exact_solution(α,β)

res = b - A * x_num

relative_res = norm(res)/norm(b)

condition_number = cond(A)

x_exact, x_num, relative_res, condition_number =  solve(α,β)

##println("exact solution: x_exact")




return x_exact, x_num, relative_res, condition_number
end


function table()
    α = 0.1
    β_values = [1, 10 , 100, 1000, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12]

    results = DataFrame(α=Float64[], x1_exact=Float64[],x1_num=Float64[], res=Float64[], condition_number=Float64[])
    for β in β_values
        x_exact, x_num, res, condition_number = solve(α,β)pus!(results,(β, x_exact[1], x_num[1], residual, condition_number))
    end
return results
end


problem 3

function internal_rate(C::Vector)
    T=length(C) - 1
    NPV(r) = sum(C[t+1]/(1+r)^t for t in 0:T)
    r_internalrate = find_zero(NPV, 0.1)
    return r_internalrate
end

problem 4
 
 
    function production_function(α, σ, x1, x2) 
        if σ == 1 
            return x1^α * x2^(1 - α) 
        else  
            return (α * x1^((σ - 1) / σ) + (1 - α) * x2^((σ - 1) / σ))^(σ / (σ - 1))
         end 
    end 
          function plot_production_function(α, σ) 
            x1 = 0:0.1:10
            x2 = 0:0.1:10 
            z = [production_function(α, σ, x1_i, x2_j) for x1_i in x1, x2_j in x2]
             contour(x1, x2, z, xlabel="x1", ylabel="x2", title="CES Production Function (σ=$σ)") 
            end 
             cost_minimization(α, σ, w1, w2, y) 
              function cost_minimization(α, σ, w1, w2, y) 
                cost(x) = w1 * x[1] + w2 * x[2] 
                constraint(x) = production_function(α, σ, x[1], x[2]) - y 
                result = optimize(cost, [1.0, 1.0], lower_bounds=[0.0, 0.0], 
                f_inequality_constraints=constraint) 
                return Optim.minimizer(result)[1], Optim.minimizer(result)[2], Optim.minimum(result) 
            end 
    
            function plot_cost_and_demand(α, y, σ_values, w2) 
                w1_range = 0.1:0.1:5 
                for σ in σ_values 
                    costs = [] 
                    x1s = []
                    x2s = [] 
                    for w1 in w1_range 
                        x1, x2, cost = cost_minimization(α, σ, w1, w2, y) 
                        push!(costs, cost) 
                        push!(x1s, x1) 
                        push!(x2s, x2) 
                    end 
                    plot(w1_range, costs, label="Cost (σ=$σ)", xlabel="w1", ylabel="Value", title="Cost and Input Demands") 
                    plot!(w1_range, x1s, label="x1 (σ=$σ)") 
                    plot!(w1_range, x2s, label="x2 (σ=$σ)") 
                end 
            end 
            
        