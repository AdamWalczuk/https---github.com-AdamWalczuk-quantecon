using LinearAlgebra
using DataFrame

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

println("exact solution: x_exact")




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