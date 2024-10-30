using Plots
task 1

function odd_or_even(n::Int)
if n % 2 ==0
    println("even")
else
    println("odd")
end
end
odd_or_even(6)

task 2

function compare_three(a::Int,b::Int,c::Int)
    if a > 0 & b > 0 & c > 0
        println("all numbers are positive")
    elseif a==0 & b== 0 & c==0
        println("all numbers are 0")
    else 
        println("at least one number is not positive")
    end
end
compare_three(0,0,0)

task 3

function my_factorial(m::Int)
    g = 1
    for i in 1:m
        g = g*i
    end
    println(g)
end

my_factorial(5)

task 4

function count_positives(arr::Array{Int})
    count = 0
    for num in arr
        if num > 0
            count +=1
        end
    end
    println("total count of positive numbers is: $count")|
end

task 5

function plot_powers(n::Int)
    plot(plot_powers,-10:0.1:10)


end
TASK 6




function count_positives_broadcasting(arr)
    sum(arr .> 0)
end
count_positives_broadcasting([1,2,3,4,-2,-6,-5,0])


