using PowerModels, LinearAlgebra, StaticArrays, IntervalArithmetic, IntervalRootFinding, JuMP, Ipopt, ProgressBars, ForwardDiff, Dates

include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/Data modifier-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/Visualisation-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/Graph properties-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/System indices-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/Very generic functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/PFequations builder.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/Matrices-the functions.jl")


global mydata = make_basic_network(PowerModels.parse_file("C:/Users/jim.delitroz/.julia/packages/PowerModels/4b7PA/test/data/matpower/case5_gap.m"))

function main(k)
solver = Ipopt.Optimizer
stopper = 0

# ------------------------------------------Initialisation------------------------------------------


#mydatab = make_basic_network(PowerModels.parse_file("C:/Users/jim.delitroz/.julia/packages/PowerModels/4b7PA/test/data/matpower/case24b.m"))
n = length(mydata["bus"])

buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])

D = starting_interval_KS(mydata, pi)

m = mid.(D)


# F = (x -> SVector{2n}([f(x) for f in return_f(mydata)]))
F = (x -> SVector{n}([f(x) for f in return_f_KS(mydata, k)]))
#J = (x -> ForwardDiff.jacobian(F, x))

rts = roots(F, D, Krawczyk)

#k = round(k, digits = 3)

if length(rts) == 0
    #println("Ratio $k lead to 0 sol")
else
    for r in rts
        m = mid(r.interval)
        if length(rts) > 2
            #plot!([m[1]], [m[2]], seriestype = :scatter, color = :green, label="", show = true)
        else
            #plot!([m[1]], [m[2]], seriestype = :scatter, color = :red, label="", show = true)
        end
    end
    if length(rts) == 1
        println("Critical value !? $k") # Closest candidat for the "jim standard 3-bus test case" k = 0.34860. Note that k = 0.34861 yields 0 solution
        println("Solution : $m")
        println("")
    end
end

return length(rts)

#K = Krawczyk(F, J)

# ------------------------------------------Recherche itérative de solutions------------------------------------------


#=
mysearch = DepthFirstSearch(D, K, 1e-2)
counter = 0
start_time = now()

println("I'm searching")


for tree in (mysearch) # The first dimension along which the subdivision happens is the 27th. (in case24.m + DepthFirst)

    if counter % 5 == 0
        for l in tree.leaves
            print("$l  ")
            println(round(maximum(diam.(l.second.data.interval)), digits = 2))
            println("")
        end
    end
      
    global final_tree = tree
    elapsed_time = now() - start_time
    start_time = now()
    #println("Iter $(counter += 1), in $elapsed_time")

end

println("Number of iters : $counter")

println("Final result : ")
solutions = []
for l in final_tree.leaves
    println(l.second.data)

    m = mid.(l.second.data.interval)
    cand = F(m) # Affiche l'évaluation de F au milieu de notre petit intervalle, pour indiquer si l'intervalle est proche de contenir une solution

    println(cand)
    if maximum(cand) < 5e-2 # Si le milieu de l'intervalle est proche d'une solution, on le garde et on va le Show. Cependant, avec une tol de 1e-2 sur case5_gap, l'intervalle qui contient C ne satisfait pas cand < 5e-2
        push!(solutions, m)
    end


    println("")
end
if stopper == 0
    println("Bizarre") # Signifierait que C n'est dans aucune subdivision finale, or C est une solution
end
=#

end

function find_critical_k_by_bisection()
    upper = 3
    lower = 0
    k = 0.3
    finalk = 0.3

    for i = 1:10
        println(k)
        Z = main(k)
        if Z > 0
            lower = k
            k = (k + upper)/2 
        else
            upper = k
            k = (k+lower)/2
        end
        finalk = k
    end
    println("im done")
    println("Final k = $finalk")
    return finalk
end
find_critical_k_by_bisection()
    

