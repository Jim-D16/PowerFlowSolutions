using PowerModels, IntervalArithmetic, IntervalRootFinding
using LinearAlgebra, StaticArrays, ForwardDiff
using JuMP
using Ipopt
using DelimitedFiles
using ProgressBars
using Plots
using PrettyTables
using Statistics

include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Data modifier-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Visualisation-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Graph properties-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/System indices-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Very generic functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Matrices-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/PFequations builder.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Solutions and losses-the functions.jl")


function main()
gr()
solver = Ipopt.Optimizer
stopper=0


# bus type (1 = PQ, 2 = PV (jaune), 3 = slack, 4 = isolated)

# ------------------------------------------Initialisation------------------------------------------


global mydata = make_basic_network(PowerModels.parse_file("C:/Users/jim.delitroz/.julia/packages/PowerModels/4b7PA/test/data/matpower/case24.m"))
#mydatab = make_basic_network(PowerModels.parse_file("C:/Users/jim.delitroz/.julia/packages/PowerModels/4b7PA/test/data/matpower/case24b.m"))
#mydata = make_basic_network(PowerModels.parse_file("C:/Users/jim.delitroz/Documents/pantagruel.json"))
#mydata = make_basic_network(PowerModels.parse_file("C:/Users/jim.delitroz/.julia/packages/PowerModels/4b7PA/test/data/matpower/case5_gap.m"))







n = length(mydata["bus"])
B = Susceptance(mydata)
#B = round.(B, digits = 1)

G = Conductance(mydata) # Non-diago terms of G are usually <0
#G = round.(G, digits = 1)





buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
#=
branches = sort(collect(mydata["branch"]), by = x -> parse(Int64,x[1]))
gens = sort(collect(mydata["gen"]), by = x -> x[2]["gen_bus"])
=#

power_losses_via_pf(mydata)


NetworkDrawing(mydata)



# ------------------------------------------Modification des données------------------------------------------
#=

T = vector_angles(7,1) 

init_voltage_angles!(mydata,T)
=#
# bus_type!(mydata,2,2) 
# bus_type!(mydata,4,2)
                
# randomize_init!(mydata)

#ShowAdjacency(mydata)

ShowMe(mydata["bus"])
#NetworkDrawing(mydata)

# ------------------------------------------Résolution des équations------------------------------------------
result = solve_ac_pf(mydata,solver)
Result_buses = nice_solution(mydata,result)

for b in Result_buses
    println(b[2]["va"])
end
for b in Result_buses
    println(b[2]["vm"])
end

ShowMe(Result_buses)

# ------------------------------------------Etude de point d'équilibre-----------------------------------------

finaldata = deepcopy(mydata)
finaldata["bus"] = Result_buses # Attention ! Ils ne sont pas ordonnés. Si je les ordonne, alors ils deviennent un Vector(Pair(String,Any)) et du coup, plus possible de calculer la jacobienne

#=L_sens = global_Lindex(finaldata)
svd_sens = SVD_index(finaldata)
red_sens = reduced_SVD_index(finaldata)

compare_two_indices(svd_sens, red_sens, "SVD sensitivity", "Reduced-SVD sensitivity")
#savefig("Case 24 svd vs reduced svd")
=#
#NetworkDrawing(finaldata["bus"],finaldata["branch"])













end

main()
