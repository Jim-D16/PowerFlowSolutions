using Plots
using PowerModels
using ProgressBars



function main()

mydata=PowerModels.parse_file("C:/Users/jim.delitroz/Documents/pantagruel.json")

buses = sort(collect(mydata["bus"]), by = x -> x.second["bus_i"])
branches=mydata["branch"]


X=[]
Y=[]

for b in (buses)
    push!(X,b[2]["coord"][1])
    push!(Y,b[2]["coord"][2])   
end

scatter!(X,Y,show=true,legend=false,markersize=0.8)
readline()








end
main()