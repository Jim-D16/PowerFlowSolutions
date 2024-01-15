# ------------------------------------------Tentative de représenter une partie du graphe de Pantagruel------------------------------------------

# Il faut inclure des packages etc.

SomeBuses_indices=[]
SomeBuses=[]
SomeBranches=[]
g=SimpleGraph(0)
i=1

for b in buses
    if b[2]["country"]=="CH"
        add_vertex!(g)
        push!(SomeBuses,b)
        push!(SomeBuses_indices,b[2]["bus_i"])
    end
end 

for B in (branches)   
    f=B[2]["f_bus"]
    t=B[2]["t_bus"]
    if in(f,SomeBuses_indices) && in(t,SomeBuses_indices)
        add_edge!(g,findfirst(isequal(f),SomeBuses_indices),findfirst(isequal(t),SomeBuses_indices))
    end
end

println(g)
graphplot(g,show=true,curves=false,names=SomeBuses_indices,markersize=0.2,markercolor=:yellow, fontsize=8, markershape=:circle)
readline()

