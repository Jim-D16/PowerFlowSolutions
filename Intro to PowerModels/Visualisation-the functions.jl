using Plots
using Graphs
using GraphRecipes

function slack(buses)
    i = 0
    y = 0
    for b in buses
        if b[2]["bus_type"] == 3
            i+=1
            y=b[2]["bus_i"]
        end
    end
    if i == 0
        println("We do not have a slack bus in this collection")
    elseif i>1
        println("We have $i slack buses in this collection")
    end
    return y
end

function circle(x0,y0,r)
    theta = LinRange(0,2*pi,500)
    y = y0.+ r*cos.(theta)
    x = x0.+ r*sin.(theta)
    return(x,y)
end

function ShowMeRoots(mydata, rts, colorvect = [:red, :blue, :yellow, :purple])
    superpose = false
    n = length(mydata["bus"])
    for (j,r) in enumerate(rts)

        m = mid.(r.interval)
    
        temp_IntervalBuses = deepcopy(mydata["bus"])
        IntervalBuses = sort(collect(temp_IntervalBuses), by = x -> x[2]["bus_i"])
    
        for (i,b) in enumerate(IntervalBuses)
            b[2]["va"] = m[i]
            b[2]["vm"] = m[i + n]
        end
        ShowMe(IntervalBuses, colorvect[j], false, superpose)
        superpose = true
    end
    readline()
end

function ShowMeLeaves(mydata,tree, colorvect = [:red, :blue, :yellow, :purple, :pink, :lime, :orange, :lavender, :khaki, :salmon2])
    superpose = false
    n = length(mydata["bus"])
    for (j,l) in enumerate(tree.leaves)

        m = mid.(l.second.data.interval)
    
        temp_IntervalBuses = deepcopy(mydata["bus"])
        IntervalBuses = sort(collect(temp_IntervalBuses), by = x -> x[2]["bus_i"])
    
        for (i,b) in enumerate(IntervalBuses)
            b[2]["va"] = m[i]
            b[2]["vm"] = m[i + n]
        end
        ShowMe(IntervalBuses, colorvect[j], false, superpose)
        superpose = true
    end
    readline()
end

function ShowMePoints(mydata, solutions, colorvect = [:red, :blue, :yellow, :purple, :pink, :lime, :orange, :lavender, :khaki, :salmon2])
    superpose = false
    n = length(mydata["bus"])

    for (j,m) in enumerate(solutions)
        temp_IntervalBuses = deepcopy(mydata["bus"])
        IntervalBuses = sort(collect(temp_IntervalBuses), by = x -> x[2]["bus_i"])
    
        for (i,b) in enumerate(IntervalBuses)
            b[2]["va"] = m[i]
            b[2]["vm"] = m[i + n]
        end
        ShowMe(IntervalBuses, colorvect[j % length(colorvect)+1], false, superpose)
        superpose = true
    end
    readline()
end


function ShowMe(buses, mycolor = :red, pause = true, superpose = false)
    slack_id = slack(buses)
    ang = []
    r = []
    slack_ang = []
    slack_r = []

    for b in buses
        if b[1]!="$slack_id"
            push!(ang,b[2]["va"])
            push!(r,b[2]["vm"])
        else
            push!(slack_ang,b[2]["va"])
            push!(slack_r,b[2]["vm"])
        end
    
    end 

    X = r.*cos.(ang)
    Y = r.*sin.(ang)
    slack_X = slack_r.*cos.(slack_ang)
    slack_Y = slack_r.*sin.(slack_ang)

    if !superpose
        plot(circle(0,0,1.15),seriestype=[:shape], show=true, fillalpha=0.2, aspect_ratio=:equal,legend=false)
        plot!(circle(0,0,0.85),seriestype=[:shape], show=true, fillalpha=1, c=:white,legend=false)
    end

    scatter!(X,Y,xlims=(-1.33,1.33),ylims=(-1.33,1.33),show=true,legend=false,markersize=2,markercolor=mycolor)
    println("The green bus is the slack bus")
    scatter!(slack_X,slack_Y,xlims=(-1.33,1.33),ylims=(-1.33,1.33),show=true,legend=false,markersize=3,markercolor=:green, markerstrokecolor=mycolor)

    if pause
        readline()
    end
end



function NetworkDrawing(buses,branches) 
    g=SimpleGraph(0)
    N=length(buses)
    vert_labels=[]
    vert_colors=[]
    for b in buses
        add_vertex!(g)
        push!(vert_labels,b[2]["bus_i"])
        if b[2]["bus_type"]==2
            push!(vert_colors,:yellow)
        else
            push!(vert_colors,:green)
        end
    end
    N=length(buses)
    for B in branches   
        if B[2]["f_bus"]<=N && B[2]["t_bus"]<= N
        add_edge!(g,B[2]["f_bus"],B[2]["t_bus"])
        end
    end   
    graphplot(g,show=true, curves=false, names=vert_labels, markersize=0.2, markercolor=vert_colors, fontsize=8, markershape=:circle)
    readline()
    return g
end

function NetworkDrawing(mydata)
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    branches=sort(collect(mydata["branch"]), by = x -> parse(Int64,x[1]))
    NetworkDrawing(buses,branches)
end


#=function NetworkDrawing(Y,Is_Y_A_Matrix)
    if Is_Y_A_Matrix != true
        println("I am not sure of what graph you wanted me to draw, but I will try my best")
    end
    g=SimpleGraph(0)
    N=size(Y,1)
    if size(Y,2)!=N
        println("Please provide a square admittance matrix")
        exit()
    end
    for i=1:N
       add_vertex!(g)
    end

    for i=2:N
        for j=1:(i-1)
            if Y[i,j]!=0
                add_edge!(g,i,j)
            end
        end
    end
    graphplot(g,show=true,curves=false,names=1:N,markersize=0.35,markercolor=:yellow, fontsize=3, markershape=:circle)
    readline()
end
=#

function ShowAdjacency(mydata)
    Y=Admittance(mydata)
    N=size(Y,1)
    for i=1:N
        for j=1:N
            print("$(1*(Y[i,j] != 0)) ")
        end
        println("")
    end
end


function area(buses)
    mang=0
    Mang=0

    for b in buses
        if (b[2]["va"] < mang)
            mang=b[2]["va"]
        end
        if (b[2]["va"] > Mang)
            Mang=b[2]["va"]
        end
    end

    Xgo_U=[0,1.15*cos(Mang)]
    Ygo_U=[0,1.15*sin(Mang)]
    Xgo_L=[0,1.15*cos(mang)]
    Ygo_L=[0,1.15*sin(mang)]
    plot!(Xgo_L,Ygo_L,legend=:false)
    plot!(Xgo_U,Ygo_U,legend=:false,show=true)

    readline()
end