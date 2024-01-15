function adjacent(b,b2,branches)
    sus=[]
    if(b==b2)
        return false
    else

        #println("Checking the adjacency of vertices $(b[2]["bus_i"]) and $(b2[2]["bus_i"])")

        t=0
        for B in branches
            #println("Observing the branch $(B[2]["f_bus"]), $(B[2]["t_bus"])")
            if (B[2]["f_bus"]==b[2]["bus_i"] && B[2]["t_bus"]==b2[2]["bus_i"]) || (B[2]["t_bus"]==b[2]["bus_i"] && B[2]["f_bus"]==b2[2]["bus_i"])
                t+=1
                sus=B[2]["br_x"]
                #println("We have a match !") ; println("")
                
            end
        end

        if t>1
            println("Error in the adjacency check, the branch $(b[2]["bus_i"]) and $(b2[2]["bus_i"]) seems to appear more than once")
            z=true
        elseif t==0
            z=false
        else
            z=true
        end
        if z==false
            return z
        else
            return (z,sus)
        end
    end
end

function clustering(buses,branches)

    C=0.
    println("Number of buses = $(length(buses))")
    for b in buses
        S=[]
        for b2 in buses
            if adjacent(b,b2,branches)[1]
                push!(S,b2[2]["bus_i"])
            end
        end
        k=length(S)
        t=0
        if k>1
            for B in branches
                if (B[2]["f_bus"] in S) && (B[2]["t_bus"] in S)
                    t+=1
                end
            end
            z=2*t/(k*(k-1))
            C+=z
            #println("Vertex $(b[2]["bus_i"]) : C=$z")
        end
    end
    C=C/length(buses)
    return C
end