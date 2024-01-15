
function PowerFlowUnknowns()
    T = []
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])

    for b in buses
        if b[2]["bus_type"] != 3
            push!(T,b[2]["va"])
            println("Theta_$(b[2]["bus_i"]) = $(b[2]["va"])")
        end
    end
    for b in loads(buses)
        if b[2]["bus_type"] != 3
            push!(T,b[2]["vm"])
            println("V_$(b[2]["bus_i"]) = $(b[2]["vm"])")
        end
    end
    println("Number of equations : $(length(T))")
    return T
end

function PowerFlowSystem(T) # A R^n->R^n fct that has to be zero for the power flow to be at equilibrium. 
    G = real(Admittance(mydata))
    B = Susceptance(mydata)
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    P=[]
    Q=[]
    theta=[]
    V=[]
    for b in buses
        push!(V,b[2]["vm"])
        push!(theta,b[2]["va"])
        push!(P,4)
    end

end




function nice_solution(mydata,result)

    Result_buses = deepcopy(mydata["bus"])
    Solution = result["solution"]["bus"]

    k=1
    for b in sort(collect(Result_buses), by = x -> x[2]["bus_i"])
        b[2]["va"] = Solution["$k"]["va"]
        b[2]["vm"] = Solution["$k"]["vm"]
        k+=1
    end
    return Result_buses
end

function power_losses(mydata) # Takes a state of the system as input and computes the active power losses. That is, it calculates how much power has been dissipated (ie generated - consumed) and divides it by the generated power, in order to output the percentile of generated power that ahs not been consumed. (bc it has been dissipated instead)
    
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    gens = sort(collect(mydata["gen"]), by = x -> x[2]["gen_bus"])
    loads = sort(collect(mydata["load"]), by = x -> x[2]["load_bus"])
    
    generated_power = 0
    consumed_power = 0

    P = []

    for i in eachindex(buses)
        Pi = 0
        for g in gens
            if g[2]["gen_bus"] == i
                Pi += g[2]["pg"]
            elseif g[2]["gen_bus"] > i
                break
            end
        end

        for l in loads
            if l[2]["load_bus"] == i
                Pi -= l[2]["pd"]
            elseif l[2]["load_bus"] > i
                break
            end
        end
        push!(P, Pi)
    end

    for p in P
        if p > 0
            generated_power += p
        else
            consumed_power -= p
        end
    end

    absolute_loss = generated_power - consumed_power
    println("Generated power that went in the network : $generated_power")
    println("Consumed power : $consumed_power")
    rel_loss = round(absolute_loss / generated_power, digits = 3)
    println("The total power loss of the system amounts for $(100*rel_loss)% of the generated power.")
    return rel_loss
end

function power_losses_via_pf(mydata) # If te power flow equations are satisfied, this is exactly equivalent to power_losses. It reaches the result by replacing Pi (defined as gen-load) by its expression in terms of voltages in the power flow equations

    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    n = length(buses)

    Y = Admittance(mydata)
    B = imag(Y)
    G = real(Y)
    P = []

    for i = 1:n
        Pi = 0
        for j = 1:n
            if Y[i,j] != 0 
                Pi += buses[j][2]["vm"]*(G[i,j]*cos(buses[i][2]["va"]-buses[j][2]["va"]) + B[i,j]*sin(buses[i][2]["va"]-buses[j][2]["va"]))
            end
        end
        Pi *= buses[i][2]["vm"]
        push!(P,Pi)
    end

    generated_power = 0
    consumed_power = 0

    for p in P
        if p > 0
            generated_power += p
        else
            consumed_power -= p
        end
    end

    absolute_loss = generated_power - consumed_power
    println("Generated power that went in the network, says PF : $generated_power")
    println("Consumed power, says PF : $consumed_power")
    rel_loss = round(absolute_loss / generated_power, digits = 3)
    println("The total power loss of the system amounts for $(100*rel_loss)% of the generated power, says PF.")
    return rel_loss

end