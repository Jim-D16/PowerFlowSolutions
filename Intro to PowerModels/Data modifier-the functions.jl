using Graphs
using GraphRecipes


function available_index(data_dict)
    t = 0
    for B in data_dict
        if parse(Int,B[1])>t
            t = parse(Int,B[1])
        end 
    end
    return t+1
end


function add_bus!(mydata,neighbours = [], bus_type = 1, angle = 0., voltage = 1.)
    N = available_index(mydata["bus"])
    push!(mydata["bus"], ("$(N)" => Dict{String, Any}("zone" => 1, "bus_i" => N, "bus_type" => bus_type, "vmax" => 1.1, "source_id" => Any["bus", N], "area" => 1, "vmin" => 0.9, "index" => N, "va" => angle, "vm" => voltage, "base_kv" => 230.0)))
    for t in neighbours
        add_branch!(mydata,N,t)
    end
end

function bus_type!(mydata,n,t) # Changes the type of bus n to type t
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    if buses[n][2]["bus_i"] == n
        buses[n][2]["bus_type"] = t
    else
        println("We did not change the type of bus $n as we could not find it. $n-th position is occupied by bus $(buses[n][2]["bus_i"])")
    end
end

function disconnect_bus!(mydata,t)
    Z = []
    for k in eachindex(mydata["branch"])
        B=mydata["branch"][k]
        if (B[2]["f_bus"] == t) || (B[2]["t_bus"] == t)
            push!(Z,k)
            #println("Im removing edge $(B[2]["f_bus"]), $(B[2]["t_bus"])")
        end
    end
    println("degree of vertex $t is $(length(Z))")
    Z = sort!(Z)
    for z in Z
        deleteat!(mydata["branch"],z)
        Z .= Z.-1
    end
end

function add_gen!(mydata, n, p = 1., q = 1.) # Attach a generator that will generate P and Q to the bus n
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    b = buses[n]
    if b[2]["bus_i"] != n
        println("We did not attach a generator to bus $n as we could not find it. $n-th position is occupied by bus $(buses[n][2]["bus_i"])")
    else
        if b[2]["bus_type"] != 2
            println("We wanted to attach a generator to bus $n, so we changed its type from $(b[2]["bus_type"]) to 2.")
            b[2]["bus_type"] = 2
        end
    M = available_index(mydata["gen"])
    push!(mydata["gen"], "$M" => Dict{String, Any}("vg" => 1., "mbase" => 100.0, "source_id" => Any["gen", M], "pg" => p, "model" => 2, "shutdown" => 0.0, "startup" => 0.0, "index" => M, "cost" => [-4000.0, 0.0], "qg" => q, "qmax" => q+.5, "gen_status" => 1, "qmin" => q-.5, "gen_bus" => n, "pmax" => p+1, "pmin" => p-1, "ncost" => 2))
    end
end

function remove_gen!(mydata,n) #Removes the generator n
    for g in mydata["gen"]
        if g[1] == "$n"
            delete!(mydata["gen"],g[1])
        end
    end
end

function remove_gen!(mydata,n,Index_Of_Bus) #Removes all the generators attached to bus n
    if Index_Of_Bus == false
        remove_gen!(mydata,n)
    else
        for g in mydata["gen"]
            if g[2]["gen_bus"] == n
                delete!(mydata["gen"],g[1])
            end
        end
    end
end

function add_branch!(mydata, starting_bus, arrival_bus, resistance = 0.05, reactance = 0.1, add_susceptance = 0.0) # La susceptance finale est calculée comme la partie imaginaire de l'inverse de (resistance + im*reactance), à laquelle on ajoute add_susceptance
    M = available_index(mydata["branch"])
    push!(mydata["branch"], ("$(M)" => Dict{String, Any}("br_r" => resistance, "rate_a" => 0, "shift" => 0.0, "rate_b" => 0, "br_x" => reactance, "rate_c" => 0, "g_to" => 0.0, "g_fr" => 0.0, "source_id" => Any["branch", starting_bus], "b_fr" => add_susceptance/2, "f_bus" => starting_bus, "br_status" => 1, "t_bus" => arrival_bus, "b_to" => add_susceptance/2, "index" => M, "angmin" => -pi/3, "angmax" => pi/3, "transformer" => false, "tap" 
    => 1.0)))
end

function remove_branch!(mydata,starting_bus,arrival_bus)

    for B in mydata["branch"]
        #println("Observing the branch $(B[2]["f_bus"]), $(B[2]["t_bus"])")
        if (B[2]["f_bus"] == starting_bus && B[2]["t_bus"] == arrival_bus) || (B[2]["t_bus"] == starting_bus && B[2]["f_bus"] == arrival_bus)         
            delete!(mydata["branch"],B[1])
            break
        end
    end
            
end

function nothing_but_a_cycle!(mydata, K) # Replaces all the data by a cycle of length K
    clear!(mydata)
    add_a_cycle!(mydata,K)
    bus_type!(mydata,1,3) # Set bus 1 as the slack bus
end

function clear!(mydata)
    empty!(mydata["bus"])
    empty!(mydata["branch"])
    empty!(mydata["gen"])
    empty!(mydata["load"])
end

function add_a_cycle!(mydata,K)
    N = length(mydata["bus"])
    add_bus!(mydata)
    for i = 1:K-1
        add_bus!(mydata,[N+i])
    end
    add_branch!(mydata,N+1,N+K)
end

function init_voltage_magnitudes!(mydata,V)

    for b in mydata["bus"]
        k = (b[2]["bus_i"])
        if k <= length(V)
            b[2]["vm"] = V[k]
        else
            println("We did not change the voltage magnitude of bus $k as only $(length(V)) angles were specified.")
        end
    end
end

function init_voltage_angles!(mydata,V)

    for b in mydata["bus"]
        k = (b[2]["bus_i"])
        if k <= length(V)
            b[2]["va"] = V[k]
        else
            println("We did not change the phase angle of bus $k as only $(length(V)) angles were specified.")
        end
    end
end

function vector_angles(L, wn = 1, noise = 0) # Length of the cycle, winding numbers, max magnitude of the noise (so that each angle diff equals 2pi*wn/L +- noise*rand())
    V = []
    for i = 1:(L-1)
        theta_i = wn*2*pi/L*(i-1) + noise*(rand()-0.5)
        push!(V, theta_i)
    end
    last_angle = (L-1)*pi*wn
    for v in V
        last_angle -= v
    end
    push!(V, last_angle)
    return V
end


function randomize_init!(mydata) #Randomiser les initialisations de toutes les variables, ie l'angle partout sauf au slack, et les magnitudes à chaque PQ bus
    for b in mydata["bus"]
        if b[2]["bus_type"] != 3
            b[2]["va"] = 2*pi*rand()
            if b[2]["bus_type"] == 1
                b[2]["vm"] += (rand()-0.5)/5
            end
        end
    end
end

function my_cyclic_susceptance_matrix(M, random = false) # Returns a matrix that describes a cycle. If random == true, it will instantiate a new matrix every time, otherwise, it is always the same

    if M > 20
        println("Sorry, my_cyclic_susceptance_matrix cannot provide a >20x20 matrix")
        return 0
    end
    B = zeros(M, M)
    if random
        for i = 1:(M-1)
            b = rand()*8
            B[i, i+1] = b
            B[i+1, i] = b
        end
        B[1, M] = rand()*8
        B[M, 1] = B[1, M]

        for i = 1:M
            sum_i = 0
            for j = 1:M
                sum_i -= B[i, j]
            end
            B[i, i] = sum_i
        end
    else
        arbitrary_values = [2.1, 4.8, 1.1, 5.2, 2.0, 0.3, 0.5, 7.0, 3.3, 6.7, 2.8, 2.8, 7.2, 0.2, 3.2, 0.9, 7.2, 2.0, 1.0, 0.8]
        for i = 1:(M-1)
            b = arbitrary_values[i]
            B[i, i+1] = b
            B[i+1, i] = b
        end
        B[1, M] = arbitrary_values[M]
        B[M, 1] = B[1, M]

        for i = 1:M
            sum_i = 0
            for j = 1:M
                sum_i -= B[i, j]
            end
            B[i, i] = sum_i
        end
    end
    return B
end


