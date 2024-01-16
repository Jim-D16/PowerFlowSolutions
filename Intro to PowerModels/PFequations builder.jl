function return_f(mydata)
    n = length(mydata["bus"])
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    gens = sort(collect(mydata["gen"]), by = x -> x[2]["gen_bus"])
    sol = []
    for i = 1:n
        if buses[i][2]["bus_type"] != 3
            push!(sol, (x -> equaPF(i,x,mydata)))
        else
            push!(sol, (x -> angconstr(i,x,mydata)))
        end
    end
    for i = 1:n
        if buses[i][2]["bus_type"] == 1
            push!(sol, (x -> equaPF(i,x,mydata, true))) 
        else
            push!(sol, (x -> voltconstr(i,x,mydata)))
        end
    end

    # Maintenant, vous pouvez appeler les fonctions stockées dans le vecteur sol avec la formulation suivante : 
    # result = [f(x) for f in return_f(mydata)] où x est un vecteur de taille 2n, qui contient tous les angles suivis de toutes les magnitudes

    return sol
end

function return_f_KS(mydata, k) # Une fct de R^(n-1)-> R^(n-1), qui construit le système d'éq. dont les theta_i (i != slack) sont les inconnues, et les n-1 éq. sont P_i pour i != slack. Formellement, return_f_KS rend quand même n équations, juste que celle du slack (j) c'est theta_j = theta_j_init
    
    n = length(mydata["bus"])
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    sol = []

    for i = 1:n
        if buses[i][2]["bus_type"] != 3
            push!(sol, (x -> equaKS(i, x, mydata, k)))
        else
            push!(sol, (x -> angconstr(i, x, mydata)))
        end
    end

    return sol
end


function voltconstr(i, T, mydata, normalized = false)
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    n = length(buses)

    if normalized
        return T[i+n] - 1.
    else
        return T[i+n] - buses[i][2]["vm"]
    end
end

function angconstr(i, T, mydata)
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    n = length(buses)
    return T[i] - buses[i][2]["va"]
end




function equaPF(i, T, mydata, reactive = false) # T est de longueur 2n, les n premiers termes sont les angles, les n suivants sont les magnitudes
    
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    gens = sort(collect(mydata["gen"]), by = x -> x[2]["gen_bus"])
    loads = sort(collect(mydata["load"]), by = x -> x[2]["load_bus"])
    n = length(buses)

    Y = Admittance(mydata)
    B = imag(Y)
    G = real(Y)

    Qi = 0
    Pi = 0
    for g in gens
        if g[2]["gen_bus"] == i
            Pi += g[2]["pg"]
            Qi += g[2]["qg"]
        elseif g[2]["gen_bus"] > i
            break
        end
    end

    for l in loads
        if l[2]["load_bus"] == i
            Pi -= l[2]["pd"]
            Qi -= l[2]["qd"]
        elseif l[2]["load_bus"] > i
            break
        end
    end

    sum = 0
    if !reactive
        for j = 1:n
            if Y[i,j] != 0 
                sum += T[n+j]*(G[i,j]*cos(T[i]-T[j]) + B[i,j]*sin(T[i]-T[j]))
            end
        end
        sum *= T[n+i]
        return Pi-sum
    else
        for j = 1:n
            if Y[i,j] != 0 
                sum += T[n+j]*(G[i,j]*sin(T[i]-T[j]) - B[i,j]*cos(T[i]-T[j]))
            end
        end
        sum *= T[n+i]
        return Qi-sum
    end

end

function equaKS(i, T, mydata, k) # T est de longueur n, les termes sont les angles
    # equaKS retourne la i-ème équation du système KS, ie Pi - sum (Bsin+Gcos) 
    
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    gens = sort(collect(mydata["gen"]), by = x -> x[2]["gen_bus"])
    loads = sort(collect(mydata["load"]), by = x -> x[2]["load_bus"])
    n = length(buses)


    Y = Admittance(mydata)
    B = imag(Y)
    #G = real(Y)
    G = -B .* k


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

    sum = 0

    for j = 1:n
        if Y[i,j] != 0 
            sum += G[i,j]*cos(T[i]-T[j]) + B[i,j]*sin(T[i]-T[j])
        end
    end

    return Pi-sum

end

function starting_interval_PF(mydata, AngTol = pi, VoltTol = 5e-2)
    n = length(mydata["bus"])
    A = []
    B = []
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])

    for b in buses
        push!(A, b[2]["va"])
        push!(B, b[2]["vm"])
    end
    C = [A;B]
    
    
    intervals = []
    
    for (i,c) in enumerate(C) # Intervalles autour d'une équation qui est une constrainte : rayon 1e-6
        if i <= n
            if buses[i][2]["bus_type"] == 3
                push!(intervals, (c - 1e-6..c + 1e-6)) # Slack angle
            else
                push!(intervals, (c - AngTol..c + AngTol)) # PV and PQ angle
            end
        else
            if buses[i-n][2]["bus_type"] > 1
                push!(intervals, c - 1e-6..c + 1e-6) # PV and slack voltage
            else
                push!(intervals, (c - VoltTol..c + VoltTol)) # PQ voltage
            end
        end
    end
    
    D = IntervalBox(intervals)
    return D
    
end

function starting_interval_KS(mydata, AngTol = pi)
    n = length(mydata["bus"])
    C = []
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])

    for b in buses
        push!(C, b[2]["va"])
    end
    intervals = []
    
    for (i,c) in enumerate(C) # Intervalles autour d'une équation qui est une constrainte : rayon 1e-6
        if buses[i][2]["bus_type"] == 3
            push!(intervals, (c - 1e-6..c + 1e-6)) # Slack angle
        else
            push!(intervals, (c - AngTol..c + AngTol)) # PV and PQ angle
        end
    end
    
    D = IntervalBox(intervals)
    return D
    
end