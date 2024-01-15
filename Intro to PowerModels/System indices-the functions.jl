function loads(buses)
    L=[]
    for b in buses
        if b[2]["bus_type"]!=2 # Cette décision implique que le slack bus fait partie des loads
            push!(L,b)
        end
    end
    return L
end

function generators(buses)
    G=[]
    for b in buses
        if b[2]["bus_type"]==2
            push!(G,b)
        end
    end
    return G
end


function Lindex(mydata, i)

    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])

    if buses[i][2]["bus_type"] == 2
        println("Im very sorry to announce you that node $i is a generator. I thereby cannot compute its L index.")
    end
    
    F = Fmatrix(mydata)
    j = 0 # Index du bus i dans le vecteur loads, ie ligne de la Fmatrix qui lui correspond
    for b in loads(buses)
        j += 1
        if b[2]["bus_i"] == i
            break
        end
    end


    sum = 0
    for (k,b) in enumerate(generators(buses))
            sum += F[j,k]*b[2]["vm"]
    end

    return abs(1-sum/buses[i][2]["vm"])

end

function global_Lindex(mydata)
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    Li = []
    L_subscripts = []

    for (i,b) in enumerate(buses)
        if b[2]["bus_type"] == 1
            push!(Li,Lindex(mydata,i))
            push!(L_subscripts,i)
        end
    end
    k = argmax(Li)

    println("L index of the system (at equilibrium) : $(Li[k])")
    println("Most critical load (wrt L-index): $(L_subscripts[k])")

    for (j,l) in enumerate(Li)
        if l>0.2
            println("Bus $(L_subscripts[j]) can be considered critical, as its L-index equals $l")
        end
    end
    #histogram(Li,show=true,bins=24,label="L_index distribution of the nodes")
    return Li
end


function SVD_index(mydata)
    J = Jacobian(mydata)
    N = size(J,1)
    dec = svd(J)
    sigma = dec.S[N]
    println("Lowest singular value of the system : $sigma")
    u = dec.U
    v = dec.V
    t_max = 0
    t_index = 0
    k = convert(Int64,N/2)
    if k != N/2
        println("How could the Jacobian be of dimension $N x $N ? An odd number of rows and columns is weird")
    end

    sensitivity_vector = []
    x_vector = []
    for i = (k+1):N

        if mydata["bus"]["$(i-k)"]["bus_type"] == 1
            push!(sensitivity_vector,abs(v[i,N]))
            push!(x_vector,(i-k))
        end    
        if abs(v[i,N])>t_max
            t_max = abs(v[i,N])
            t_index = i
        end
    end
    t_index -= k
    println("Weakest PQ bus (wrt SVD method) : $(t_index)") 

    return sensitivity_vector
end

function mult_SVD_index(mydata) # En chantier

    J = Jacobian(mydata)
    N = size(J,1)
    dec = svd(J)
    S = dec.S
    v = dec.V
end

function reduced_SVD_index(mydata)
    J = reduced_Jacobian(mydata)

    N = size(J,1)
    dec = svd(J)
    sigma = dec.S[N]
    println("Lowest singular value of the reduced system : $sigma")
    u = dec.U
    v = dec.V
    t_max = 0
    t_index = 0

    sensitivity_vector = []
    x_vector = []
    for i = 1:N

        if mydata["bus"]["$i"]["bus_type"] == 1
            push!(sensitivity_vector,abs(v[i,N]))
            push!(x_vector,i)
        end    
        if abs(v[i,N]) > t_max
            t_max = abs(v[i,N])
            t_index = i
        end
    end

    println("Weakest PQ bus (wrt reduced SVD method) : $(t_index)") 
    return sensitivity_vector
end

function compare_two_indices(t1,t2,name1 = "First index",name2 = "Second index") #t2 will be rescaled, so if the scale of an index matters more, let it be t1
    N = length(t1)
    x_vector = 1:N
    c = mean(t1)/mean(t2)
    t2 .= t2*c
    t = round(cor(t1,t2), digits = 3)
    plot(x_vector, t1, seriestype=:scatter, color=:red, label="")
    plot!(x_vector, t1, seriestype=:line, label = name1, show=true, color=:red)
    plot!(x_vector, t2, seriestype=:scatter, color=:green, label="")
    plot!(x_vector, t2, seriestype=:line, label = name2*"(x$(round(c, digits =2)))", show=true, color=:green)
    plot!([0.9], [0], seriestype=:line, label = "Pearson correlation : $t", color=:white)
    
    readline()
end