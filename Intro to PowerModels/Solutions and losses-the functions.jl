function Powers(V, T, B, G = B.*0)
    P = []
    n = length(V)
    for i = 1:n
        Pi = 0
        for j = 1:n
            if (B[i,j] + G[i,j]) != 0 
                Pi += V[j]*(G[i,j]*cos(T[i]-T[j]) + B[i,j]*sin(T[i]-T[j]))
            end
        end
        Pi *= V[i]
        push!(P, Pi)
    end
    return P
end


function Reactive_Powers(V, T, B, G = B.*0)
    Q = []
    n = length(V)
    for i = 1:n
        Qi = 0
        for j = 1:n
            if (B[i,j] + G[i,j]) != 0 
                Qi += V[j]*(G[i,j]*sin(T[i]-T[j]) - B[i,j]*cos(T[i]-T[j]))
            end
        end
        Qi *= V[i]
        push!(Q, Qi)
    end
    return Q
end

function return_f_KS(V, B, P, slackbus_id = 1)
    sol = []
    t = slackbus_id - 1 
    t2 = slackbus_id + 1
    n = length(V)

    for i = 1:t
        eq = (x -> Powers(V, [0; collect(x)], B)[i] - P[i])
        push!(sol, eq)
    end

    for i = t2:n
        eq = (x -> Powers(V, [0; collect(x)], B)[i] - P[i])
        push!(sol, eq)
    end
    return sol
end

function random_B(scale = 12, dim = 3, showit = false)

    B = zeros(dim, dim)

    for i = 1:(dim-1)
        for j = (i+1):dim 
            b = scale*rand()
            B[i,j] = b
            B[j,i] = b
        end    
    end

    for i = 1:dim
        sum_of_line = 0
        for j = 1:dim
            sum_of_line += B[i,j]
        end
        B[i,i] = -sum_of_line
    end
    
    if showit
        display(B)
    end

    return B
end


function next_interval(root, r = 0.4) #Provided one (multi-dim) root, generates an interval of radius r centered at this root
    # The choice of r = 0.4 as a default value culd be argued against, it probably makes sense to perform tests to see if 0.4 is relevant. Shrinking r would accelerate the program, but may lead to the loss of some solutions
    # I tried r = 0.2 on a random 3-bus, and every solution died instantly, but they worked very fine with r = 0.4. This hirondelle makes the 0.4 printemps imo
    mybox = []
    m = mid(root.interval)
    for mc in m
        I = interval(mc - r, mc + r)
        push!(mybox, I)
    end

    return IntervalBox(mybox)
end

function main_yielding_a_contradiction() # 3 bus system that contradicts the assumption (losses increase -> solutions diminish)
    #The accelerated method (ie the use of next_interval) supposes that the assumption is true, in particular, using it here will not yield a contradiction (although there is one)
    t2 = 6.06785081987959
    t3 = 4.028180014959005
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    B = [-7.53 4.33 3.2; 4.33 -5.49 1.16; 3.2 1.16 -4.36]
    println("B susceptance matrix")
    display(B)
    n = length(T)

    G0 = [-8.07466   4.95559   3.11907; # Contradicts the monotony of number of solutions with k : k=0 admits 2 sol; 0.2<=k<= admits 4
    4.95559  -6.77213   1.81654;
    3.11907   1.81654  -4.93561]   
    
    println("G0 loss matrix")
    display(G0)

    D1 = (-pi..pi)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) 

    L = collect(LinRange(0, 0.3, 14))
    for k in (L)
        println("k = $k")
        perform(V, T, B, G0, D, k)
    end
    println("im done")
    readline()
end