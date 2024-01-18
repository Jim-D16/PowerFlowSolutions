using PowerModels, LinearAlgebra, StaticArrays, IntervalArithmetic, IntervalRootFinding, JuMP, Ipopt, ProgressBars
using Colors, ForwardDiff, Dates, Plots

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

function perform(V, T, B, G0, D, k = 0, i = 1)
    n = length(V)
    G = -G0 .* k
    
    P = Powers(V, T, B, G)

    #=
    generated_power = 0
    lost_power = 0


    for p in [p1, p2, p3]
        lost_power += p
        if p > 0
            generated_power += p
        end
    end

    if generated_power == 0
        rel_loss = 0
    else
        rel_loss = lost_power / generated_power
    end 
    println("Total loss _ $lost_power")
    println("Relative loss : $rel_loss") =#
    
    
    F = (x -> SVector{n-1}([f(x) for f in return_f_KS(V, B, P)]))

    rts = roots(F, D, Krawczyk, 1e-4)

    display(rts)
    
    THD = false
    if i > 15
        i = 15
    end
    mycolors = cgrad([:tomato2, :green], 15, categorical = true)
    if n == 3
        for r in rts
            m = mid(r.interval)
            if k > 0
                plot!([m[1]], [m[2]], marker = :circle, color = mycolors[i], label="", show = true)
            else
                plot!([m[1]], [m[2]], marker = :circle, markersize = 5, color = :red, label="", show = true)
            end

        end
    elseif n == 4
        for r in rts
            m = mid(r.interval)
            if k > 0
                plot3d!([m[1]], [m[2]], [m[3]], marker = :circle, color = mycolors[i], label="", show = true)
            else
                plot3d!([m[1]], [m[2]], [m[3]], marker = :circle, color = :red, label="", show = true)
            end
        end
    else
        THD = true
    end

    if THD
        println("Too high dimension -> we cannot visualize the roots :/ ")
    end

    return rts

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
    #The choice of r = 0.4 as a default value culd be argued against, it probably makes sense to perform tests to see if 0.4 is relevant. Shrinking r would accelerate the program, but may lead to the loss of some solutions
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

function main()

    t2 = 4.8537128789937585 
    t3 = 2.315641991877846
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    n = length(T)

    println("B susceptance matrix")

    #B = random_B(10, 3, true)
    B = [-12.4998    8.05323    4.44654;
    8.05323  -8.77734    0.724114;
    4.44654   0.724114  -5.17065]

    println("G0 loss matrix")
    #G0 = random_B(6, 3, true)
    G0 = [-2.56908    2.43434    0.134738;
    2.43434   -3.17263    0.738285;
    0.134738   0.738285  -0.873022]

    D1 = (-pi..pi)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) 

    println("k = 0")
    rts = perform(V, T, B, G0, D, 0)
    already_here = false
    i = 1
    L = collect(LinRange(0.025, 0.3, 12))
    for k in (L)
        println("k = $(round(k, digits = 3))")
        new_rts = []

        for root in rts
            d = next_interval(root)
            new_root = perform(V, T, B, G0, d, k, i)
            if length(new_root) > 1
                println("$root made babies")
            elseif length(new_root) == 0
                println("$root died")
            end
            if length(new_root) > 0
                for h in new_root
                    for g in new_rts
                        if !isempty(g.interval ∩ h.interval)
                            already_here = true
                        end
                    end
                    if !already_here
                        push!(new_rts, h)
                    end
                    already_here = false
                end
            end
        end
        rts = new_rts
        i += 1 
    end
    println("im done")
    readline()
end

function main_4D()

    t2 = rand()*2*pi
    t3 = rand()*2*pi
    t4 = rand()*2*pi
    println("Initial angles = $t2, $t3, $t4")
    T = [0., t2, t3, t4]

    V = [1., 1., 1., 1.]
    n = length(T)

    println("B susceptance matrix")
    B = random_B(12, 4, true)


    println("G0 loss matrix")
    G0 = random_B(8, 4, true)

    D1 = (-pi..pi)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) 

    println("k = 0")
    rts = perform(V, T, B, G0, D, 0)
    already_here = false
    i = 1
    L = collect(LinRange(0.025, 0.3, 12))
    for k in (L)
        println("k = $(round(k, digits = 3))")
        new_rts = []

        for root in rts
            d = next_interval(root)
            new_root = perform(V, T, B, G0, d, k, i)
            if length(new_root) > 1
                println("$root made babies")
            elseif length(new_root) == 0
                println("$root died")
            end
            if length(new_root) > 0
                for h in new_root
                    for g in new_rts
                        if !isempty(g.interval ∩ h.interval)
                            already_here = true
                        end
                    end
                    if !already_here
                        push!(new_rts, h)
                    end
                    already_here = false
                end
            end
        end
        rts = new_rts
        i += 1 
    end
    println("im done")
    readline()

end

main()