using PowerModels, LinearAlgebra, StaticArrays, IntervalArithmetic, IntervalRootFinding, JuMP, Ipopt, ProgressBars, ForwardDiff, Dates, Plots

function P1(V, T, B, G = B.*0)
    V[1]*(V[2]*(B[1,2]*sin(T[1]-T[2]) + G[1,2]*cos(T[1]-T[2])) + V[3]*(B[1,3]*sin(T[1]-T[3]) + G[1,3]*cos(T[1]-T[3])) + V[1]*G[1,1])
end

function P2(V, T, B, G = B.*0)
    V[2]*(V[1]*(B[2,1]*sin(T[2]-T[1]) + G[2,1]*cos(T[2]-T[1])) + V[3]*(B[2,3]*sin(T[2]-T[3]) + G[2,3]*cos(T[2]-T[3])) + V[2]*G[2,2])
end

function P3(V, T, B, G = B.*0)
    V[3]*(V[2]*(B[3,2]*sin(T[3]-T[2]) + G[3,2]*cos(T[3]-T[2])) + V[1]*(B[3,1]*sin(T[3]-T[1]) + G[3,1]*cos(T[3]-T[1])) + V[3]*G[3,3])
end

function Q2(V, T, B, G = B.*0)
    -V[2]*(V[1]*(B[2,1]*cos(T[2]-T[1]) - G[2,1]*sin(T[2]-T[1])) + V[3]*(B[2,3]*cos(T[2]-T[3]) - G[2,3]*sin(T[2]-T[3])) + V[2]*B[2,2])
end

function Q3(V, T, B, G = B.*0)
    -V[3]*(V[1]*(B[3,1]*cos(T[3]-T[1]) - G[3,1]*sin(T[3]-T[1])) + V[2]*(B[3,2]*cos(T[3]-T[2]) - G[3,2]*sin(T[3]-T[2])) + V[3]*B[3,3])
end

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

function truncated_f(F, I) # For a fct F that returns a vector of size n, builds a form that will only return a vector of size |I|, by choosing the coords of F that corresp to elements of I
    t = length(I)
    G = (x -> SVector{t}([F(x)[i] for i in I]))
    return G
end


function perform(V, T, B, G0, k = 0)
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


    D1 = (-pi..pi)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) #The slack bus is assumed to have index 1

    rts = roots(F, D, Krawczyk, 1e-6)

    
    println("k = $k")
    display(rts)

   

    #k = round(k, digits = 3)
    THD = false
    if length(rts) == 0
        println("Ratio $k lead to 0 sol")
    elseif n == 3
        for r in rts
            m = mid(r.interval)
            if length(rts) > 2
                plot!([m[1]], [m[2]], seriestype = :scatter, color = :green, label="", show = true)
            else
                plot!([m[1]], [m[2]], seriestype = :scatter, color = :red, label="", show = true)
            end
        end
        if length(rts) == 1
            println("Critical value !? $k") 
            println("Solution : $m")
            println("")
        end
    elseif n == 4
        for r in rts
            m = mid(r.interval)
            if length(rts) > 2
                plot3d!([m[1]], [m[2]], [m[3]], marker = :circle, color = :green, label="", show = true)
            else
                plot3d!([m[1]], [m[2]], [m[3]], marker = :circle, color = :red, label="", show = true)
            end
        end
        if length(rts) == 1
            println("Critical value !? $k") 
            println("Solution : $m")
            println("")
        end
        println("Solution for k = $k is 3D-plotted")
    else
        THD = true
    end
    if THD
        println("Too high dimension -> we cannot visualize the roots :/ ")
    end
    #println("k = $(round(k, digits = 3))")
    #display(rts)

    return length(rts)

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


function main_yielding_a_contradiction()
    t2 = 6.06785081987959
    t3 = 4.028180014959005
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    B = [-7.53 4.33 3.2; 4.33 -5.49 1.16; 3.2 1.16 -4.36]


    G0 = [-8.07466   4.95559   3.11907; # Contradicts the monotony of number of solutions with k : k=0 admits 2 sol; 0.2<=k<= admits 4
    4.95559  -6.77213   1.81654;
    3.11907   1.81654  -4.93561]
    
    #G0 = random_B(6)
    
    println("G0 loss matrix")
    display(G0)
    L = collect(LinRange(0, 1, 11))
    for k in (L)
        perform(V, T, B, G0, k)
    end
    println("im done")
    readline()
end

function main()

    t2 = rand()*2*pi
    t3 = rand()*2*pi
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]

    println("B susceptance matrix")

    B = random_B(10, 3, true)
    #B = [-7.53 4.33 3.2; 4.33 -5.49 1.16; 3.2 1.16 -4.36]

    println("G0 loss matrix")
    G0 = random_B(6, 3, true)

    L = collect(LinRange(0, 1.5, 25))
    for k in (L)
        perform(V, T, B, G0, k)
    end
    println("im done")
    readline()
end

function main_4D()

    t2 = rand()*2*pi
    t3 = rand()*2*pi
    t4 = rand()*2*pi
    #println("Initial angles = $t2, $t3, $t4")
    #T = [0., t2, t3, t4]
    T = [0., 4.745265081210639, 3.5352214556297867, 0.6953053086777519]
    V = [1., 1., 1., 1.]

    println("B susceptance matrix")
    #B = random_B(12, 4, true)

    B = [-11.8462     6.40909   4.40631     1.03082;
    6.40909  -18.7585    2.56043     9.78902;
    4.40631    2.56043  -7.48118     0.514437;
    1.03082    9.78902   0.514437  -11.3343]

    println("G0 loss matrix")
    #G0 = random_B(8, 4, true)

    G0 = [-9.2289    4.34688   1.95079   2.93123;
    4.34688  -8.93819   1.85527   2.73604;
    1.95079   1.85527  -7.33457   3.52851;
    2.93123   2.73604   3.52851  -9.19578]

    L = collect(LinRange(0, 0.3, 14))
    for k in (L)
        perform(V, T, B, G0, k)
    end
    println("im done")
    readline()

end

main_4D()