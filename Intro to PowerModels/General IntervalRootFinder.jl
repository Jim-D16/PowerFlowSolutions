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

function return_f(V, B, p1, p2, p3, q2, q3)
    sol = []
    eq1 = (x -> P1(V, [0, x[1], x[2]], B) - p1)
    eq2 = (x -> P2(V, [0, x[1], x[2]], B) - p2)
    eq3 = (x -> P3(V, [0, x[1], x[2]], B) - p3)
    eq4 = (x -> Q2(V, [0, x[1], x[2]], B) - q2)
    eq5 = (x -> Q3(V, [0, x[1], x[2]], B) - q3)
    push!(sol, eq1)
    push!(sol, eq2)
    push!(sol, eq3)
    push!(sol, eq4)
    push!(sol, eq5)
    return sol
end

function truncated_f(F, I) # For a fct F that returns a vector of size n, builds a form that will only return a vector of size |I|, by choosing the coords of F that corresp to elements of I
    t = length(I)
    G = (x -> SVector{t}([F(x)[i] for i in I]))
    return G
end


function perform(V, T, B, G0, k = 0)
    G = -G0 .* k
    
    p1 = P1(V, T, B, G)
    p2 = P2(V, T, B, G)
    p3 = P3(V, T, B, G)
    q2 = Q2(V, T, B, G)
    q3 = Q3(V, T, B, G)

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
    
    #println("G / B ratio : $k")
    
    F = (x -> SVector{5}([f(x) for f in return_f(V, B, p1, p2, p3, q2, q3)]))
    G23 = truncated_f(F, [2,3])

    D1 = (-pi..pi)
    D = IntervalBox(D1, D1)

    rts = roots(G23, D, Krawczyk)
   

    #k = round(k, digits = 3)

    if length(rts) == 0
        println("Ratio $k lead to 0 sol")
    else
        for r in rts
            m = mid(r.interval)
            if length(rts) > 2
                plot!([m[1]], [m[2]], seriestype = :scatter, color = :green, label="", show = true)
            else
                plot!([m[1]], [m[2]], seriestype = :scatter, color = :red, label="", show = true)
            end
        end
        if length(rts) == 1
            println("Critical value !? $k") # Closest candidat for the "jim standard 3-bus test case" k = 0.34860. Note that k = 0.34861 yields 0 solution
            println("Solution : $m")
            println("")
        end
    end
    println("k = $(round(k, digits = 3))")
    display(rts)
    readline()

    return length(rts)

end

function random_B(scale = 12)
    b12 = scale*rand()
    b13 = scale*rand()
    b23 = scale*rand()

    return [-(b12+b13) b12 b13; b12 -(b12+b23) b23; b13 b23 -(b13+b23)]
end
#=
function give_me_k(b23)

    V = [1., 1., 1.]
    #t2 = 2*rand()-2
    #t3 = 2*rand()-2
    
    t2 = pi/6
    t3 = -pi/10
    T = [0., t2, t3]



    #B = random_B()
    B = [-7.52 4.33 3.2; 4.33 -4.48 b23; 3.2 b23 -3.36]
    #B = [-27 12 15; 12 -32 20; 15 20 -35]

    B = round.(B, digits = 2)
    #display(B)
    #println("t2 = $t2, t3 = $t3")

    finalk = 0
    k = 0.3
    upper = 1
    lower = 0

    for i = 1:7
        Z = perform(V, T, B, k)

        if Z > 2
            lower = k
            k = (k + upper)/2 
        else
            upper = k
            k = (k+lower)/2
        end
        finalk = k
    end
    #println("im done")
    return finalk

#annotate!((0.95, 0.95), text("ratio = $k", :black, :right, 8))
end =#

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
    B = random_B(10)
    println("B susceptance matrix")
    display(B)

    #B = [-7.53 4.33 3.2; 4.33 -5.49 1.16; 3.2 1.16 -4.36]

    G0 = random_B(6)
    println("G0 loss matrix")
    display(G0)

    L = collect(LinRange(0, 1.5, 25))
    for k in (L)
        perform(V, T, B, G0, k)
    end
    println("im done")
    readline()
end

main()