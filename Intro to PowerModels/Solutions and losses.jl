using PowerModels, LinearAlgebra, StaticArrays, IntervalArithmetic, IntervalRootFinding, JuMP, Ipopt, ProgressBars
using Colors, ForwardDiff, Dates, Plots

include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Solutions and losses-the functions.jl")

function perform(V, T, B, G0, D, k = 0, i = 1) #Every computation made from this script had the "correct?" sign of G, that is, G was the opposite of B which is GREAT.
    n = length(V)
    G = -G0 .* k
    
    P = Powers(V, T, B)

    
    F = (x -> SVector{n-1}([f(x) for f in return_f_KS(V, P, B, G)]))


    rts = roots(F, D, Krawczyk, 1e-4)
    
    numberofcolors = 28
    THD = false
    if i > numberofcolors
        i = numberofcolors
    end
    mycolors = cgrad([:tomato2, :green], numberofcolors, categorical = true)
    if n == 3
        for r in rts
            
            m = mid(r.interval)
            T2 = [[0]; m]
            Pnow = Powers(V, T2, B, G)
            total_loss = 0
            for p in Pnow
                total_loss += p
            end
            println("Total loss : $(round(total_loss, digits = 3))")
            if total_loss < -0.01
                println("$r Not physically meaningful")
            else
                display(r)
                if k > 0
                    plot!([m[1]], [m[2]], marker = :circle, color = mycolors[i], label="", show = true)
                else
                    plot!([m[1]], [m[2]], marker = :circle, markersize = 5, color = :red, label="", show = true)
                end
            end

        end
    elseif n == 4
        for r in rts
            m = mid(r.interval)
            T2 = [[0]; m]
            Pnow = Powers(V, T2, B, G)
            total_loss = 0
            for p in Pnow
                total_loss += p
            end
            println("Total loss : $total_loss")
            if total_loss < -0.01
                println("$r Not physically meaningful")
            else
                display(r)
                if k > 0
                    plot3d!([m[1]], [m[2]], [m[3]], marker = :circle, color = mycolors[i], label="", show = true)
                else
                    plot3d!([m[1]], [m[2]], [m[3]], marker = :circle, markersize = 5, color = :red, label="", show = true)
                end
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


function main()
    initialize_background()

    t2 = rand()*2*pi
    t3 = rand()*2*pi
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    n = length(T)

    println("B susceptance matrix")

    B = random_B(10, 3, true)
    #=
    B = [-2.3     1    1.3;
    1  -2.8     1.8;
    1.3    1.8  -3.1] 
    =#


    println("G0 loss matrix")
    G0 = random_B(6, 3, true)
    #=
    G0 = [-12.5  6    6.5;
    6  -6.5  0.5;
    6.5   0.5  -7] 
    =#

    P = Powers(V, T, B)
    println("P2 = $(P[2])")
    println("P3 = $(P[3])") 
    
    D1 = (-pi/2..pi/2)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) 

    rts = perform(V, T, B, G0, D, 0)
    println("Every lossless solution has been found ;-)")
    valid = 0
    for root in rts
        m = mid(root.interval)
        if (cos(m[1]) >= 0 && cos(m[2]) >= 0 && cos(m[1]-m[2]) >= 0)
            valid += 1
        end
    end
    if valid > 1
        readline()
    end

#=


    already_here = false
    i = 1
    L = collect(LinRange(0, 1, 11))
    for k in (L)
        valid = 0
  
        rts = perform(V, T, B, G0, D, k, i)
        println("k = $(round(k, digits = 3))  Number of sol : $(length(rts))")
        for root in rts
            m = mid(root.interval)
            if (cos(m[1]) >= 0 && cos(m[2]) >= 0 && cos(m[1]-m[2]) >= 0)
                valid += 1
            end
        end
        #=new_rts = []

        for root in rts
            d = next_interval(root, 0.4)
            new_root = perform(V, T, B, G0, d, k, i) # 
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
        i += 1 =#

    end
    println("im done")
    readline() 
=#
end

function main_4D()

    t2 = 4.540823929597235
    t3 = 4.952715884631364
    t4 = 5.824477533260155
    println("Initial angles = $t2, $t3, $t4")
    T = [0., t2, t3, t4]

    V = [1., 1., 1., 1.]
    n = length(T)

    println("B susceptance matrix")
    #B = random_B(12, 4, true)
    B = [ -11.5312       10.4033     1.11996    0.00793174;
    10.4033      -21.758      5.50708    5.84768;
     1.11996       5.50708  -16.6064     9.97934;
     0.00793174    5.84768    9.97934  -15.8349]


    println("G0 loss matrix")
    #G0 = random_B(8, 4, true)
    G0 = [-18.018      3.01572    7.80488    7.19738;
    3.01572  -15.3365     5.13042    7.1904;
    7.80488    5.13042  -18.1905     5.25517;
    7.19738    7.1904     5.25517  -19.6429]

  
    D2 = IntervalBox([interval(-1.47327, -1.47325), interval(-1.07745, -1.07743), interval(-0.20256, -0.202551)])
    D1 = IntervalBox([interval(-1.74241, -1.74231), interval(-1.33051, -1.33043), interval(-0.458743, -0.458672)])

    #=mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) =#

    println("k = 0")
    rts = perform(V, T, B, G0, D1, 0)
    println(typeof(rts))
    rts2 = perform(V, T, B, G0, D2, 0)
    rts = [rts; rts2]

    already_here = false
    i = 1
    L = collect(LinRange(0.025, 0.2, 25))
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

function main_yielding_a_contradiction_bis() # 3 bus system that USED TO contradict the assumption (losses increase -> solutions diminish), but does not actually
    #The accelerated method (ie the use of next_interval) supposes that the assumption is true, in particular, using it here will not yield a contradiction (although there is one)
    initialize_background()
    t2 = 0
    t3 = 3/2*pi
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    B = [-7 4 3; 4 -5 1; 3 1 -4]
    println("B susceptance matrix")
    display(B)
    n = length(T)

    G0 = [-8   5   3; # Contradicts the monotony of number of solutions with k : k=0 admits 2 sol; 0.2<=k<= admits 4
    5  -7   2;
    3   2  -5]   
    
    println("G0 opposite of (unscaled) loss matrix")
    display(G0)


    P = Powers(V, T, B)
    println("P2 = $(P[2])")
    println("P3 = $(P[3])") 

    D1 = (-pi..pi)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) 

    i = 1
    L = collect(LinRange(0, 2.5, 26))
    for k in (L)
        println("k = $(round(k, digits = 2))")
        perform(V, T, B, G0, D, k, i)
        i += 1
    end
    println("im done")
    readline()
end

function main_searching_a_contradiction()
    #The accelerated method (ie the use of next_interval) supposes that the assumption is true, in particular, using it here will not show a contradiction (even if there is one)
    initialize_background()
    t2 = rand()*2*pi
    t3 = rand()*2*pi
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    println("B susceptance matrix")
    B = random_B(12, 3, true)

    n = length(T)

    spot = 100
    
    println("G0 opposite of (unscaled) loss matrix")
    G0 = random_B(10, 3, true)


    P = Powers(V, T, B)
    println("P2 = $(P[2])")
    println("P3 = $(P[3])") 

    D1 = (-pi..pi)
    mybox = []
    for i = 2:n
        push!(mybox, D1)
    end

    D = IntervalBox(mybox) 

    i = 1
    L = collect(LinRange(0, 2.5, 26))
    for k in (L)
        println("k = $(round(k, digits = 2))")
        if (l = length(perform(V, T, B, G0, D, k, i))) > spot
            println("k reached value $k. Consequently, our number of solutions jumped from $spot to $l")
            readline()
        elseif l >= 4
            break
        end
        spot = l
        i += 1
        
    end
    println("im done")
end
for i = 1:15
    main()
end
exit()

