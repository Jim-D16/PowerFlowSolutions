using PowerModels, LinearAlgebra, StaticArrays, IntervalArithmetic, IntervalRootFinding, JuMP, Ipopt, ProgressBars
using Colors, ForwardDiff, Dates, Plots

include("C:/Users/jim.delitroz/Documents/Julia scripts/Intro to PowerModels/Solutions and losses-the functions.jl")

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
    mycolors = cgrad([:tomato2, :green], 52, categorical = true)
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


function main()

    t2 = 4.475179362142027
    t3 = 0.4610013805421309
    println("Initial angles = $t2, $t3")
    T = [0., t2, t3]
    V = [1., 1., 1.]
    n = length(T)

    println("B susceptance matrix")

    #B = random_B(10, 3, true)
    B = [-12.4522     5.70476    6.74743;
    5.70476  -12.4265     6.72173;
    6.74743    6.72173  -13.4692]


    println("G0 loss matrix")
    #G0 = random_B(6, 3, true)
    G0 = [-7.9739    4.99974    2.97416;
    4.99974  -5.92632    0.926577;
    2.97416   0.926577  -3.90074]

    
   

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
    L = collect(LinRange(0, 2, 50))
    for k in (L)
        println("k = $(round(k, digits = 3))")
        new_rts = []

        for root in rts
            d = next_interval(root, 0.4)
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

main_4D()