using IntervalRootFinding, LinearAlgebra, Plots, Colors

include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Solutions and losses-the functions.jl")
include("C:/Users/jim.delitroz/Documents/Julia scripts/src/Intro to PowerModels/Data modifier-the functions.jl")

function ShowTheCycle(L, angles)

    # Define polygon vertices and angles
    n = L # Number of vertices
    vertices = [(cos(2π*i/n)-1, sin(2π*i/n)+1) for i in 0:n-1]
    
    # Calculate endpoints of vectors
    vectors = [(0.4*cos(angle), 0.4*sin(angle)) for angle in angles]
    
    # Plot polygon and vectors
    plot(vertices, seriestype = :shape, fillalpha = 0, label="", aspect_ratio = :equal)
    for i in 1:n
        quiver!([vertices[i]], quiver=([vectors[i]]), color=:red, label="", show = true)
    end
end

function main()
    M = 17
    T0 = vector_angles(M, 0, 0.2)
    T1 = vector_angles(M, 1, 0.2)
    T2 = vector_angles(M, 2, 0.4)

    plot1 = ShowTheCycle(M, T2)

    colors = [:blue, :green, :red]
    B = my_cyclic_susceptance_matrix(M) # Add a true in the call of this fct to generate a new random B 
    display(B)
    V = []
    for i = 1:M
        push!(V, 0.95+rand()/10)
    end
    number_of_steps = 2
    mycolors = cgrad([:tomato2, :green], number_of_steps, categorical = true)
    i = 1
    for k in collect(LinRange(0., 0.9, number_of_steps))
         
        G = -B .* k
        for j in [2] #Any buses to compute their V-Q sensitivities
    
            x_axis_is_voltage = []
            y_axis_is_Q = []
            y_axis_is_Qwn = []
            y_axis_is_Qwn2 = []

            for v in collect(LinRange(0., 1.3, 27))
                V[j] = v
                Q = Reactive_Powers(V, T0, B, G)
          
                Q_wn = Reactive_Powers(V, T1, B, G)
                Q_wn2 = Reactive_Powers(V, T2, B, G)
                push!(x_axis_is_voltage, v)
                push!(y_axis_is_Q, Q[j])
                push!(y_axis_is_Qwn, Q_wn[j])
                push!(y_axis_is_Qwn2, Q_wn2[j])

            end
            plot!(x_axis_is_voltage, y_axis_is_Q, color = mycolors[i], label = "Q$j", legend = :topleft, show = true)
            plot!(x_axis_is_voltage, y_axis_is_Qwn, color = mycolors[i], linestyle = :dash, label = "Q$j wn", legend = :topleft)
            plot!(x_axis_is_voltage, y_axis_is_Qwn2, color = colors[i], linestyle = :dashdot, label = "Q$j double wn", legend = :topleft)
            i += 1
        end
    end 
    readline()
    #savefig("VQ curves.png")
end
main()