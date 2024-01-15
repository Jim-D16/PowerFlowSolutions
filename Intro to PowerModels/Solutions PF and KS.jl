using LinearAlgebra, StaticArrays, IntervalArithmetic, IntervalRootFinding, Dates, ForwardDiff

global k = 0.2


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

function PF_model()
    B = [-27 12 15; 12 -32 20; 15 20 -35] 
    #B = [-20 12 8; 12 -35 23; 8 23 -31] 
    
    G = -B .* k

    V = [1., 1., 1.]
    T = [0., -0.728, 0.923]

    p1 = P1(V, T, B, G)
    p2 = P2(V, T, B, G)
    p3 = P3(V, T, B, G)
    q2 = Q2(V, T, B, G)
    q3 = Q3(V, T, B, G)
    
    println("P1 = $(p1)")
    println("P2 = $(p2)")
    println("P3 = $(p3)")
    println("Q2 = $(q2)")
    println("Q3 = $(q3)")
    println("Losses : $(p1 + p2 + p3)")
    
    
end

function return_f(V, B, p1, p2, p3, q2, q3)
    #k = 0.2
    G = -B .* k
    sol = []
    eq1 = (x -> P1(V, [0, x[1], x[2]], B, G) - p1)
    eq2 = (x -> P2(V, [0, x[1], x[2]], B, G) - p2)
    eq3 = (x -> P3(V, [0, x[1], x[2]], B, G) - p3)
    eq4 = (x -> Q2(V, [0, x[1], x[2]], B, G) - q2)
    eq5 = (x -> Q3(V, [0, x[1], x[2]], B, G) - q3)
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


function KS_model(v2,v3,t2,t3)
    #k = 0.2


    B = [-27 12 15; 12 -32 20; 15 20 -35]
    G = -B .* k
    V = [1., v2, v3]
    T = [0., t2, t3]

    p1 = P1(V, T, B, G)
    p2 = P2(V, T, B, G)
    p3 = P3(V, T, B, G)
    q2 = Q2(V, T, B, G)
    q3 = Q3(V, T, B, G)

   # V = [1., 1., 1.] # This line constrains every voltage at 1, ie reduces the PF equations to some KS equations

    F = (x -> SVector{5}([f(x) for f in return_f(V, B, p1, p2, p3, q2, q3)]))
    G23 = truncated_f(F, [2,3])

    D1 = (-pi..pi)
    D = IntervalBox(D1, D1)
   
    rts23 = roots(G23, D, Krawczyk)
    #rts13 = roots(G13, D, Krawczyk)
    display(rts23)
    println("")

    return rts23

    #=
    for i in rts12
        for j in rts13
            k = intersect(i.interval,j.interval)
            push!(A,k)
        end
    end
    println("Any solution to the system will lie in $A")
    =#

    

end
PF_model()
KS_model(1., 1., -0.728, 0.923)

# Plusieurs tests sur le KS model pour espérer en trouver un qui n'a pas de solution
#=
for i = 1:5

    v2 = rand()/100 + 0.95
    v3 = rand()/100 + 0.95
    t2 = rand()*2*pi - pi
    t3 = rand()*2*pi - pi

    if length(KS_model(v2, v3, t2, t3)) == 0
        println("The setup yielded 0 solution : ")
        println("V2 = $v2, V3 = $v3, T2 = $t2, T3 = $t3")
    end
    println("Iter $i over")
end

for i = 6:10

    v2 = rand()/100 + 1.04
    v3 = rand()/100 + 0.95
    t2 = rand()*2*pi - pi
    t3 = rand()*2*pi - pi

    if length(KS_model(v2, v3, t2, t3)) == 0
        println("The setup yielded 0 solution : ")
        println("V2 = $v2, V3 = $v3, T2 = $t2, T3 = $t3")
    end
    println("Iter $i over")
end

for i = 11:15

    v2 = rand()/100 + 0.95
    v3 = rand()/100 + 1.04
    t2 = rand()*2*pi - pi
    t3 = rand()*2*pi - pi

    if length(KS_model(v2, v3, t2, t3)) == 0
        println("The setup yielded 0 solution : ")
        println("V2 = $v2, V3 = $v3, T2 = $t2, T3 = $t3")
    end
    println("Iter $i over")
end

for i = 16:20

    v2 = rand()/100 + 1.04
    v3 = rand()/100 + 1.04
    t2 = rand()*2*pi - pi
    t3 = rand()*2*pi - pi

    if length(KS_model(v2, v3, t2, t3)) == 0
        println("The setup yielded 0 solution : ")
        println("V2 = $v2, V3 = $v3, T2 = $t2, T3 = $t3")
    end
    println("Iter $i over")
end
=#
