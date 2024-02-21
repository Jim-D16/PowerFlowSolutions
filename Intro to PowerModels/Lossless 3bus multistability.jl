# The point of this program is to find parameters so that my 3bus system admits multiple lossless equilibrium points, with <pi/2 angular differences. 
# With the change of variables, my 2 power equations (P2 and P3) each define a plane, that intersect in a line.
# I'm looking for parameters so that this line intersects the arcsine Pringles in two points
# Translation : a = B_12, b = B_13, c = B_23, x = sin(theta_2), y = sin(theta_3), z = sin(theta_2 - theta_3)

# Answer !
# a = 0.48
# b = 0.49
# c = -0.31
# p1 = -0.08
# p2 = -0.04

# Better answer !
# There is no solution with c\geq 0 (see overleaf), therefore no meaningful solution


function main()

    a = round(rand(), digits = 2)
    b = round(rand(), digits = 2)
    c = round(rand(), digits = 2)
    p1 = round(rand() - 0.5, digits = 2)
    p2 = round(rand() - 0.5, digits = 2)

    #=
    a = 0.48
    b = 0.49
    c = -0.31
    p1 = -0.08
    p2 = -0.04
    =#
    
  
    tol = 1e-3
    temoin = -1
    Z = []

    for z in collect(LinRange(-1, 1, 201))
        if abs((p1 - c*z)/a) <= 1 && abs((p2 + c*z)/b) <= 1
            if abs(asin(z) - asin((p1 - c*z)/a) + asin((p2 + c*z)/b)) < tol
                #=println("x = $((p1 - c*z)/a)")
                println("y = $((p2 + c*z)/b)")
                println("z = $z")
                println("asin(x) = $(asin((p1 - c*z)/a))")
                println("asin(y) = $(asin((p2 + c*z)/b))")
                println("asin(z) = $(asin(z))") =#
                push!(Z, z)
            end
        end
    end
    if length(Z) > 1
        println("Values of z that are in the Pringles : $(round.(Z, digits = 2))")
        println(" a = $a")
        println(" b = $b")
        println(" c = $c")
        println(" p1 = $p1")
        println(" p2 = $p2")
        println("")
    end


end

for i = 1:15
    main()
end

