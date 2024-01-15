using Printf

function f(x)
    y=x^4+2x^3+3x^2-x-2.5
    return y
end

function f2(x)
    y=4x^3+6x^2+6x-1
    return y
end


#=
function AffineRoot(g,y) #Donnez une fonction g et un nombre y tq g(x) est proche de 0 for some x<y. La fonction return a tq g(a)<0.2

p=0.1 #step
t=0.2 #tolerance
#y=1 #Plus grand entier après le zéro de la fonction que l'on recherche
z=1
global a=y

while z>t

    global z=abs(g(a))
    if z>t
    global a-=p
    else
        break
    end
end

println("Final result : g($a)=$(g(a))")
return a
end

function AffineRoot(g)
    return AffineRoot(g,2)
end 
=#

function NewtonRaphson(g,h,x) #Bien que la fct rende un résultat pour n'importe quelle paire de fcts g,g' , ce résultat ne fait sens que si g' est effectivement la dérivée de g
    x=convert(AbstractFloat, x)
    X=[x]
    k=1
    z=3.14
    while abs(g(z))>0.01
        
         z=X[k]-g(X[k])/h(X[k])
        push!(X,z)
        
        println("Iteration $(k-1)") 
        k+=1

        if k==30
            @printf("Il semble impossible de converger : g(%.3f) = %.3f \n", z, g(z))
            break
        end 

    end
    k-=1
    z=X[k]-g(X[k])/h(X[k])
@printf("J'ai quitté la boucle car g(%.3f) = %.5f \n", z, g(z))    
println("Nombre d'itérations : $(k)")
@printf("x_ %d = %.3f \n", k, z)

return z
end


NewtonRaphson(f, f2, 13)
