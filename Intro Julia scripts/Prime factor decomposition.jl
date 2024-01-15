function Decompo(n)
    # my first comment
    m=n
    j=2
    s=floor(n/2)
    X=[]

    if m%1 != 0
        println("You just tried to fool me didn't you ?")
    else

        while j<=s  
            if m%j==0
                push!(X,j)
                m=m/j
                j-=1
            end
            j+=1
        end
        
    end
    return X

end

    n=4.5
    D=Decompo(n)
    println(n)

    if isempty(D)
        println(" is a prime number .")
    else
    println(D)
    end
    println("")
