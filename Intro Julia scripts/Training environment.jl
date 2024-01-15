using Printf

function Max(a,b)
    if abs(a)>abs(b)
        return a
    elseif abs(b)>abs(a)
        return b
    elseif angle(b)>angle(a)
        return b
    else 
        return a
    end
end


function Max(a)
    m=0
    n=length(a)
    for i=1:n
        m=Max(m,a[i])
    end
    return m
end

function iMax(a)
    m=0
    j=0
    n=length(a)
    for i=1:n
        if a[i]==Max(a)
            j=i
        end
    end
    return j
end

function randPM(x) #returns a random real number between -x and x
    return 2x*(-0.5+rand())
end

function jround(A, n) #rounds every entry of an array A to n digits
    B=[]
    k=length(A)
    for i=1:k
        t=round(A[i], digits=n)
        push!(B, t)
    end
    return B
end

 function Sort(A)
    B=[]
    C=A
    n=length(A)
    for i=1:n
        pushfirst!(B, Max(C))
        splice!(C,iMax(C))
    end
    return B
end 



V=Array{Complex}(undef,5)
for j in 1:5
    V[j]=randPM(10)+ randPM(10)*im
end

 W=jround(V, 3)

p=iMax(V)

println(W)
println("And now ! $(Sort(W))")



