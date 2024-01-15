#We want to solve a system of linear eqns, that is find the vector X st for A, Y given we have AX=y


import Pkg 
import RowEchelon

function Mprod(A,B)
    n=size(B,1)
    p=size(B,2)
    m=size(A,1)
    Y=zeros(m,p)
    if size(A,2)!=n
        println("Dimension issue, sorry")
    else
        for a=1:m*p
            t=0
            q=floor(Int, a/m)
            r=a%m
            println("Je suis en train de calculer le $a i eme terme")
            if r==0
                q-=1
                r+=m
            end
            for k=1:n
                t+=A[r+(k-1)*m]*B[(q)*n+k]
            end
            Y[a]=t
        end
    end
    return Y
end



A=[0 2 1
-1 0.5im 2im]

X=[2+im 4-3im
-1 5
6 3]

AA=imag(A)

println(AA)
                                                                                                                                                                                          

        

