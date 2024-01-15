# Let's work on N buses indexed by i 
# What do we need to input ? An NxN matrix Y, and the actual numerical values of P_i and Q_i for every i


# The way to recollect P and Q once the |V_i| and delta_i are known (as well as the Y_ij of course)

function PowerFlowP(delta, V, Y)
    N=length(delta)
    if(N==length(V) && size(Y)==(N,N))
        #println("Dimensions : check")
    else exit()
    end

    G=real(Y)
    Theta=zeros(N,N)
    for l=1:N^2
        Theta[l]=angle(Y[l])
    end
    println(Theta)
    Ym=zeros(N,N)
    for l=1:N^2
        Ym[l]=abs(Y[l])
    end
    P=[]
    for i=1:N
        psum_i=0
        for n=1:N
            if n!=i
                psum_i+=V[i]*V[n]*Ym[i,n]*cos(Theta[i,n]+delta[n]-delta[i])
                #println("psum_$i = $psum_i")
            end
        end
        push!(P,(V[i]^2*G[i,i]+psum_i))
        #println("P_$i = $(V[i]^2*G[i,i]+psum_i)")
    end
    return P

end

function PowerFlowQ(delta, V, Y)
    N=length(delta)
    if(N==length(V) && size(Y)==(N,N))
        #println("Dimensions : check")
    else exit()
    end
    B=imag(Y)
    Theta=zeros(N,N)
    for l=1:N^2
        Theta[l]=angle(Y[l])
    end
    Ym=zeros(N,N)
    for l=1:N^2
        Ym[l]=abs(Y[l])
    end
    Q=[]

    for i=1:N
        qsum_i=0
        for n=1:N
            if n!=i
                qsum_i+=V[i]*V[n]*Ym[i,n]*sin(Theta[i,n]+delta[n]-delta[i])
            end
        end
        push!(Q,-V[i]^2*B[i,i]-qsum_i)
        #println("Q_$i = $(-V[i]^2*B[i,i]-qsum_i)")
    end
    return Q

end

function Mismatches(delta,V,Y,P,Q) # The output vector contains wP_1,...,wP_N,wQ_1,...,wQ_n where wP_i is DELTA P_i, the diff between our last estimate and the truth (Mismatches)
    N=length(V)
    Pinit=PowerFlowP(delta,V,Y)
    Qinit=PowerFlowQ(delta,V,Y)
    #println("En direct de la fct Mismatches : P= $P")
    #println("En direct de la fct Mismatches : Pinit= $Pinit")

    wP=zeros(N)
    wQ=zeros(N)
    wP=complex(wP)
    wQ=complex(wQ)
    #println(typeof(wP))

    for i=1:N
       wP[i]=complex(P[i]-Pinit[i])
        wQ[i]=complex(Q[i]-Qinit[i])
    end
    T=vcat(wP,wQ)
    #println(T)
    return T
end

function Jacobian(delta,V,Y) # The output will be the Jacobian evaluated at the delta_i, |V_i| values that we input, in order to create the system that will give the wDelta_i and wV_i/V_i
# j indicates the row, k the column so that (2,4) will be on the 2nd row, 4th column

    N=length(V)
    G=real(Y)
    B=imag(Y)
    Theta=zeros(N,N)
    for l=1:N^2
        Theta[l]=angle(Y[l])
    end
    Ym=zeros(N,N)
    for l=1:N^2
        Ym[l]=abs(Y[l])
    end



    J=zeros(2N,2N)
    for j=1:N 
        for k=1:N # Top-left NxN quadrant of the Jacobian
            if k!=j
                J[j,k]=-V[j]*V[k]*Ym[j,k]*sin(Theta[k,j]+delta[k]-delta[j])
            
            else
                sumj=0
                for n=1:N
                    sumj+=V[j]*V[n]*Ym[j,n]*sin(Theta[k,j]+delta[n]-delta[j])
                end
                J[j,j]=sumj
            end
        end
        for k=(N+1):2N # Top-right quadrant
            k2=k-N
            if k2!=j
                J[j,k]=V[k2]*V[j]*Ym[j,k2]*cos(Theta[k2,j]+delta[k2]-delta[j])
            else
                sumj=2*V[j]*G[j,j]
                for n=1:N
                    sumj+=V[n]*Ym[j,n]*cos(Theta[n,j]+delta[n]-delta[j])
                end
                sumj-=V[j]*Ym[j,j]
                J[j,k]=V[j]*sumj
            end 
        end
    end

    for j=(N+1):2N
        j2=j-N
        for k=1:N # Bottom-left quadrant
            if k!=j2
                J[j,k]=complex(-V[j2]*V[k]*Ym[j2,k]*cos(Theta[k,j2]+delta[k]-delta[j2]))
            else
                sumj=0
                for n=1:N
                    sumj+=V[j2]*V[n]*Ym[j2,n]*cos(Theta[n,j2]+delta[n]-delta[j2])
                end
                sumj-=V[j2]^2*Ym[j2,j2]
                J[j,k]=complex(sumj)
            end
        end
        for k=(N+1):2N # Bottom-right quadrant
            j2=j-N
            k2=k-N
            if k2!=j2
                J[j,k]=complex(-V[k2]*V[j2]*Ym[j2,k2]*sin(Theta[k2,j2]+delta[k2]-delta[j2]))
            else
                sumj=2*V[j2]*B[j2,j2]
                for n=1:N
                    sumj+=V[n]*Ym[j2,n]*sin(Theta[n,j2]+delta[n]-delta[j2])
                end
                J[j,k]=complex(-V[j2]*sumj)
            end 

        end
    end
    return J
end


function NewtonRaphson(delta,V,Y,P,Q)
    t=0.01
    N=length(delta)
    t0=10
    delta0=complex(delta)
    V0=complex(V)
    Iter=0

    while(t0>t && Iter<12)
        Iter+=1
        #P0=PowerFlowP(delta0,V0,Y)
        #Q0=PowerFlowQ(delta0,V0,Y)
        T0=Mismatches(delta0,V0,Y,P,Q)
        J=Jacobian(delta0,V0,Y)

        Z0=J\T0
        #println(Z0)
        print("The mismatches vector looks like $T0 ")
        println("after $Iter iterations")

        for i=1:N
            delta0[i]+=Z0[i] # C'est complètement faux. Z0 contient les corrections, et je le considère comme les vraies valeurs
            V0[i]+=V0[i]*Z0[i+N]
            #println("Iteration $i of finding the solutions is ok")
        end
        
    end

    if Iter==12
        println("12 iterations have passed, which is a lot. Here are our current solutions :")
        println("V = $V0")
        println("delta = $delta0")
    else
        println("The method converged in $Iter iterations to the following :")
        println("V = $V0")
        println("delta = $delta0")
    end

end





#Data from Grainger94, p338, example 9.2
Y=[8.99-im*44.84 -3.82+im*19.08 -5.17+im*25.85 0 
-3.82+im*19.08 8.99-im*44.84 0 -5.17+im*25.85
-5.17+im*25.85 0 8.19-im*40.86 -3.02+im*15.12
0 -5.17+im*25.85 -3.02+im*15.12 8.19-im*40.86]

P=[50,170,200,-238]
Q=[30.99,-05.35,123.94,49.58]


deltaInit=[1 1.5 0.8 0.3]
VInit=[1.5, 1.7, 1.2, 1.4]

NewtonRaphson(deltaInit,VInit,Y,P,Q)