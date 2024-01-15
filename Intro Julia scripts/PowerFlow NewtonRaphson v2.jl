
using LinearAlgebra
using ForwardDiff

# Data input
global V_1=1
global delta_1=0 #Phase and voltage of the slack bus


global Pd=[170 200 -238] #Power of every bus (apart from the slack bus). Positive sign indicates the load, so that the generator actually is in the last bus here
global Qd=[105.25 123.94] #Reactive power for every bus that is not voltage controlled, ie every PQ bus
global Vd=[1.02] #Voltage of every voltage controlled bus, ie PV bus

global N=1+length(Pd) # Number of buses
global M=length(Qd)+1 # Number of unknown voltages, +1 (so that the unknown voltages are V_2,V_3,...,V_m)

global Y=[8.99-im*44.84 -3.82+im*19.08 -5.17+im*25.85 0 
-3.82+im*19.08 8.99-im*44.84 0 -5.17+im*25.85
-5.17+im*25.85 0 8.19-im*40.86 -3.02+im*15.12
0 -5.17+im*25.85 -3.02+im*15.12 8.19-im*40.86] #Admittance matrix

V=ones(length(Qd))
pushfirst!(V,V_1)

V=vcat(V,Vd) # V now contains : the voltage of the slack bus, M-1 voltages, temporarily set to 1, that we want to determine, and finally the controlled voltages

delta=zeros(length(Pd))
pushfirst!(delta,delta_1) # delta now contains N terms, currently all equal to 0, but that we want to determine (except for the first term that will remain 0)


# We want to find the voltage corresponding to the PQ buses, as well as the phase angle for every bus (except the slack bus)


# Change of plans : have ForwardDiff compute the Jacobian for us


# Définition de la fonction g
function g(x)
    delta_here=[delta_1]
    for i=1:N-1
        push!(delta_here,x[i])
    end

    V_here=[V_1]
    for i=N:N+M-2
        push!(V_here,x[i])
    end
    V_here=vcat(V_here,Vd)

    Zp=PowerFlowP(delta_here,V_here,Y)
    Zq=PowerFlowQ(delta_here,V_here,Y)
    Z=vcat(Zp,Zq)
    return Z
end


function Jakobian(g,delta,V)
    u=unkn(delta,V)
    J=ForwardDiff.jacobian(g,u)
    return J
end

function unkn(delta,V)
    deltaX=copy(delta)
    VX=copy(V)
    splice!(deltaX,1)
    splice!(VX,1)
    if length(Vd)!=0
        for i=1:length(Vd)
            splice!(VX,length(VX))
        end
    end
    T=vcat(deltaX,VX)
    return T
end

#=
function Jacobian(delta,V,Y) # The output will be the Jacobian evaluated at the delta_i, |V_i| values that we input, in order to create the system that will give the wDelta_i and wV_i/V_i
    # j indexes the row, k the column 
    
        G=real(Y)
        B=imag(Y)
        
        J=zeros(N+M-2,N+M-2)
        for j=2:N 
            for k=2:N # Top-left quadrant of the Jacobian
                if k!=j
                    J[j-1,k-1]=V[j]*V[k]*(G[j,k]*sin(delta[k]-delta[j])-B[j,k]*cos(delta[k]-delta[j]))
                
                else
                    sumj=-V[j]^2*B[j,j]
                    for n=1:N
                        sumj+=V[j]*V[n]*(-G[j,n]*sin(delta[n]-delta[j])+B[j,n]*cos(delta[n]-delta[j]))
                    end
                    J[j-1,j-1]=sumj
                end
            end
            for k=2:M # Top-right quadrant
                k2=k+N-2
                if k!=j
                    J[j-1,k2]=V[j]*(G[j,k]*cos(delta[k]-delta[j])+B[j,k]*sin(delta[k]-delta[j]))
                else
                    sumj=V[j]*G[j,j]
                    for n=1:N
                        sumj+=V[n]*(G[j,n]*cos(delta[n]-delta[j])+B[j,n]*sin(delta[n]-delta[j]))
                    end
                    J[j-1,k-1]=sumj
                end 
            end
        end
    
        for j=2:M
            
            for k=2:N # Bottom-left quadrant
                if k!=j
                    J[j+N-2,k-1]=V[j]*V[k]*(-G[j,k]*cos(delta[k]-delta[j])-B[j,k]*sin(delta[k]-delta[j]))
                else
                    sumj=-V[j]^2*G[j,j]
                    for n=1:N
                        sumj+=V[j]*V[n]*(G[j,n]*cos(delta[n]-delta[j])+B[j,n]*sin(delta[n]-delta[j]))
                    end
                    J[j+N-2,k-1]=sumj
                end
            end
            for k=2:M # Bottom-right quadrant
                if k!=j
                    J[j+N-2,k+N-2]=V[j]*(G[j,k]*sin(delta[k]-delta[j])-B[j,k]*cos(delta[k]-delta[j]))
                else
                    sumj=-V[j]*B[j,j]
                    for n=1:N
                        sumj+=V[n]*(G[j,n]*sin(delta[n]-delta[j])-B[j,n]*cos(delta[n]-delta[j]))
                    end
                    J[j+N-2,k+N-2]=sumj
                end 
    
            end
        end
        return J
end

=#
    
function PowerFlowP(delta, V, Y)

    G=real(Y)
    B=imag(Y)
  
    P=[]
    for i=2:N
        psum_i=0
        for n=1:N
                psum_i+=V[i]*V[n]*(G[i,n]*cos(delta[n]-delta[i])+B[i,n]*sin(delta[n]-delta[i]))
                #println("psum_$i = $psum_i")
        end
        push!(P,psum_i)
        #println("P_$i = $psum_i)")
    end
    return P

end

function PowerFlowQ(delta, V, Y)

    G=real(Y)
    B=imag(Y)

    Q=[]

    for i=2:M
        qsum_i=0
        for n=1:N
                qsum_i+=V[i]*V[n]*(G[i,n]*sin(delta[n]-delta[i])-B[i,n]*cos(delta[n]-delta[i]))
        end
        push!(Q,qsum_i)
        #println("Q_$i = $qsum_i)")
    end
    return Q

end





function Mismatches(delta,V,Y,P,Q) # The output vector contains wP_1,...,wP_N,wQ_1,...,wQ_n where wP_i is DELTA P_i, the diff between our last estimate and the truth (Mismatches)
    Piter=PowerFlowP(delta,V,Y)
    Qiter=PowerFlowQ(delta,V,Y)
    #println("En direct de la fct Mismatches : P= $P")
    #println("En direct de la fct Mismatches : Piter= $Piter")

    wP=zeros(N-1)
    wQ=zeros(M-1)

    for i=1:N-1
       wP[i]=complex(Piter[i]-P[i])
    end

    for i=1:M-1
        wQ[i]=complex(Qiter[i]-Q[i])
    end

    T=vcat(wP,wQ)
    return T
end

function NewtonRaphson(delta,V,Y,P,Q)
    t=0.5
    t0=10
    delta0=copy(delta)
    V0=copy(V)
    Iter=0

    while(t0>t && Iter<5)
        Iter+=1
        #P0=PowerFlowP(delta0,V0,Y)
        #Q0=PowerFlowQ(delta0,V0,Y)
        T0=Mismatches(delta0,V0,Y,P,Q)
        J=Jakobian(g,delta0,V0)
        X=J\T0
        #println(X)
        #print("The mismatches vector looks like $T0 ")
        print("The norm of the mismatches is $t0 ")
        println("after $Iter iterations")

        for i=2:N
            delta0[i]+=X[i-1] 
            #println("Iteration $i of finding the solutions is ok")
        end
        q=2
        for i=N:length(X)
            V0[q]+=X[i]
            q+=1
        end

        t0=norm(T0)
        
    end

    if Iter==5
        println("$Iter iterations have passed, which is a lot. Here are our current solutions :")
        println("V = $V0")
        println("delta = $delta0")
    else
        println("The method converged in $Iter iterations to the following :")
        println("V = $V0")
        println("delta = $delta0")
    end

end

NewtonRaphson(delta,V,Y,Pd,Qd) 
