function Admittance(mydata) # La matrice est "ordonnée" selon les index des bus, ie 1e colonne bus 1, 2e colonne bus 2 etc.
    # L'admittance est calculée comme suit : si z est la branche entre i et j, Y_ij = inv(z.br_r + im*z.br_x) + (z.g + im*z.b), où g,b sont spécifiés dans le .m- Dans la data, ils apparaissent comme g = g_to + g_fr, b = b_to + b_fr
    Y = Matrix(calc_basic_admittance_matrix(mydata))
end

function Susceptance(mydata) # Retourne qqch de différent de Matrix(calc_basic_susceptance_matrix(mydata)). Cependant, la diff est seulement sur les termes diagonaux, et sur les i,j tq la branche ij a une valeur ratio != 0. La différence dépend d'ailleurs de ce ratio
    B = imag(Matrix(calc_basic_admittance_matrix(mydata)))
end

function Conductance(mydata) 
    G = real(Matrix(calc_basic_admittance_matrix(mydata)))
end

function LoadGens_Admittance(mydata) # La matrice est ordonnée en commençant par les loads aux lignes & colonnes 1,...,m puis les générateurs en m+1,...,n
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    Y = Admittance(mydata)
    V=[]

    for b in loads(buses)
        push!(V,b[2]["bus_i"])
    end
    for b in generators(buses)
        push!(V,b[2]["bus_i"])
    end
    T = permute_matrix(Y,V)
end

function Jacobian(mydata) # First n rows = (P1,...,Pn) ; last n rows = (Q1,...,Qn) ; first n columns = theta ; last n columns = V
    J = Matrix(calc_basic_jacobian_matrix(mydata))
end

function reduced_Jacobian(mydata)
    J = Jacobian(mydata)
    N = size(J,1)
    n = convert(Int64,N/2)
    if n != N/2
        println("How could the Jacobian be of dimension $N x $N ? An odd number of rows and columns is weird")
    end
    J1 = zeros(n,n)
    J2 = zeros(n,n)
    J3 = zeros(n,n)
    J4 = zeros(n,n)

    for i = 1:n
        for j = 1:n
            J1[i,j] = J[i,j]
            J2[i,j] = J[i,j+n]
            J3[i,j] = J[i+n,j]
            J4[i,j] = J[i+n,j+n]
        end
    end

    J1 = inv(J1)

    return J4 - J3*J1*J2

end


function Fmatrix(mydata) 
    buses = sort(collect(mydata["bus"]), by = x -> x[2]["bus_i"])
    Y = LoadGens_Admittance(mydata)
    m = length(loads(buses))
    k = length(generators(buses))
    Y_LL = Array{ComplexF64}(undef, m, m)
    Y_LG = Array{ComplexF64}(undef, m, k)

    for i=1:m
        for j=1:m
            Y_LL[i,j]=Y[i,j]
        end
    end
    for i=1:m
        for j=1:k
            Y_LG[i,j]=Y[i,m+j]
        end
    end

   if (det(Y_LL) != 0)
        return -inv(Y_LL)*Y_LG
   else
        println("The determinant of Y_LL is zero, so the whole operation is compromised.")
   end
end

function permute_matrix(Y, V) # Permutes the rows and columns of matrix Y (ie keeps it an adjacency matrix) following the permutation sigma st sigma(1,...,n) = V.
    n = size(Y, 1)
    permuted_matrix = similar(Y)

    for i in 1:n
        for j in 1:n
            permuted_matrix[i, j] = Y[V[i], V[j]]
        end
    end

    return permuted_matrix
end


