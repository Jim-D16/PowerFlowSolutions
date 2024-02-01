using Plots, Colors

function main()
V = [1, 2, 3]

for i = 1:5
    push!(V, 0)
    if (n=length(V)) > 5
        println("dingz $n")
    end
    println(n)
end

end
main()
