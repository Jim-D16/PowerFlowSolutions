using Plots


function main()
    m = [1]
    n = [3]
    o = [0]
    plot3d(m, n, o, marker = :circle, color = :red, show = true)
    readline()

end
main()