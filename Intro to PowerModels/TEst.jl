using Plots, Colors

mycolors = cgrad([:orange, :blue], 30, categorical = true)
println(mycolors[3])



L = collect(LinRange(0.025, 0.3, 12))
for i = 1:100
    plot!([i], [i], marker = :circle, color = mycolors[i], show = true, label = "")
    
end
readline()