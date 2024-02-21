using Plots, Colors



function ShowTheCycle(L, angles)

# Define polygon vertices and angles
n = L # Number of vertices
vertices = [(cos(2π*i/n), sin(2π*i/n)) for i in 0:n-1]

# Calculate endpoints of vectors
vectors = [(0.4*cos(angle), 0.4*sin(angle)) for angle in angles]


# Plot polygon and vectors
plot(vertices, seriestype = :shape, fillalpha = 0, label="", aspect_ratio = :equal)
for i in 1:n
    quiver!([vertices[i]], quiver=([vectors[i]]), color=:red, label="", show = true)
end
readline()
end
main()
