using CairoMakie

#figure = Makie.arrows([(1.0, 2.0), (3.0, 4.0)], [(1.0, 2.0), (3.0, 4.0)])
figure = Makie.arrows([1.0, 2.0], [3.0, 4.0], [1.0, 2.0], [3.0, 4.0])
save("./arrows_test.png", figure)
