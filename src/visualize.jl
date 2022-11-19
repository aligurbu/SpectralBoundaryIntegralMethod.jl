

function visualizeGeometry(Xi; Color=:red)
    fig = Figure(resolution = (1200, 800), fontsize = 22)
    ax = Axis3(fig[1, 1], aspect = :data,
            #xgridvisible = false, ygridvisible = false, zgridvisible = false,
            xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
            xlabelvisible = false, ylabelvisible = false, zlabelvisible = false,
            xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false,
            xticksvisible = false, yticksvisible = false, zticksvisible = false
            )
    sm = surface!(ax, Xi[:,:,1], Xi[:,:,2], Xi[:,:,3];
                  colormap = [Color, Color], transparency = true,
                  lightposition = Vec3f(0, 0, 0.8), ambient = Vec3f(0.6, 0.6, 0.6),
                  backlight = 1.0f0)
    wireframe!(ax, Xi[:,:,1], Xi[:,:,2], Xi[:,:,3];
               overdraw = true, transparency = true, color = (:black, 0.7))
    scatter!(ax, Xi[:,:,1][:], Xi[:,:,2][:], Xi[:,:,3][:];
             color = :black, transparency = true, overdraw = true)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    fig
end