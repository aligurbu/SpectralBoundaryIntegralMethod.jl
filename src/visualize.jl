

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

"""
   This encloses gap in the north and south poles,
   as well as the gap between the first and last longitudinal points.
"""
function visualizeGeometry(Xi, cXi, nlat, nlon,
                           Eq_Pn0, Eq_Pnm_cos_m_phi, Eq_Pnm_sin_m_phi;
                           Color=:red)
    Eq_Xi = zeros(nlat, nlon, size(Xi,3));
    Eq_Xi = sphericalHarmonicSynthesis!(Eq_Xi, cXi, nlat, nlon,
                        Eq_Pn0, Eq_Pnm_cos_m_phi, Eq_Pnm_sin_m_phi)

    Xi_poles = [reshape(Eq_Xi[1,:,:], 1, nlon, 3); Xi; reshape(Eq_Xi[end,:,:],1,nlon,3)]

    Xiplot = zeros(size(Xi,1)+2, size(Xi,2)+1, size(Xi,3))
    [Xiplot[:,:,kk] = hcat(Xi_poles[:,:,kk], Xi_poles[:,1,kk]) for kk = 1:3]

    fig = Figure(resolution = (1200, 800), fontsize = 22)
    ax = Axis3(fig[1, 1], aspect = :data,
        #xgridvisible = false, ygridvisible = false, zgridvisible = false,
        xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
        xlabelvisible = false, ylabelvisible = false, zlabelvisible = false,
        xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false,
        xticksvisible = false, yticksvisible = false, zticksvisible = false
        )
    sm = surface!(ax, Xiplot[:,:,1], Xiplot[:,:,2], Xiplot[:,:,3];
                  colormap = [Color, Color], transparency = true,
                  lightposition = Vec3f(0, 0, 0.8), ambient = Vec3f(0.6, 0.6, 0.6),
                  backlight = 1.0f0)
    wireframe!(ax, Xiplot[:,:,1], Xiplot[:,:,2], Xiplot[:,:,3];
               overdraw = true, transparency = true, color = (:black, 0.7))
    scatter!(ax, Xiplot[:,:,1][:], Xiplot[:,:,2][:], Xiplot[:,:,3][:];
             color = :black, transparency = true, overdraw = true, markersize = 10)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    fig
end