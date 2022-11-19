
"""
    Construct a sphere geometry on the given (thet, phi) grid and Radius

"""
function sphereGeometry(thet, phi, Radius)
    nlat = length(thet)
    nlon = length(phi)
    Xi = zeros(nlat, nlon, 3)
    Xi[:,:,1] = Radius*sin.(thet).*cos.(phi)
    Xi[:,:,2] = Radius*sin.(thet).*sin.(phi)
    Xi[:,:,3] = Radius*cos.(thet).*ones(1,nlon)
    return Xi
end

"""
    Construct a sphere geometry on the given (thet, phi) grid and Radius

    Position    = Position of center of mass
    OrientVec = Orientation around 'OrientVec'
"""
function sphereGeometry(thet, phi, Radius;
                        Position = [0; 0; 0],
                        OrientVec = [0; 0; 0])

    nlat = length(thet)
    nlon = length(phi)
    Xi = zeros(nlat, nlon, 3)
    Xi[:,:,1] = Radius*sin.(thet).*cos.(phi)
    Xi[:,:,2] = Radius*sin.(thet).*sin.(phi)
    Xi[:,:,3] = Radius*cos.(thet).*ones(1,nlon)

    Orientation = exp(hat(OrientVec*pi/2))
    for m = 1:nlat
        for n = 1:nlon
            x = Position + Orientation*Xi[m,n,:]
            Xi[m,n,:] = x
        end
    end
    return Xi
end

"""
    Undeformed shape of red blood cell (RBC) taken from
    Pozrikidis, C. (2005). "Axisymmetric motion of a file of red blood cells
    through capillaries." Physics of Fluids, 17(3) [Equation (31)].
    Position    = Position of center of mass of RBC
    OrientVec = Orientation around 'OrientVec' of RBC
    (Position,InitOrient) give the initial position and orientation of the RBC
"""
function RBCInitialGeometry(thet, phi;
                            Position = [0; 0; 0],
                            OrientVec = [0; 0; 0])

    nlat = length(thet)
    nlon = length(phi)
    alpha = 1.3858189

    Xi = zeros(nlat, nlon, 3)
    Xi[:,:,1] = alpha*(sin.(thet).*cos.(phi))
    Xi[:,:,2] = alpha*(sin.(thet).*sin.(phi))
    Xi[:,:,3] = (0.5*alpha)*
                (0.207 .+ 2.003*sin.(thet).^2 .- 1.123*sin.(thet).^4).*
                                                    cos.(thet).*ones(1,nlon)

    Orientation = exp(hat(OrientVec*pi/2))
    for m = 1:nlat
        for n = 1:nlon
            x = Position + Orientation*Xi[m,n,:]
            Xi[m,n,:] = x
        end
    end
    return Xi
end