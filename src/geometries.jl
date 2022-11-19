
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
    Deformed spherical geometry = twist about x axis
"""
function sphereGeometry(thet,phi,Radius,twist)
    # Initial undeformed spherical shape
    Xi = sphereGeometry(thet, phi, Radius)
    # Twisting
    xi = zeros(size(Xi))
    for m = 1:size(Xi,1)
        for n = 1:size(Xi,2)
            XX = Xi[m,n,:]
            rotvec = [0 0 0;
                      0 0 -twist;
                      0 twist 0]*XX[1]
            xi[m,n,:] = exp(rotvec)*XX
        end
    end
    return xi
end

"""
    Undeformed shape of red blood cell (RBC) taken from
    Pozrikidis, C. (2005). "Axisymmetric motion of a file of red blood cells
    through capillaries." Physics of Fluids, 17(3) [Equation (31)].
    Position    = Position of center of mass of RBC
    OrientVec = Orientation around 'OrientVec' of RBC
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

"""
    Deformed shape red blood cell shape
    Bending about y-axis followed by twist about x axis
    rho = radius of curvature
    twist = twist
"""
function RBCDeformedGeometry(thet,phi,rho,twist)

    # Initial undeformed shape of RBC
    Xi = RBCInitialGeometry(thet, phi)
    alpha = 1.3858189

    xi = zeros(size(Xi))
    for m = 1:size(Xi,1)
        for n = 1:size(Xi,2)
            XX = Xi[m,n,:]
            vartheta = XX[1]/rho
            xx = [(rho-XX[3])*sin(vartheta);
                  XX[2];
                  rho-(rho-XX[3])*cos(vartheta)]
            rotvec = [0 0 0;
                      0 0 -twist;
                      0 twist 0]*(XX[1]/alpha)
            xi[m,n,:] = exp(rotvec)*xx
        end
    end
    return xi
end

"""
    Construct the ellipsoid geometry on the given (thet,phi) grid
"""
function ellipsoidalGeometry(thet, phi, a, b, c;
                             Position = [0; 0; 0],
                             OrientVec = [0; 0; 0])

    nlat = length(thet)
    nlon = length(phi)
    Xi = zeros(nlat, nlon, 3)

    Xi[:,:,1] = a*sin.(thet).*cos.(phi)
    Xi[:,:,2] = b*sin.(thet).*sin.(phi)
    Xi[:,:,3] = c*cos.(thet).*ones(1,nlon)

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
    Egg shape taken from http://www.mathematische-basteleien.de/eggcurves.htm
    from the section "From oval to the egg shape"
    z^2 + x^2*t(z) = 1;
    t(z) = 1 + bet*z
    becomes the following when parametrized in spherical coordinates
"""
function eggGeometry(thet, phi;
                     Position = [0; 0; 0],
                     OrientVec = [0; 0; 0])

    nlat = length(thet)
    nlon = length(phi)
    Xi = zeros(nlat, nlon, 3)

    alph = 0.2
    Xi[:,:,1] = sin.(thet).*cos.(phi)
    Xi[:,:,2] = sin.(thet).*sin.(phi)
    Xi[:,:,3] = (1 .+ alph*cos.(thet)).*cos.(thet).*ones(1,nlon)

    Orientation = exp(hat(OrientVec*pi/2))
    for m = 1:nlat
        for n = 1:nlon
            x = Position + Orientation*Xi[m,n,:]
            Xi[m,n,:] = x
        end
    end
    return Xi
end