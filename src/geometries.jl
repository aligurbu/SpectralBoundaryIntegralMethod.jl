
"""
    Construct a sphere geometry on the given (thet, phi) grid and Radius

"""
function getSphereGeometry(thet, phi, Radius)
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
function getSphereGeometry(thet, phi, Radius;
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
