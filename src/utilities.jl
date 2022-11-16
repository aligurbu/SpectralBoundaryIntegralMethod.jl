

"""
    Produce a skew-symmetric matrix using the components of x vector
"""
function hat(x::Vector)
    X = [   0  -x[3]  x[2]
          x[3]    0  -x[1]
         -x[2]  x[1]    0 ]
end

"""
    The quadrature grid points on the unit sphere
"""
function gridOnSphere(N)
    nlat, nlon = N+1, 2N+1
    thet_, weight = gausslegendre(nlat)
    thet = acos.(thet_)
    reverse!(thet)
    phi = transpose([(0:1/nlon:(1-1/nlon))*2π;])
    return nlat, nlon, thet, phi, weight
end
