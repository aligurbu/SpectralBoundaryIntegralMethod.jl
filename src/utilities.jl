

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

"""
    Associated Legendre function
    Pnm[n, m](θ)
"""
function associatedLegendreFun(n::Int64, θ::Float64)
    CondonShortleyPhase = (-1).^[0:n;]'
    legendre(LegendreOrthoNorm(), 0:n, 0:n, cos(θ)).*CondonShortleyPhase
end
"""
    Associated Legendre function
    Pnm[n, m](π)
"""
function associatedLegendreFun(n::Int64, θ::Irrational{:π})
    θ = convert(Float64,θ)
    associatedLegendreFun(n, θ)
end
"""
    Associated Legendre function
    Pnm[m][θ, n]
"""
function associatedLegendreFun(n::Int64, θ::Vector{Float64})
    Pnm = legendre(LegendreOrthoNorm(), 0:n, 0:n, cos.(θ))
    [Pnm[:,:,m+1].*(-1)^m for m in 0:n]
end

"""
    Spherical harmonic basis functions and weighted basis functions
"""
function sphericalHarmonicBasisFun(nlat, nlon, thet, phi, weight)
    # Trigonometric functions
    cos_m_phi = [cos.((m-1)*phi) for m=2:nlat]
    sin_m_phi = [sin.((m-1)*phi) for m=2:nlat]

    Pnm = associatedLegendreFun(N, thet)
    Pn0 = Pnm[1]
    Pn0_wg = Pn0.*weight
    Pnm_cos_m_phi = [Pnm[m][:,n].*cos_m_phi[m-1] for m=2:N+1 for n=m:N+1]
    Pnm_sin_m_phi = [Pnm[m][:,n].*sin_m_phi[m-1] for m=2:N+1 for n=m:N+1]
    Pnm_cos_m_phi_wg = [Pnm_cos_m_phi[k].*weight for k = 1:Int(N*(N+1)/2)]
    Pnm_sin_m_phi_wg = [Pnm_sin_m_phi[k].*weight for k = 1:Int(N*(N+1)/2)]
    return Pnm, Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi,
                Pn0_wg, Pnm_cos_m_phi_wg, Pnm_sin_m_phi_wg
end