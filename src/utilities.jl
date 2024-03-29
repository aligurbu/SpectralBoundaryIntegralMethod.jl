"""
    hat(x::Vector)

Produce a skew-symmetric matrix using the components of x vector
"""
function hat(x::Vector)
    X = [   0  -x[3]  x[2]
          x[3]    0  -x[1]
         -x[2]  x[1]    0 ]
end

"""
    etaFun(N, I_thet)

The regularization function for the integrals on a unit sphere
"""
function etaFun(N, I_thet)
    I_eta = 2*sin.(I_thet/2).*sum(Plm(0:N, 0, cos.(I_thet)), dims = 2)
end

"""
    upSampling(cG, N, N_up)

Up-sampling the coefficients of the spherical harmonics expansions
"""
function upSampling(cG, N, N_up)
    lengthofCoeffVector = Int(N*(N+1)/2)
    an0 = view(cG, 1:N+1, :)
    anm = view(cG, N+1 .+ (1:lengthofCoeffVector), :)
    bnm = view(cG, N+1 + lengthofCoeffVector .+ (1:lengthofCoeffVector), :)

    up_lengthofCoeffVector = Int(N_up*(N_up+1)/2)
    up_cG = zeros((N_up+1)^2,size(cG,2))
    up_cG[1:N+1, :] = an0
    up_cG[N_up+1 .+ (1:lengthofCoeffVector), :] = anm
    up_cG[N_up+1 + up_lengthofCoeffVector .+ (1:lengthofCoeffVector), :] = bnm
    return up_cG
end

"""
    gridOnSphere(N)

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
    integrationGridOnSphere(N, NGSphere)

The quadrature grid points on the unit sphere for the integration of the
weakly boundary integrals equations.
Also computing the relevant quantities on this quadratures such as the
spherical harmonic basis functions and the derivatives of the associated
Legendre functions.
"""
function integrationGridOnSphere(N, NGSphere)
    I_nlat, I_nlon, I_thet, I_phi, I_weight = gridOnSphere(NGSphere);

    I_eta = etaFun(N, I_thet)

    I_Pnm, I_Pn0, I_Pnm_cos_m_phi, I_Pnm_sin_m_phi =
                    sphericalHarmonicBasisFun(I_nlat, I_nlon, I_thet, I_phi);

    I_DPnm, I_D2Pnm = derivativesAssociatedLegendreFun(I_Pnm);

    return I_nlat, I_nlon, I_thet, I_phi, I_weight, I_eta,
           I_Pnm, I_Pn0, I_Pnm_cos_m_phi, I_Pnm_sin_m_phi,
           I_DPnm, I_D2Pnm
end


"""
    associatedLegendreFun(n::Int64, θ::Float64)

Associated Legendre function
Pnm[n, m](θ)
"""
function associatedLegendreFun(n::Int64, θ::Float64)
    CondonShortleyPhase = (-1).^[0:n;]'
    legendre(LegendreOrthoNorm(), 0:n, 0:n, cos(θ)).*CondonShortleyPhase
end
"""
    associatedLegendreFun(n::Int64, θ::Irrational{:π})

Associated Legendre function
Pnm[n, m](π)
"""
function associatedLegendreFun(n::Int64, θ::Irrational{:π})
    θ = convert(Float64,θ)
    associatedLegendreFun(n, θ)
end
"""
    associatedLegendreFun(n::Int64, θ::Vector{Float64})

Associated Legendre function
Pnm[m][θ, n]
"""
function associatedLegendreFun(n::Int64, θ::Vector{Float64})
    Pnm = legendre(LegendreOrthoNorm(), 0:n, 0:n, cos.(θ))
    [Pnm[:,:,m+1].*(-1)^m for m in 0:n]
end

"""
    sphericalHarmonicBasisFun(N, thet, phi, weight)

Spherical harmonic basis functions and weighted basis functions
"""
function sphericalHarmonicBasisFun(N, thet, phi, weight)
    # Trigonometric functions
    cos_m_phi = [cos.((m-1)*phi) for m=2:N+1]
    sin_m_phi = [sin.((m-1)*phi) for m=2:N+1]

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
"""
    sphericalHarmonicBasisFun(N, thet, phi)

Spherical harmonic basis functions
"""
function sphericalHarmonicBasisFun(N, thet, phi)
    # Trigonometric functions
    cos_m_phi = [cos.((m-1)*phi) for m=2:N+1]
    sin_m_phi = [sin.((m-1)*phi) for m=2:N+1]

    Pnm = associatedLegendreFun(N, thet)
    Pn0 = Pnm[1]
    Pnm_cos_m_phi = [Pnm[m][:,n].*cos_m_phi[m-1] for m=2:N+1 for n=m:N+1]
    Pnm_sin_m_phi = [Pnm[m][:,n].*sin_m_phi[m-1] for m=2:N+1 for n=m:N+1]
    return Pnm, Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi
end

"""
    sphericalHarmonicAnalysis(G::Array{Float64, 3},
                              nlat::Int64, nlon::Int64,
                              Pn0_wg::Matrix{Float64},
                              Pnm_cos_m_phi_wg::Vector{Matrix{Float64}},
                              Pnm_sin_m_phi_wg::Vector{Matrix{Float64}})

Spherical harmonic analysis on the components of the position, velocity, and
force fields.
The analysis is performed on a gaussian grid in latitude and an equally
spaced grid in longitude.
The spherical harmonics basis functions are precomputed and stored.

In this version of the analysis function, the coefficients are stored in the format we used in Matlab.
"""
function sphericalHarmonicAnalysis(G::Array{Float64, 3},
                                   nlat::Int64, nlon::Int64,
                                   Pn0_wg::Matrix{Float64},
                                   Pnm_cos_m_phi_wg::Vector{Matrix{Float64}},
                                   Pnm_sin_m_phi_wg::Vector{Matrix{Float64}})
    N = nlat - 1
    numDimension = size(G,3)
    aG = zeros(Int((N+1)*(N+2)/2),numDimension)
    bG = zeros(Int(N*(N+1)/2),numDimension)
    k, ℓ = 1, 1
    for n=1:nlat
        aG[k,:] = sum(G.*Pn0_wg[:,n], dims =[1 2])
        k = k + 1
        for m=2:n
            aG[k,:] = sum(G.*Pnm_cos_m_phi_wg[ℓ], dims = [1 2])
            bG[ℓ,:] = sum(G.*Pnm_sin_m_phi_wg[ℓ], dims = [1 2])
            k, ℓ = k + 1, ℓ + 1
        end
    end
    aG = (2/nlon)*aG
    bG = -(2/nlon)*bG
    cG = [aG[:]; bG[:]]
    return cG
end

"""
    sphericalHarmonicAnalysis!(cG::Matrix{Float64}, G::Array{Float64, 3},
                               nlat::Int64, nlon::Int64,
                               Pn0_wg::Matrix{Float64},
                               Pnm_cos_m_phi_wg::Vector{Matrix{Float64}},
                               Pnm_sin_m_phi_wg::Vector{Matrix{Float64}})

Spherical harmonic analysis on the components of the position, velocity, and
force fields.
The analysis is performed on a gaussian grid in latitude grid and an equally
spaced grid in longitude.
The spherical harmonics basis functions are precomputed and stored.
"""
function sphericalHarmonicAnalysis!(cG::Matrix{Float64}, G::Array{Float64, 3},
                                    nlat::Int64, nlon::Int64,
                                    Pn0_wg::Matrix{Float64},
                                    Pnm_cos_m_phi_wg::Vector{Matrix{Float64}},
                                    Pnm_sin_m_phi_wg::Vector{Matrix{Float64}})
    N = nlat - 1
    lengthofCoeffVector = Int(N*(N+1)/2)
    an0 = view(cG, 1:nlat, :)
    anm = view(cG, nlat .+ (1:lengthofCoeffVector), :)
    bnm = view(cG, nlat + lengthofCoeffVector .+ (1:lengthofCoeffVector), :)
    for n=1:nlat
        an0[n,:] = (2/nlon)*sum(G.*Pn0_wg[:,n], dims =[1 2])
    end
    for ℓ =1:lengthofCoeffVector
        anm[ℓ,:] = (2/nlon)*sum(G.*Pnm_cos_m_phi_wg[ℓ], dims = [1 2])
        bnm[ℓ,:] = -(2/nlon)*sum(G.*Pnm_sin_m_phi_wg[ℓ], dims = [1 2])
    end
    return cG
end


"""
    sphericalHarmonicSynthesis!(G::Array{Float64, 3},
                                cG::Matrix{Float64},
                                nlat::Int64, nlon::Int64,
                                Pn0::Matrix{Float64},
                                Pnm_cos_m_phi::Vector{Matrix{Float64}},
                                Pnm_sin_m_phi::Vector{Matrix{Float64}})

Spherical harmonic synthesis on the spherical harmonics coefficients.
The synthesis is performed on a gaussian latitude grid and an equally
spaced grid longitude.
The spherical harmonic basis functions are precomputed and stored.
"""
function sphericalHarmonicSynthesis!(G::Array{Float64, 3},
                                     cG::Matrix{Float64},
                                     nlat::Int64, nlon::Int64,
                                     Pn0::Matrix{Float64},
                                     Pnm_cos_m_phi::Vector{Matrix{Float64}},
                                     Pnm_sin_m_phi::Vector{Matrix{Float64}})
    N = nlat - 1
    lengthofCoeffVector = Int(N*(N+1)/2)::Int64
    dimension = size(cG,2)
    an0 = view(cG, 1:nlat, :)
    anm = view(cG, nlat .+ (1:lengthofCoeffVector), :)
    bnm = view(cG, nlat + lengthofCoeffVector .+ (1:lengthofCoeffVector), :)

    G += 0.5*sum(reshape(an0,1,nlat,dimension).*Pn0, dims = 2).*ones(1,nlon)
    for kk = 1:lengthofCoeffVector
        G += reshape(anm[kk,:], 1, 1, dimension).*Pnm_cos_m_phi[kk]
        G -= reshape(bnm[kk,:], 1, 1, dimension).*Pnm_sin_m_phi[kk]
    end
    return G
end

"""
    derivativesAssociatedLegendreFun(Pnm::Vector{Matrix{Float64}})

"""
function derivativesAssociatedLegendreFun(Pnm::Vector{Matrix{Float64}})
    nlat = size(Pnm,1)
    DPnm = [zeros(nlat, nlat) for k=1:nlat]
    D2Pnm = [zeros(nlat, nlat) for k=1:nlat]

    # Derivative of associated legendre function wrt theta
    # DP00 = 0
    for n = 1:nlat-1
        # m == 0
        # DPn0 = - sqrt(n*(n+1))*Pn1
        DPnm[1][:,n+1] = - sqrt(n*(n+1))*Pnm[2][:,n+1]
        for m = 1:n
            if m == n
                # DPnn = sqrt(n/2)*Pn,n-1
                DPnm[n+1][:,n+1] = sqrt(n/2)*Pnm[n][:,n+1]
            else
                # DPnm = 0.5*sqrt((n+m)*(n-m+1))Pn,m-1 -
                #        0.5*sqrt((n+1+m)*(n-m))*Pn,m+1
                fac1 = 0.5*sqrt((n+m)*(n-m+1))
                fac2 = 0.5*sqrt((n+1+m)*(n-m))
                DPnm[m+1][:,n+1] = fac1*Pnm[m][:,n+1] - fac2*Pnm[m+2][:,n+1]
            end
        end
    end

    # Second derivative of associated legendre function wrt theta
    # D2P00 = 0
    for n = 1:nlat-1
        # m == 0
        if n == 1
            # D2P10 = - P10
            D2Pnm[1][:,n+1] = - Pnm[1][:,n+1]
        else
            # D2Pn0 = -0.5*(n+1)*n*Pn0 + 0.5*sqrt((n-1)*n*(n+1)*(n+2))*Pn2
            D2Pnm[1][:,n+1] = -0.5*(n+1)*n*Pnm[1][:,n+1] +
                               0.5*sqrt((n-1)*n*(n+1)*(n+2))*Pnm[3][:,n+1]
        end
        # m == 1
        if n == 1
            # D2P11 = -P11
            D2Pnm[2][:,n+1] = -Pnm[2][:,n+1]
        elseif n == 2
            # D2P21 = - 4*P21
            D2Pnm[2][:,n+1] = -4*Pnm[2][:,n+1]
        else
            # D2Pn1 = -0.5*(n*(n+1) + 0.5*(n+2)*(n-1))*Pn1
            #           + 0.25*sqrt((n-2)*(n-1)*(n+2)*(n+3))*Pn3
            D2Pnm[2][:,n+1] = -0.5*(n*(n+1) + 0.5*(n+2)*(n-1))*Pnm[2][:,n+1] +
                               0.25*sqrt((n-2)*(n-1)*(n+2)*(n+3))*Pnm[4][:,n+1]
        end
        for m = 2:n
            if m == n
                # D2Pnn = 0.25*sqrt(4*n*(2*n-1))*Pn,n-2 - 0.5*n*Pnn
                D2Pnm[n+1][:,n+1] = 0.25*sqrt(4*n*(2*n-1))*Pnm[n-1][:,n+1] -
                                    0.5*n*Pnm[n+1][:,n+1]
            elseif m == n-1
                # D2Pn,n-1 = 0.5*sqrt(3*(n-1)*(2*n-1))*Pn,n-3 - 0.5*(3*n-1)*Pn,n-1
                D2Pnm[n][:,n+1] = 0.5*sqrt(3*(n-1)*(2*n-1))*Pnm[n-2][:,n+1] -
                                  0.5*(3*n-1)*Pnm[n][:,n+1]
            else
                # D2Pnm = 0.25*sqrt( (n+m)*(n+m-1)*(n-m+1)*(n-m+2) )*Pn,m-2
                #         -0.25*( (n+m)*(n-m+1) + (n+m+1)*(n-m) )*Pnm
                #         +0.25*sqrt((n-m)*(n-m-1)*(n+m+1)*(n+m+2))*Pn,m+2
                D2Pnm[m+1][:,n+1] =
                    0.25*sqrt((n+m)*(n+m-1)*(n-m+1)*(n-m+2))*Pnm[m-1][:,n+1] -
                    0.25*((n+m)*(n-m+1) + (n+m+1)*(n-m))*Pnm[m+1][:,n+1] +
                    0.25*sqrt((n-m)*(n-m-1)*(n+m+1)*(n+m+2))*Pnm[m+3][:,n+1]
            end
        end
    end
    return DPnm, D2Pnm
end


"""
    gradient(cG, nlat, nlon, thet, phi, Pnm, DPnm, D2Pnm)

The first and second order gradient of a scalar field, G,

    ∇_θ G = ∂G/∂θ;      ∇_ϕ G = 1/sinθ ∂G/∂ϕ
    ∇_θθ G = ∂²G/∂θ²;   ∇_θϕ G = 1/sinθ ∂²G/∂θ∂ϕ;   ∇_ϕϕ G = 1/(sinθ)² ∂²G/∂ϕ²

using its spherical harmonics coefficients cG, precomputed and stored.
"""
function gradient(cG, nlat, nlon, thet, phi, Pnm, DPnm, D2Pnm)
    dimension = size(cG, 2)
    N = nlat - 1
    lengthofCoeffVector = Int(N*(N+1)/2)
    an0 = reshape(cG[1:nlat, :], 1, nlat, dimension)
    anm = cG[nlat .+ (1:lengthofCoeffVector), :]
    bnm = cG[nlat + lengthofCoeffVector .+ (1:lengthofCoeffVector), :]

    # Computing the derivative wrt theta and phi
    derivative_thet_G = zeros(nlat, nlon, dimension)
    derivative_phi_G = zeros(nlat, nlon, dimension)

    derivative_thet_thet_G = zeros(nlat, nlon, dimension)
    derivative_phi_phi_G = zeros(nlat, nlon, dimension)
    derivative_thet_phi_G = zeros(nlat, nlon, dimension)

    derivative_thet_G += 0.5*sum(an0.*DPnm[1], dims = 2).*ones(1,nlon)
    derivative_thet_thet_G += 0.5*sum(an0.*D2Pnm[1], dims = 2).*ones(1,nlon)
    kk = 1
    for n=2:nlat
        for m=2:n
            anm_k = reshape(anm[kk,:], 1, 1, dimension)
            bnm_k = reshape(bnm[kk,:], 1, 1, dimension)
            derivative_thet_G += anm_k.*DPnm[m][:,n].*cos.((m-1)*phi) -
                                 bnm_k.*DPnm[m][:,n].*sin.((m-1)*phi)
            derivative_phi_G += (m-1)*Pnm[m][:,n].*
                              (-anm_k.*sin.((m-1)*phi) - bnm_k.*cos.((m-1)*phi))
            derivative_thet_thet_G += (anm_k.*D2Pnm[m][:,n].*cos.((m-1)*phi) -
                                       bnm_k.*D2Pnm[m][:,n].*sin.((m-1)*phi))
            derivative_phi_phi_G += (m-1)^2*Pnm[m][:,n].*
                              (-anm_k.*cos.((m-1)*phi) + bnm_k.*sin.((m-1)*phi))
            derivative_thet_phi_G += (m-1)*DPnm[m][:,n].*
                              (-anm_k.*sin.((m-1)*phi) - bnm_k.*cos.((m-1)*phi))
            kk += 1
        end
    end

    grad_thet_G = derivative_thet_G
    grad_phi_G = derivative_phi_G./sin.(thet)

    grad_thet_thet_G = derivative_thet_thet_G
    grad_phi_phi_G = derivative_phi_phi_G./(sin.(thet).^2)
    grad_thet_phi_G = derivative_thet_phi_G./sin.(thet)
    return grad_thet_G, grad_phi_G,
           grad_thet_thet_G, grad_phi_phi_G, grad_thet_phi_G
end

"""
    gradient(cG, nlat, nlon, thet, phi, Pnm, DPnm)

Gradient of a scalar field, G,

    ∇_θ G = ∂G/∂θ;      ∇_ϕ G = 1/sinθ ∂G/∂ϕ

using its spherical harmonics coefficients cG, precomputed and stored.
"""
function gradient(cG, nlat, nlon, thet, phi, Pnm, DPnm)
    dimension = size(cG, 2)
    N = nlat - 1
    lengthofCoeffVector = Int(N*(N+1)/2)
    an0 = reshape(cG[1:nlat, :], 1, nlat, dimension)
    anm = cG[nlat .+ (1:lengthofCoeffVector), :]
    bnm = cG[nlat + lengthofCoeffVector .+ (1:lengthofCoeffVector), :]

    # Computing the derivative wrt theta and phi
    derivative_thet_G = zeros(nlat, nlon, dimension)
    derivative_phi_G = zeros(nlat, nlon, dimension)

    derivative_thet_G += 0.5*sum(an0.*DPnm[1], dims = 2).*ones(1,nlon)
    kk = 1
    for n=2:nlat
        for m=2:n
            anm_k = reshape(anm[kk,:], 1, 1, dimension)
            bnm_k = reshape(bnm[kk,:], 1, 1, dimension)
            derivative_thet_G += anm_k.*DPnm[m][:,n].*cos.((m-1)*phi) -
                                 bnm_k.*DPnm[m][:,n].*sin.((m-1)*phi)
            derivative_phi_G += (m-1)*Pnm[m][:,n].*
                              (-anm_k.*sin.((m-1)*phi) - bnm_k.*cos.((m-1)*phi))
            kk += 1
        end
    end

    grad_thet_G = derivative_thet_G
    grad_phi_G = derivative_phi_G./sin.(thet)
    return grad_thet_G, grad_phi_G
end