using DrWatson
quickactivate(@__DIR__, "SpectralBoundaryIntegralMethod")
using SpectralBoundaryIntegralMethod

# This is a file to test coverage.

## see srcdir("utilities.jl")
hat([1; 2; 3])

N = 4
nlat, nlon, thet, phi, weight = gridOnSphere(N)

Xi = ellipsoidalGeometry(thet, phi, 1, 2, 0.5); # see srcdir("geometries.jl")
Xi[:,:,1] # x-components
Xi[:,:,2] # y-components
Xi[:,:,3] # z-components

# Pnm[n, m](θ)
Pnm = associatedLegendreFun(N, thet[1])
# Pnm[n, m](π)
Pnm = associatedLegendreFun(N, π)
# Pnm[m][θ, n]
Pnm = associatedLegendreFun(N, thet);
Pnm[1]
Pnm[2]
Pnm[3]
Pnm[4]
Pnm[5]


# Spherical harmonic basis functions
Pnm, Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi = 
                                        sphericalHarmonicBasisFun(N, thet, phi);
Pnm # Pnm[m][θ, n]
Pn0 # = Pnm[1]
Pnm_cos_m_phi[1] ≈ Pnm[2][:,2].*cos.((2-1)*phi)
Pnm_cos_m_phi[2] ≈ Pnm[2][:,3].*cos.((2-1)*phi)
Pnm_cos_m_phi[3] ≈ Pnm[2][:,4].*cos.((2-1)*phi)
Pnm_cos_m_phi[4] ≈ Pnm[2][:,5].*cos.((2-1)*phi)
Pnm_cos_m_phi[5] ≈ Pnm[3][:,3].*cos.((3-1)*phi)
Pnm_cos_m_phi[6] ≈ Pnm[3][:,4].*cos.((3-1)*phi)
Pnm_cos_m_phi[7] ≈ Pnm[3][:,5].*cos.((3-1)*phi)
Pnm_cos_m_phi[8] ≈ Pnm[4][:,4].*cos.((4-1)*phi)
Pnm_cos_m_phi[9] ≈ Pnm[4][:,5].*cos.((4-1)*phi)
Pnm_cos_m_phi[10] ≈ Pnm[5][:,5].*cos.((5-1)*phi)
# Similarly
Pnm_sin_m_phi[10] ≈ Pnm[5][:,5].*sin.((5-1)*phi)

Pnm, Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi,
     Pn0_wg, Pnm_cos_m_phi_wg, Pnm_sin_m_phi_wg =
                                sphericalHarmonicBasisFun(N, thet, phi, weight);
Pn0_wg # = Pn0.*weight
Pnm_cos_m_phi_wg # = Pnm_cos_m_phi.*weight
Pnm_sin_m_phi_wg # = Pnm_sin_m_phi.*weight


cXi = sphericalHarmonicAnalysis(Xi, nlat, nlon,
                                Pn0_wg, Pnm_cos_m_phi_wg, Pnm_sin_m_phi_wg);
cXi # = [aXi[:]; bXi[:]] where aXi = [(N+1)*(N+2)/2), 3] and bXi = [N*(N+1)/2), 3]
# aXi[:, 1] and bXi[:, 1] -> x-component
# aXi[:, 2] and bXi[:, 2] -> y-component
# aXi[:, 3] and bXi[:, 3] -> z-component


# or to avoid creating extra memory allocation for the coefficients cXi
cXi = zeros((N+1)^2,3)
cXi = sphericalHarmonicAnalysis!(cXi, Xi, nlat, nlon,
                                 Pn0_wg, Pnm_cos_m_phi_wg, Pnm_sin_m_phi_wg)

sphericalHarmonicSynthesis!(Xi, cXi, nlat, nlon,
                            Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi)

cXi = zeros((N+1)^2,3)

# see srcdir("geometries.jl")
N = 16
nlat, nlon, thet, phi, weight = gridOnSphere(N)
Radius = 5
Xi = sphereGeometry(thet, phi, Radius)
Xi = sphereGeometry(thet, phi, Radius;
                       Position = [-0.5; 0.5; 0],
                       OrientVec = [0; 1; 0])
visualizeGeometry(Xi)

N = 16;
nlat, nlon, thet, phi, weight = gridOnSphere(N);
function testTime(kk)
    for _ in 1:kk
        Pn = associatedLegendreFun(N, thet)
    end
end
@time testTime(10000)


Pnm, Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi,
     Pn0_wg, Pnm_cos_m_phi_wg, Pnm_sin_m_phi_wg =
                                sphericalHarmonicBasisFun(N, thet, phi, weight);

Xi = ellipsoidalGeometry(thet, phi, 1, 0.5, 2);

thet_equal = [range(0,pi,nlat);];
Eq_Pnm, Eq_Pn0, Eq_Pnm_cos_m_phi, Eq_Pnm_sin_m_phi =
                                  sphericalHarmonicBasisFun(N, thet_equal, phi);

cXi = zeros((N+1)^2,3);
cXi = sphericalHarmonicAnalysis!(cXi, Xi, nlat, nlon, Pn0_wg,
                                 Pnm_cos_m_phi_wg, Pnm_sin_m_phi_wg);

cG = cXi;
G = zeros(size(Xi));
G = sphericalHarmonicSynthesis!(G, cG, nlat, nlon,
                                     Pn0, Pnm_cos_m_phi, Pnm_sin_m_phi);

maximum(abs.(G - Xi))

visualizeGeometry(Xi, cXi, nlat, nlon,
                  Eq_Pn0, Eq_Pnm_cos_m_phi, Eq_Pnm_sin_m_phi;
                  Color=:red)
