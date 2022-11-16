using DrWatson
quickactivate(@__DIR__, "SpectralBoundaryIntegralMethod")
using SpectralBoundaryIntegralMethod

# This is a file to test coverage.
N = 16
nlat, nlon, thet, phi, weight = gridOnSphere(N)
Radius = 5
Xi = getSphereGeometry(thet, phi, Radius)
Xi = getSphereGeometry(thet, phi, Radius;
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