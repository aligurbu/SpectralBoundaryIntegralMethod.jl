

"""
    The first fundamental form coefficients
"""
function firstFundamentalFormCoeff(grad_thet_G, grad_phi_G)
    EE = sum(grad_thet_G.*grad_thet_G, dims = 3)
    FF = sum(grad_thet_G.*grad_phi_G, dims = 3)
    GG = sum(grad_phi_G.*grad_phi_G, dims = 3)
    WW = EE.*GG - FF.^2
    Jbrev = sqrt.(WW) # Jacobian determinant
    return EE, FF, GG, WW, Jbrev
end

"""
    The unit normal to the undeformed surface
"""
function unitNormalVector(grad_thet_G, grad_phi_G, JGbrev)
    unitNormal = zeros(size(grad_thet_G))
    unitNormal[:,:,1] = (grad_thet_G[:,:,2].*grad_phi_G[:,:,3] -
                           grad_thet_G[:,:,3].*grad_phi_G[:,:,2])./JGbrev
    unitNormal[:,:,2] = (grad_thet_G[:,:,3].*grad_phi_G[:,:,1] -
                           grad_thet_G[:,:,1].*grad_phi_G[:,:,3])./JGbrev
    unitNormal[:,:,3] = (grad_thet_G[:,:,1].*grad_phi_G[:,:,2] -
                           grad_thet_G[:,:,2].*grad_phi_G[:,:,1])./JGbrev
    return unitNormal
end