

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
    The second fundamental form coefficients
"""
function secondFundamentalFormCoeff(grad_thet_thet_G,
                                    grad_thet_phi_G,
                                    grad_phi_phi_G,
                                    unitNormalG)
    LL = -sum(grad_thet_thet_G.* unitNormalG, dims = 3)
    MM = -sum(grad_thet_phi_G.* unitNormalG, dims = 3)
    NN = -sum(grad_phi_phi_G .* unitNormalG, dims = 3)
    return LL, MM, NN
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

function derivativesUnitNormalVector(unitNormalG,
                                     grad_thet_G, grad_phi_G,
                                     grad_thet_thet_G, grad_phi_phi_G,
                                     grad_thet_phi_G,
                                     JGbrev)
    PP = zeros(size(unitNormalG, 1), size(unitNormalG, 2), 6)
    PP[:,:,1] = 1 .- unitNormalG[:,:,1].*unitNormalG[:,:,1]
    PP[:,:,2] = 1 .- unitNormalG[:,:,2].*unitNormalG[:,:,2]
    PP[:,:,3] = 1 .- unitNormalG[:,:,3].*unitNormalG[:,:,3]
    PP[:,:,4] = - unitNormalG[:,:,2].*unitNormalG[:,:,3]
    PP[:,:,5] = - unitNormalG[:,:,1].*unitNormalG[:,:,3]
    PP[:,:,6] = - unitNormalG[:,:,1].*unitNormalG[:,:,2]

    grad_thet_thet_G_X_grad_phi_G = zeros(size(unitNormalG))
    grad_thet_thet_G_X_grad_phi_G[:,:,1] =
                            (grad_thet_thet_G[:,:,2].*grad_phi_G[:,:,3] -
                             grad_thet_thet_G[:,:,3].*grad_phi_G[:,:,2])
    grad_thet_thet_G_X_grad_phi_G[:,:,2] =
                            (grad_thet_thet_G[:,:,3].*grad_phi_G[:,:,1] -
                             grad_thet_thet_G[:,:,1].*grad_phi_G[:,:,3])
    grad_thet_thet_G_X_grad_phi_G[:,:,3] =
                            (grad_thet_thet_G[:,:,1].*grad_phi_G[:,:,2] -
                             grad_thet_thet_G[:,:,2].*grad_phi_G[:,:,1])

    grad_thet_G_X_grad_thet_phi_G = zeros(size(unitNormalG))
    grad_thet_G_X_grad_thet_phi_G[:,:,1] =
                            (grad_thet_G[:,:,2].*grad_thet_phi_G[:,:,3] -
                             grad_thet_G[:,:,3].*grad_thet_phi_G[:,:,2])
    grad_thet_G_X_grad_thet_phi_G[:,:,2] =
                            (grad_thet_G[:,:,3].*grad_thet_phi_G[:,:,1] -
                             grad_thet_G[:,:,1].*grad_thet_phi_G[:,:,3])
    grad_thet_G_X_grad_thet_phi_G[:,:,3] =
                            (grad_thet_G[:,:,1].*grad_thet_phi_G[:,:,2] -
                             grad_thet_G[:,:,2].*grad_thet_phi_G[:,:,1])

    grad_thet_phi_G_X_grad_phi_G = zeros(size(unitNormalG))
    grad_thet_phi_G_X_grad_phi_G[:,:,1] =
                            (grad_thet_phi_G[:,:,2].*grad_phi_G[:,:,3] -
                             grad_thet_phi_G[:,:,3].*grad_phi_G[:,:,2])
    grad_thet_phi_G_X_grad_phi_G[:,:,2] =
                            (grad_thet_phi_G[:,:,3].*grad_phi_G[:,:,1] -
                             grad_thet_phi_G[:,:,1].*grad_phi_G[:,:,3])
    grad_thet_phi_G_X_grad_phi_G[:,:,3] =
                            (grad_thet_phi_G[:,:,1].*grad_phi_G[:,:,2] -
                             grad_thet_phi_G[:,:,2].*grad_phi_G[:,:,1])

    grad_thet_G_X_grad_phi_phi_G = zeros(size(unitNormalG))
    grad_thet_G_X_grad_phi_phi_G[:,:,1] =
                            (grad_thet_G[:,:,2].*grad_phi_phi_G[:,:,3] -
                             grad_thet_G[:,:,3].*grad_phi_phi_G[:,:,2])
    grad_thet_G_X_grad_phi_phi_G[:,:,2] =
                            (grad_thet_G[:,:,3].*grad_phi_phi_G[:,:,1] -
                             grad_thet_G[:,:,1].*grad_phi_phi_G[:,:,3])
    grad_thet_G_X_grad_phi_phi_G[:,:,3] =
                            (grad_thet_G[:,:,1].*grad_phi_phi_G[:,:,2] -
                             grad_thet_G[:,:,2].*grad_phi_phi_G[:,:,1])

    DiffUnitNormalG_thet = (grad_thet_thet_G_X_grad_phi_G +
                            grad_thet_G_X_grad_thet_phi_G)./JGbrev

    grad_thet_UnitNormalG = zeros(size(unitNormalG))
    grad_thet_UnitNormalG[:,:,1] = PP[:,:,1].*DiffUnitNormalG_thet[:,:,1] +
                                    PP[:,:,6].*DiffUnitNormalG_thet[:,:,2] +
                                    PP[:,:,5].*DiffUnitNormalG_thet[:,:,3]
    grad_thet_UnitNormalG[:,:,2] = PP[:,:,6].*DiffUnitNormalG_thet[:,:,1] +
                                    PP[:,:,2].*DiffUnitNormalG_thet[:,:,2] +
                                    PP[:,:,4].*DiffUnitNormalG_thet[:,:,3]
    grad_thet_UnitNormalG[:,:,3] = PP[:,:,5].*DiffUnitNormalG_thet[:,:,1] +
                                    PP[:,:,4].*DiffUnitNormalG_thet[:,:,2] +
                                    PP[:,:,3].*DiffUnitNormalG_thet[:,:,3]

    DiffUnitNormalG_phi = (grad_thet_phi_G_X_grad_phi_G +
                            grad_thet_G_X_grad_phi_phi_G)./JGbrev

    grad_phi_UnitNormalG = zeros(size(unitNormalG))
    grad_phi_UnitNormalG[:,:,1] = PP[:,:,1].*DiffUnitNormalG_phi[:,:,1] +
                                PP[:,:,6].*DiffUnitNormalG_phi[:,:,2] +
                                PP[:,:,5].*DiffUnitNormalG_phi[:,:,3]
    grad_phi_UnitNormalG[:,:,2] = PP[:,:,6].*DiffUnitNormalG_phi[:,:,1] +
                                PP[:,:,2].*DiffUnitNormalG_phi[:,:,2] +
                                PP[:,:,4].*DiffUnitNormalG_phi[:,:,3]
    grad_phi_UnitNormalG[:,:,3] = PP[:,:,5].*DiffUnitNormalG_phi[:,:,1] +
                                PP[:,:,4].*DiffUnitNormalG_phi[:,:,2] +
                                PP[:,:,3].*DiffUnitNormalG_phi[:,:,3]
    return grad_thet_UnitNormalG, grad_phi_UnitNormalG
end