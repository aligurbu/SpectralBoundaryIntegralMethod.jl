

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