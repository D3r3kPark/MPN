__author__ = 'Derek'

#The purpose of this program is to perform an evolutionary invasion analysis on the model I have built for HSCs.


import numpy as np
import scipy
import matplotlib
test = 0

if test:
    #This is to test the process of finding the eigenvalues of a matrix:


    a = np.matrix([[1,2],[2,4]])


    w,v =  np.linalg.eig(a)
    print w
    print v
    print np.linalg.det(a)



def jacobian(H, PL, PM, IL, IM, L, M, N, current_KN, current_rN):

    rH = 0.1
    rT = 0.1
    rN = current_rN

    KH = 100
    KM = 0.7*10**6
    KL = 0.3*10**6
    KN = current_KN

    B = 0.25
    a = 0.001
    alpha = 10**2

    SL = 1.2
    SM = 1.2


    dH = 0.1
    dPL = 0.01
    dPM = 0.01
    dL = 0.25
    dM = 0.25
    dN = 0.1
    dIL = 0.1
    dIM = 0.1

    fMax = 2.0
    fL = min(fMax*(1 - L / KL), fMax)
    fM = min(fMax*(1 - M / KM), fMax)

    dfL = -1.0 * fMax / KL
    dfM = -1.0 * fMax / KM


    # The elements of the Jacobian

    J_HH = (rT + rH * (1 - (2 * H)/KH)) * (fL*fM - fL*SL - fM*SM) - dH
    J_HPL = 0
    J_HPM = 0
    J_HIL = 0
    J_HIM = 0
    J_HL = (rT + rH  * (1 - H/KH)) * H * (dfL * (fM - SL))
    J_HM = (rT + rH  * (1 - H/KH)) * H * (dfM * (fL - SM))
    J_HN = 0

    J_PLH = (rT + rH * (1 - (2*H)/KH)) * fL * (1 + 2*SL)
    J_PLPL = -1 * a - dPL
    J_PLPM = 0
    J_PLIL = 0
    J_PLIM = 0
    J_PLL = (rT + rH * ( 1 - H/KH)) * H * (1 + 2*SL) * dfL
    J_PLM  = 0
    J_PLN = 0

    J_PMH = (rT + rH * ( 1 - (2*H)/KH))  * fM * (1 + 2*SM)
    J_PMPL = 0
    J_PMPM = -1 * a - dPM
    J_PMIL = 0
    J_PMIM = 0
    J_PML = 0
    J_PMM = (rT + rH * (1 - H/KH)) * H * (1 + 2*SM) * dfM
    J_PMN = 0

    J_ILH = 0
    J_ILPL = alpha
    J_ILPM = 0
    J_ILIL = -1 * fL * B - dIL
    J_ILIM = 0
    J_ILL = -1 * dfL * B * IL
    J_ILM = 0
    J_ILN = 0

    J_IMH = 0
    J_IMPL = 0
    J_IMPM = alpha
    J_IMIL = 0
    J_IMIM = -1 * fM * B - dIM
    J_IML = 0
    J_IMM = -1 * dfM * B * IM
    J_IMN = 0

    J_LH = 0
    J_LPL = 0
    J_LPM = 0
    J_LIL = fL * B
    J_LIM = 0
    J_LL = dfL * B * IL - dL
    J_LM = 0
    J_LN = 0

    J_MH = 0
    J_MPL = 0
    J_MPM = 0
    J_MIL = 0
    J_MIM = fM * B
    J_ML = 0
    J_MM = dfM * B * IM - dM
    J_MN = 0

    J_NH = rN * N * (-1/KN)
    J_NPL = rN * N * (-1/KN)
    J_NPM = rN * N * (-1/KN)
    J_NIL = rN * N * (-1/KN)
    J_NIM = rN * N * (-1/KN)
    J_NL = rN * N * (-1/KN)
    J_NM = rN * N * (-1/KN)
    J_NN = rN * (1 - (H + PL + PM + IL + IM + L + M)/KN) - dN

    J = np.matrix([[J_HH,   J_HPL,  J_HPM,  J_HIL,  J_HIM,  J_HL,  J_HM,  J_HN],
                   [J_PLH, J_PLPL, J_PLPM, J_PLIL, J_PLIM, J_PLL, J_PLM, J_PLN],
                   [J_PMH, J_PMPL, J_PMPM, J_PMIL, J_PMIM, J_PML, J_PMM, J_PMN],
                   [J_ILH, J_ILPL, J_ILPM, J_ILIL, J_ILIM, J_ILL, J_ILM, J_ILN],
                   [J_IMH, J_IMPL, J_IMPM, J_IMIL, J_IMIM, J_IML, J_IMM, J_IML],
                   [J_LH, J_LPL, J_LPM, J_LIL, J_LIM, J_LL, J_LM, J_LN],
                   [J_MH, J_MPL, J_MPM, J_MIL, J_MIM, J_ML, J_MM, J_MN],
                   [J_NH, J_NPL, J_NPM, J_NIL, J_NIM, J_NL, J_NM, J_NN]])

    return J


# Take a matrix and return its eigenvalues
def eigvals(J):

    vals, vects = np.linalg.eig(J)

    leadingVal = vals.flat[abs(vals).argmax()]

    return vals, vects, leadingVal
