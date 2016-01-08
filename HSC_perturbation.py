__author__ = 'Derek'
import HSC_differentiation_model as hmod
import HSC_invasion as hinv
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import os

H_0 = 10
P_L_0 = 10
P_M_0 = 10
IL_0 = 10
IM_0 = 10
L_0 = 10
M_0 = 10

initial_pops = [H_0, P_L_0, P_M_0, IL_0, IM_0, L_0, M_0]
m = hmod.HSC_model(initial_pops)

valid_homeostasis_test = 0

if valid_homeostasis_test:
    var_start = 100
    var_end = 1000

    var_range = np.linspace(var_start, var_end, 100)


    valid_param = []
    i = 0

    while i < len(var_range):
        curr_var = var_range[i]

        H1, PL1, PM1, IL1, IM1, L1, M1 = m.run(1000, 0, curr_var)

        hem_sum = H1 + PL1 + PM1 + IL1 + IM1 + L1 + M1

        myeloid_frac = (PM1 + IM1 + M1)/hem_sum
        lymphoid_frac = (PL1 + IL1 + L1)/hem_sum

        correct_homeostasis = 0

        if H1/hem_sum < 0.001:
            correct_homeostasis += 1

        if PL1 > H1 and PM1 > H1:
            if PL1 < IL1 and PM1 < IM1:

                correct_homeostasis += 1

        if myeloid_frac > 0.6 and myeloid_frac < 0.9:
            if lymphoid_frac> 0.1 and lymphoid_frac< 0.4:

                correct_homeostasis += 1

        if correct_homeostasis == 3:
            valid_param.append(1)
        else:
            valid_param.append(0)

        i += 1

    matplotlib.rcParams.update({'font.size': 22})
    plt.figure(1)
    ax = plt.subplot(111)

    ax.plot(var_range, valid_param, label = "Valid Homeostases", linewidth = 5.0)

    #plt.yscale("log")
    plt.ylabel("Validity")
    plt.xlabel("r_T")
    plt.title("Valid Homeostases")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc = 'center left', bbox_to_anchor = (1,0.5))
    plt.show(block = True)
    plt.close()

two_param_perturbation = 1

if two_param_perturbation:

    N1 = 10

    var_start = 2
    var_end = 3
    #Assumes that the homeostases generated over this range are valid
    var_range = np.logspace(var_start, var_end, 200)

    rN_start = 0.1
    rN_end = 10.0
    rN_range = np.linspace(rN_start, rN_end, 100)

    i = 0
    leading_eigvals = []
    while i < len(var_range):

        curr_var = var_range[i]

        H1, PL1, PM1, IL1, IM1, L1, M1 = m.run(1000, 0, curr_var)

        curr_eigvals = []

        k = 0

        while k < len(rN_range):

            curr_rN = rN_range[k]

            jac = hinv.jacobian(H1, PL1, PM1, IL1, IM1, L1, M1, N1, curr_var, curr_rN)

            vals, vects, leadingVal = hinv.eigvals(jac)

            curr_eigvals.append(leadingVal)
            print 'Running simulation ' + str(i) + ":" + str(k)
            k += 1

        leading_eigvals.append(curr_eigvals)
        i += 1

    matplotlib.rcParams.update({'font.size': 22})
    plt.figure(1)
    #ax = plt.subplot(111)
    plt.contourf(rN_range, var_range, leading_eigvals, cmap=plt.cm.coolwarm)
    plt.yscale('log')
    plt.ylabel("Alpha")
    plt.xlabel("r_N")
    plt.title("Leading Eigenvalue vs Alpha and r_N")
    plt.colorbar()
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #ax.legend(loc = 'center left', bbox_to_anchor = (1,0.5))
    plt.show(block=True)
