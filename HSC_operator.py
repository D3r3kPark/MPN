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

H1, PL1, PM1, IL1, IM1, L1, M1 = m.run(1000, 1)



run_invasion_two_params = 0

if(run_invasion_two_params):


    print "Starting population is: " + str(H1 + PL1 + PM1 + IL1 + IM1 + L1 + M1)

    N1 = 10


    # Evaluate the eigenvalue with a range of r_N parameters and K_N parameters
    rN_val_start = 2
    rN_val_end = 0.1
    rN_val_range = np.linspace(rN_val_start, rN_val_end, 100)

    KN_val_start = 5
    KN_val_end = 7
    KN_val_range = np.logspace(KN_val_start, KN_val_end, 100)

    i = 0

    leading_eig_vals = []

    while i < len(rN_val_range):

        current_rN = rN_val_range[i]
        curr_eig_vals = []

        k = 0

        while k < len(KN_val_range):

            current_KN = KN_val_range[k]

            jac = hinv.jacobian(H1, PL1, PM1, IL1, IM1, L1, M1, N1, current_KN, current_rN)

            vals, vects, leadingVal = hinv.eigvals(jac)

            curr_eig_vals.append(leadingVal)

            k += 1

        leading_eig_vals.append(curr_eig_vals)

        i += 1


    matplotlib.rcParams.update({'font.size': 22})
    plt.figure(1)
    #ax = plt.subplot(111)
    plt.contourf(KN_val_range, rN_val_range, leading_eig_vals, cmap=plt.cm.coolwarm)
    plt.xscale('log')
    plt.ylabel("r_N")
    plt.xlabel("K_N")
    plt.title("Leading Eigenvalue vs K_N and r_N")
    plt.colorbar()
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #ax.legend(loc = 'center left', bbox_to_anchor = (1,0.5))
    plt.show(block=True)


run_invasion_one_param = 0

if(run_invasion_one_param):

    print "Starting population is: " + str(H1 + PL1 + PM1 + IL1 + IM1 + L1 + M1)

    N1 = 10

    #Evaluate the model with a range of just one parameter

    param_start =4.5
    param_end = 9
    param_range = np.logspace(param_start, param_end, 1000)

    i = 0
    leading_eig_vals = []
    while i < len(param_range):

        curr_param = param_range[i]
        jac = hinv.jacobian(H1, PL1, PM1, IL1, IM1, L1, M1, N1, curr_param, 1.0)
        vals, vects, leadingVal = hinv.eigvals(jac)

        leading_eig_vals.append(leadingVal)

        i += 1

    matplotlib.rcParams.update({'font.size': 22})
    plt.figure(1)
    plt.plot(param_range, leading_eig_vals, linewidth = 3)
    plt.ylabel("Leading Eigenvalue")
    plt.xlabel("K_N")
    plt.xscale('log')
    plt.title('Leading Eigenvalue vs K_N for r_N = 1.0')
    plt.show(block=True)
