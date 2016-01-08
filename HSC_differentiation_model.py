__author__ = 'Derek Park'
#The purpose of this program is to model the HSC transitions into the Common myeloid progenitors and common lymphoid progenitors.

import sys
import os
import csv
import pandas




print sys.executable
import numpy
import scipy
from scipy.integrate import odeint
import matplotlib

from matplotlib import pyplot as plt

plt.ion()


class HSC_model:

    #This class has all of the global parameters for the model as well as the initialized starting cell values

    # Define the global parameters:

    #Carrying capacities
    K_H = 1.0*10*2  # HSC carrying capacity
    K_L = 0.1*10**6 #Lymphocyte carrying capacity
    K_M = 0.9*10**6 #Myelocyte carrying capacity
    K_N = 10**6 # Tumor carrying capacity


    #Cell cycle rates
    r_T = 0.1
    r_H = 0.1  # HSC birth rate
    r_N = 0.1   # Cancer cell growth rate

    # Expansion and maturation proportions:

    B = 1
    a = 0.01
    alpha = 10**2

    # Symmetric differentiation feedback strengths
    s_L = 0.1
    s_M = s_L

    # Max feedback strength
    phi_max = 2

    #Death rates
    d_H = 0.1
    d_PL = 0.01
    d_PM = 0.01
    d_L = 0.25
    d_M = 0.25
    d_N = 0.1
    d_IL = 0.1
    d_IM = 0.1

    def __init__(self, initial_vals = None, initial_params = None):

        start_H = 10
        start_P_L = 10
        start_P_M = 10
        start_IL = 10
        start_IM = 10
        start_L = 10
        start_M = 10

        if initial_vals != None:

            start_H = initial_vals[0]
            start_P_L = initial_vals[1]
            start_P_M = initial_vals[2]
            start_IL = initial_vals[3]
            start_IM = initial_vals[4]
            start_L = initial_vals[5]
            start_M = initial_vals[6]

        self.H = start_H
        self.IL = start_IL
        self.IM = start_IM
        self.P_L = start_P_L
        self.P_M = start_P_M
        self.L = start_L
        self.M = start_M

        #If a custom parameter list was presented to the model, then use that instead of the default
        if initial_params != None:

            self.set_all_params(initial_params)


    def set_all_params(self, initial_params):

        #Carrying capacities
        self.K_H = initial_params[0]  # HSC carrying capacity
        self.K_L = initial_params[1]  # Lymphocyte carrying capacity
        self.K_M = initial_params[2]  # Myelocyte carrying capacity
        self.K_N = initial_params[3]  # Tumor carrying capacity


        #Cell cycle rates
        self.r_T = initial_params[4]  # HSC turnover rate
        self.r_H = initial_params[5]  # HSC birth rate
        self.r_N = initial_params[6]  # Cancer cell growth rate

        # Expansion and maturation proportions:

        self.B = initial_params[7]
        self.a = initial_params[8]
        self.alpha = initial_params[9]

        # Symmetric differentiation feedback strengths
        self.s_L = initial_params[10]
        self.s_M = self.s_L

        # Max feedback strength
        self.phi_max = initial_params[11]

        #Death rates
        self.d_H = initial_params[12]
        self.d_PL = initial_params[13]
        self.d_PM = initial_params[14]
        self.d_L = initial_params[15]
        self.d_M = initial_params[16]
        self.d_N = initial_params[17]
        self.d_IL = initial_params[18]
        self.d_IM = initial_params[19]


    # Return the parameters of this model
    def get_params(self):

        return [self.K_H, self.K_L, self.K_M, self.K_N,
                self.r_T, self.r_H, self.r_N,
                self.B, self.a, self.alpha,
                self.s_L, self.phi_max,
                self.d_H, self.d_PL, self.d_PM, self.d_L, self.d_M, self.d_N, self.d_IL, self.d_IM]

    # Set a specific parameter, identified by index of assignment in set_all_params()
    def set_specific_param(self, value, index):

        curr_params = self.get_params()

        curr_params[index] = value

        self.set_all_params(curr_params)


    # Simulate the model given initial parameters
    def growth(self, start_vals, time_series):

        # The current variable values
        H = start_vals[0]
        PL = start_vals[1]
        PM = start_vals[2]
        IL = start_vals[3]
        IM = start_vals[4]
        L = start_vals[5]
        M = start_vals[6]

        # Calculate the feedbacks
        phi_L = max( self.phi_max * (1 - L / self.K_L), 0)
        phi_M = max( self.phi_max * (1 - M / self.K_M), 0)

        dHdt = (self.r_T + self.r_H *  (1 - H/self.K_H)) * H * (phi_L * phi_M  - self.s_L * phi_L - self.s_M * phi_M) - self.d_H * H
        dPLdt = (self.r_T + self.r_H * (1 - H/self.K_H)) * H * phi_L * ( 1 + 2 * self.s_L) - self.a * PL - self.d_PL * PL
        dPMdt = (self.r_T + self.r_H * (1 - H/self.K_H)) * H * phi_M * ( 1 + 2 * self.s_M) - self.a * PM - self.d_PM * PM
        dILdt = PL * self.alpha - phi_L * self.B * IL - self.d_IL * IL
        dIMdt = PM * self.alpha - phi_M * self.B * IM - self.d_IM * IM
        dLdt = phi_L * self.B * IL - self.d_L * L
        dMdt = phi_M * self.B * IM - self.d_M * M

        return [dHdt, dPLdt, dPMdt, dILdt, dIMdt, dLdt, dMdt]

    # Check to see whether a set of population values satisfies what we consider to be an acceptable homeostasis
    def isValid(self, parameters):

        # First of all, test to see if the population ever drops below 1 cell.
        # If so, then failure.

        if (min(parameters[0]) < 1 or
            min(parameters[1]) < 1 or
            min(parameters[2]) < 1 or
            min(parameters[3]) < 1 or
            min(parameters[4]) < 1 or
            min(parameters[5]) < 1 or
            min(parameters[6]) < 1 ):

            return -1

        # So if all cell populations remained alive, then we test to see if the homeostasis we get is valid
        H1 = parameters[0][-1]
        PL1 = parameters[1][-1]
        PM1 = parameters[2][-1]
        IL1 = parameters[3][-1]
        IM1 = parameters[4][-1]
        L1 = parameters[5][-1]
        M1 = parameters[6][-1]

        hem_sum = H1 + PL1 + PM1 + IL1 + IM1 + L1 + M1

        myeloid_frac = (PM1 + IM1 + M1)/hem_sum
        lymphoid_frac = (PL1 + IL1 + L1)/hem_sum

        correct_homeostasis = 0
        # Is the HSC fraction less than 0.1%?
        if H1/hem_sum < 0.001:
            correct_homeostasis += 1

        # Are progenitors more common than HSCs, but less common than immature cells?
        if PL1 > H1 and PM1 > H1:
            if PL1 < IL1 and PM1 < IM1:

                correct_homeostasis += 1

        # Are myeloid and lymphoid fractions in the appropriate ranges?
        if myeloid_frac > 0.6 and myeloid_frac < 0.9:
            if lymphoid_frac> 0.1 and lymphoid_frac< 0.4:

                correct_homeostasis += 1

        if correct_homeostasis == 3:
            return 1
        else:
            return -1


    #Run the model and return the results as well as whether it is a valid model.
    def run(self, time_end, plot):

        init_conds = [self.H, self.IL, self.IM, self.P_L, self.P_M, self.L, self.M]

        #Create the time axis
        time_seq = numpy.linspace(0,time_end, int(time_end*500), endpoint=False)


        #Run the model

        results = odeint(self.growth, init_conds, time_seq)

        #Get the resulting series
        H = results[:,0]
        PL = results[:,1]
        PM = results[:,2]
        IL = results[:,3]
        IM = results[:,4]
        L = results[:,5]
        M = results[:,6]



        if plot == 1:

            matplotlib.rcParams.update({'font.size': 22})

            plt.figure(1)
            ax = plt.subplot(111)

            ax.plot(time_seq, H, label = "HSCs", linewidth = 5.0)
            ax.plot(time_seq, IL, label = "IL", linewidth = 5.0)
            ax.plot(time_seq, IM, label = "IM", linewidth = 5.0)
            ax.plot(time_seq, PL, label = "CLPs", linewidth = 5.0)
            ax.plot(time_seq, PM, label = "CMPs", linewidth = 5.0)
            ax.plot(time_seq, L, label = "Lymphocytes", linewidth = 5.0)
            ax.plot(time_seq, M, label = "Myelocytes", linewidth = 5.0)

            plt.yscale("log")
            plt.ylabel("Cell number")
            plt.xlabel("Time")
            plt.title("Hematopoietic Growth")
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc = 'center left', bbox_to_anchor = (1,0.5))



            # Now print the final values of the system
            print "Final values:"

            sum = H[-1] + PL[-1] + PM[-1] + L[-1] + M[-1] + IL[-1] + IM[-1]
            print "Total number of cells: " + str(sum)
            print "H: " + str(H[-1])
            print "CLP: " + str(PL[-1])
            print "CMP: " + str(PM[-1])
            print "IL: " + str(IL[-1])
            print "IM: " + str(IM[-1])
            print "L: " + str(L[-1])
            print "M: " + str(M[-1])

            print "Final Percentages"
            print "H: " + str(H[-1]/sum*100) + " %"
            print "CLP: " + str(PL[-1]/sum*100) + " %"
            print "CMP: " + str(PM[-1]/sum*100) + " %"
            print "IL: " + str(IL[-1]/sum*100) + " %"
            print "IM: " + str(IM[-1]/sum*100) + " %"
            print "L: " + str(L[-1]/sum*100) + " %"
            print "M: " + str(M[-1]/sum*100) + " %"
            print "Total Myeloid: " + str((PM[-1] + IM[-1] + M[-1])/sum*100)
            print "Total Lymphoid: " + str((PL[-1] + IL[-1] + L[-1])/sum*100)

            print "Final Feedbacks:"
            print "f_L: " + str(max( self.phi_max - L[-1] / self.K_L, 0))
            print "f_M: " + str(max( self.phi_max - M[-1] / self.K_M, 0))

            # plt.figure(3)
            # # The slices will be ordered and plotted counter-clockwise.
            # labels = 'HSCs', 'CLPs', 'L', 'CMPs', 'M'
            # sizes = [H[-1], PL[-1], L[-1], PM[-1],M[-1]]
            # colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'red']
            # explode = (0, 0, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')
            #
            # plt.pie(sizes, labels=labels, colors=colors,autopct='%1.1f%%', shadow=True, startangle=90)
            # # Set aspect ratio to be equal so that pie is drawn as a circle.
            # plt.axis('equal')
            plt.show(block = True)
            plt.close()

        # Lastly, evaluate whether or the final population numbers represent a biologically realistic homeostasis:
        validity = self.isValid([H, PL, PM, IL, IM, L, M])

        return time_seq, H, PL, PM, IL, IM, L, M, validity

