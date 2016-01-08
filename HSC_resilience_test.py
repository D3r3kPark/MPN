# The objective of this program is to test the resilience of the equilibria generated
import sys
import os
import csv
import datetime
import math
import pandas
operation_start_time = datetime.datetime.now()



print sys.executable
import numpy as np
import scipy
import matplotlib
import HSC_differentiation_model as hmod
from matplotlib import pyplot as plt


class resilience_test:

    default_params = []


    def __init__(self):
        self.default_params = hmod.HSC_model().get_params()

    # Simulate recovery of the HSC populations from start_pops to end_pops.
    def simulate_recovery(self, start_pops, end_pops, param_set):

        #Identify the target values.
        target_H = end_pops[0]
        target_PL = end_pops[1]
        target_PM = end_pops[2]
        target_IL = end_pops[3]
        target_IM = end_pops[4]
        target_L = end_pops[5]
        target_M = end_pops[6]

        #Create the HSC model with the given parameters
        model = hmod.HSC_model(start_pops, param_set)

        #Simulate it for 1000 timesteps and get the results

        time, H, PL, PM, IL, IM, L, M, isRealistic = model.run(500, 0)

        #Test to see at what time it reaches the end populations (>=)

        index_of_recovery = -1
        i = 0

        while i < len(time):

            curr_H = H[i]
            curr_PL = PL[i]
            curr_PM = PM[i]
            curr_IL = IL[i]
            curr_IM = IM[i]
            curr_L = L[i]
            curr_M = M[i]



            # Check to see if the populations have recovered
            if (curr_H >= target_H and
                curr_PL >= target_PL and
                curr_PM >= target_PM and
                curr_IL >= target_IL and
                curr_IM >= target_IM and
                curr_L >= target_L and
                curr_M >= target_M):

                index_of_recovery = i

                i = len(time)
            i += 1

        recovery_time = time[index_of_recovery]

        return recovery_time, [time, H, PL, PM, IL, IM, L, M]

    # Scan two parameter spaces and return the resiliency
    def scan(self, range1, var1_index, var1_label, range2, var2_index, var2_label):

        if len(range1) != len(range2):
            print 'Warning: This is not a square parameter space'

        # The header for the full output file
        full_output = [[var1_label, var2_label, 'Recovery from 10 percent', 'Biologically reasonable homeostasis? 1 = Yes']]

        recovery_times = [[-1 for x in range1] for x in range2]
        valid_runs = [[0 for x in range1] for x in range2]


        # Get the time the search was started:

        start_time = datetime.datetime.now()
        # Scan through the parameter space
        i1 = 0

        while i1 < len(range1):

            i2 = 0

            while i2 < len(range2):
                curr_sim_number = i1*len(range2) + i2 + 1

                if curr_sim_number % 10 == 0:
                    # Get the current time:
                    curr_time = datetime.datetime.now()

                    time_elapsed = curr_time - start_time
                    total_sim_number = len(range1)*len(range2)

                    time_per_sim = time_elapsed/curr_sim_number

                    proj_time_remain = time_per_sim * (total_sim_number-curr_sim_number)
                    proj_finish = curr_time + proj_time_remain
                    proj_total_time = total_sim_number*time_per_sim
                    print ('Simulating ' + str(curr_sim_number) + ' of ' + str(total_sim_number) +
                           '. Time elapsed: ' + str(time_elapsed) +
                           '. Time remaining: ' + str(proj_time_remain) +
                           '. Projected Finish: ' + str(proj_finish) +
                           '. Total est. time: ' + str(proj_total_time))

                curr_param_set = self.default_params

                # Choose the parameters for the new set.
                curr_param_set[var1_index] = range1[i1]
                curr_param_set[var2_index] = range2[i2]

                # Create the model based on this parameter set and run it for 1000 timesteps
                starting_cell_pops = [10,10,10,10,10,10,10]

                curr_model = hmod.HSC_model(starting_cell_pops,curr_param_set)

                # Run the model for 1000 timesteps to generate the homeostasis
                t, H, PL, PM, IL, IM, L, M, valid = curr_model.run(500, 0)


                # Store the biological validity of this parameter set
                valid_runs[i1][i2] = valid


                # Only do any analysis if the model is biologically valid:
                curr_recov_time = -1
                if valid == 1:

                    # Determine the final homeostatic population counts
                    final_pop = [H[-1], PL[-1], PM[-1], IL[-1], IM[-1], L[-1], M[-1]]

                    # Determine the post-disturbance population counts
                    disturb_pop = [H[-1]*0.9, PL[-1]*0.9, PM[-1]*0.9, IL[-1]*0.9, IM[-1]*0.9, L[-1]*0.9, M[-1]*0.9]

                    # Determine the recovery time for this disturbance
                    curr_recov_time, curr_recov_data = self.simulate_recovery(disturb_pop,final_pop,curr_param_set)

                recovery_times[i1][i2] = curr_recov_time
                full_output.append([range1[i1], range2[i2], curr_recov_time, valid])


                print_output = False

                if print_output and valid == 1:

                    # First we need to concatenate the data sets:
                    t = t.tolist()
                    H = H.tolist()
                    PL = PL.tolist()
                    PM = PM.tolist()
                    IL = IL.tolist()
                    IM = IM.tolist()
                    L = L.tolist()
                    M = M.tolist()
                    


                    t2_start = max(t)
                    #print t2_start
                    t2 = curr_recov_data[0]

                    print min(t2)
                    print max(t2)

                    H2 = curr_recov_data[1]
                    PL2 = curr_recov_data[2]
                    PM2 = curr_recov_data[3]
                    IL2 = curr_recov_data[4]
                    IM2 = curr_recov_data[5]
                    L2 = curr_recov_data[6]
                    M2 = curr_recov_data[7]

                    print len(H2)

                    for i in range(len(t2)):

                        t.append(t2[i] + t2_start)
                        H.append(H2[i])
                        PL.append(PL2[i])
                        PM.append(PM2[i])
                        IL.append(IL2[i])
                        IM.append(IM2[i])
                        L.append(L2[i])
                        M.append(M2[i])

                    matplotlib.rcParams.update({'font.size': 22})

                    plt.figure(1)
                    ax = plt.subplot(111)

                    ax.plot(t, H, label = "HSCs", linewidth = 5.0)
                    ax.plot(t, IL, label = "IL", linewidth = 5.0)
                    ax.plot(t, IM, label = "IM", linewidth = 5.0)
                    ax.plot(t, PL, label = "CLPs", linewidth = 5.0)
                    ax.plot(t, PM, label = "CMPs", linewidth = 5.0)
                    ax.plot(t, L, label = "Lymphocytes", linewidth = 5.0)
                    ax.plot(t, M, label = "Myelocytes", linewidth = 5.0)

                    plt.yscale("log")
                    plt.ylabel("Cell number")
                    plt.xlabel("Time")
                    plt.title("Hematopoietic Growth for " +
                              var1_label + ' = ' + str(range1[i1]) + ', ' +
                              var2_label + ' = ' + str(range2[i2]))
                    box = ax.get_position()
                    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                    ax.legend(loc = 'center left', bbox_to_anchor = (1,0.5))
                    file_name = var1_label + '_' + str(math.floor(range1[i1]*1000)/1000) + '_' + var2_label + '_' + str(math.floor(range2[i2]*1000)/1000)
                    plt.savefig(file_name + '.png')
                    plt.show(block = True)
                    plt.close()




                i2 += 1
            i1 += 1


        # Now print the output to file as a csv file
        outputfile = open(var1_label + '_vs_' + var2_label + '.csv', 'w')
        outputWriter = csv.writer(outputfile)
        for x in full_output:
            outputWriter.writerow(x)
        outputfile.close()

        # Return the results as well as the original ranges
        return recovery_times, valid_runs

    # Run a full resilience analysis over a parameter space or for a specific simulation.
    def run(self, axis_size, specific=0):
        # Set up the two variables and the ranges to be scanned
        var1_label = 'r_T'
        var1_index =4

        var2_label = 'alpha'
        var2_index = 9

        # If a specific area is not desired, then scan an area
        if specific == 0:

            var1_range = np.linspace(0.01, 2, axis_size)

            var2_range = np.logspace(2, 3, axis_size)
            # Scan the parameter area
            recovery_times, valid_runs = self.scan(var1_range,var1_index, var1_label, var2_range, var2_index, var2_label)
        # Show a specific simulation
        else:

            spec_var1 = float(input(var1_label + " value? \n"))
            spec_var2 = float(input(var2_label + " value? \n"))

            recovery_times, valid_runs = self.scan([spec_var1],var1_index, var1_label, [spec_var2], var2_index, var2_label)


    # Graph the results as loaded from a CSV
    def graph(self, axis_size):

        # Reconstruct the array
        data_file = open('r_T_vs_alpha.csv')
        data_reader = csv.reader(data_file)

        x_axis = [[-1 for x in range(axis_size)] for x in range(axis_size)]
        y_axis = [[-1 for x in range(axis_size)] for x in range(axis_size)]
        resilience_vals = [[-1 for x in range(axis_size)] for x in range(axis_size)]
        valid_vals = [[-1 for x in range(axis_size)] for x in range(axis_size)]
        x_index = 0
        y_index = 0
        for row in data_reader:
            if data_reader.line_num != 1:
                x_axis[x_index][y_index] = row[0]
                y_axis[x_index][y_index] = row[1]
                resilience_vals[x_index][y_index] = row[2]
                valid_vals[x_index][y_index] = row[3]

                y_index += 1
                y_index = y_index % axis_size

                # Check to see if we need to increment x
                if y_index == 0:

                    x_index += 1

                x_index = x_index % axis_size

        # Plot the recovery times and the valid runs
        # contourf(X,Y,Z) python documentation states: X and Y must both be 2-D with the same shape as Z,
        # or they must both be 1-D such that len(X) is the number of columns in Z and len(Y) is the number of rows in Z.
        # Matrices are always (row,column)
        matplotlib.rcParams.update({'font.size': 22})
        plt.figure(1)
        plt.contourf(x_axis, y_axis, resilience_vals, cmap=plt.cm.coolwarm)
        #plt.yscale('log')
        plt.xlabel('r_T')
        plt.ylabel('alpha')
        plt.title('Recovery times')# for ' + str(var1_label) + ' vs. ' + str(var2_label))
        plt.colorbar()

        plt.figure(2)
        plt.contourf(x_axis, y_axis, valid_vals)
        #plt.yscale('log')
        plt.xlabel('r_T')
        plt.ylabel('alpha')
        plt.title('Valid homeostasis')# in ' + str(var1_label) + ' vs. ' + str(var2_label))
        plt.colorbar()

        plt.show(block=True)


res = resilience_test()
dimension = 25
#res.run(dimension)
res.graph(dimension)
look_specific = 'n'

while look_specific == 'y':
    res.run(specific=1)
    look_specific = str(input('Run again? (y/n) \n'))

operation_finish_time = datetime.datetime.now()
print 'Operation started at: ' + str(operation_start_time) + ' and finished at: ' + str(operation_finish_time) + '. Total time: ' + str(operation_finish_time-operation_start_time)
