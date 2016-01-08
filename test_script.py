# The objective of this program is to test things quickly
import sys
import os
import csv
import datetime
import HSC_differentiation_model as hmod
import threading
print 'start'
start = datetime.datetime.now()
import multiprocessing
from multiprocessing import Pool

num_cores = multiprocessing.cpu_count()

print 'Number of detected CPU cores: ' + str(num_cores)



def func(bound, threadName):
    h = hmod.HSC_model()

    for i in range(bound):

        h.run(1000,0)

        if i % 10 == 0 and i != 0:


            curr_time = datetime.datetime.now()
            elapsed = curr_time-start

            time_per_run = elapsed/i

            num_left = bound - i
            time_left = time_per_run * num_left

            proj_finish = curr_time + time_left
            total_time = bound*time_per_run

            print (threadName + ' - ' +
                  'Simulation ' + str(i) + ' of ' + str(bound) +
                  '. TPR: ' + str(time_per_run) +
                  '. Time left: ' + str(time_left) +
                  '. Est. finish: ' + str(proj_finish) +
                  '. Est. run time is: ' + str(total_time))


if __name__ == '__main__':
    start = datetime.datetime.now()
    jobs = []
    total_sims_needed = 1000
    sims_per_process = total_sims_needed/num_cores
    for i in range(num_cores):

        p = multiprocessing.Process(target=func, args=(sims_per_process,'Process ' + str(i)))
        jobs.append(p)
        p.start()


    # Join the processes
    for i in jobs:
        i.join()

    end = datetime.datetime.now()
    elapsed = end - start


    print str(total_sims_needed) + ' simulations took ' + str(elapsed) + ' on ' + str(num_cores) + ' cores.'
