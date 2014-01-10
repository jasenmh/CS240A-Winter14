#!/usr/bin/python

import commands
import re

MAXTIME = 600.0

f = open("hw1data.csv", "w")

for algo in range(1, 8):
        print('Algorithm ' + str(algo) +':')
        matdim = 2
        multtime = 0
        while multtime < MAXTIME:
            print('\t' + str(matdim) + 'x' + str(matdim) + 'matrices')
            cmd = "./matrix_multiply -n" + str(matdim) + " -a" + str(algo)
            o = commands.getoutput(cmd)
            multtime = float(re.search('Time \= (.*) sec, .*', o).group(1))
            f.write(str(multtime) + " ")
            print('\t\trun time ' + str(multtime))
            matdim = matdim << 1

        f.write('\n')

print('\nAnalysis complete.')
