#!/usr/bin/python

import commands
import re

MAXTIME = 600.0
ITERATIONS = 4
DEBUG = False

f = open("hw1data.csv", "w")

for algo in range(1, 8):
        if DEBUG:
          print('Algorithm ' + str(algo) +':')

        matdim = 2
        multtime = 0
        tottime = 0

        while multtime < MAXTIME:

          itercount = 0
          while itercount < ITERATIONS:

            if DEBUG:
              print('\t' + str(matdim) + 'x' + str(matdim) + 'matrices, run ' + str(itercount + 1))

            cmd = "./matrix_multiply -n" + str(matdim) + " -a" + str(algo)
            o = commands.getoutput(cmd)
            multtime = float(re.search('Time \= (.*) sec, .*', o).group(1))
            tottime = tottime + multtime

            if DEBUG:
              print('\t\trun time ' + str(multtime))

            itercount = itercount + 1

          f.write(str(tottime/float(ITERATIONS)) + ", ")
          matdim = matdim << 1

        f.write('\n')

if DEBUG:
  print('\nAnalysis complete.')
