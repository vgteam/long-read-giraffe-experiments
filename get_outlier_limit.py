# This is used for making split plots for outliers that are bigger than everything else
# For plotting runtimes and unmapped bases
# Takes a file of numbers with the numbers in the second column, and a string "big" or "small" for if we want the lower limit for plotting the big ones or the upper limit for plotting the small ones.

import sys
import numpy as np

runtimes = []
with open(sys.argv[1]) as in_file:
    for line in in_file:
        runtimes.append(float(line.split()[1]))


iqr = np.percentile(runtimes, 75) - np.percentile(runtimes, 25)
cutoff = np.percentile(runtimes, 75) + (1.5 * iqr)

bigs = list(filter(lambda x : x > cutoff, runtimes))
smalls = list(filter(lambda x : x <= cutoff, runtimes))

if not len(bigs) == 0 and sys.argv[2] == "big":

    print(min(bigs) / 60.0 * 0.80)

if not len(smalls) == 0 and sys.argv[2] == "small":

    print(max(smalls) * 1.1)


