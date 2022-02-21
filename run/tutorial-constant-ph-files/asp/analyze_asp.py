#!/usr/bin/env python
from __future__ import print_function
from sys import argv
from math import sqrt
from collections import OrderedDict

from numpy import loadtxt, int32, concatenate


def read_cphlog(cphlogfilename):
    # Read the header and determine the pH value.
    pH = None
    for line in open(cphlogfilename, 'r'):
        if not line.startswith('#'):
          break
        tokens = line.lstrip('#').strip().split()
        try:
            pH = float(tokens[tokens.index('pH')+1])
        except ValueError:
            pass
    if pH is None:
        raise ValueError('No pH info in header of %s'%cphlogfilename)

    # Read the lambda_1 and lambda_2 columns (ignore the time stamp).
    lambda_ = loadtxt(cphlogfilename, int32, usecols=(1, 2))

    return lambda_, pH 

# Read and concatenate any number of files, sorting by the pH (use the pH value
# as a dictionary key). Since standard dict objects cannot be sorted, used the
# special OrderedDict class to sort the pH values in ascending order.
#
tmp = OrderedDict()
for cphlogfilename in argv[1:]:
    lambda_, pH = read_cphlog(cphlogfilename) 
    if pH not in tmp:
        tmp[pH] = lambda_
    else:
        tmp[pH] = concatenate((tmp[pH], lambda_))

pHs = tmp.keys()
pHs.sort()
cphlogs = OrderedDict()
for pH in pHs:
    cphlogs[pH] = tmp[pH]
del tmp    

# Each value in cphlogs is now an ndarray with shape (nsamples, 2), where
# nsamples is the total number of simulation cycles at that pH. The second axis
# of the array (axis=1) corresponds to the two protonation sites.
#
# Compute chi as the joint probability of _either_ site 1 or site 2 being
# occupied:
#
print('#  pH   frac')
for pH, lambda_ in cphlogs.iteritems():
    chi = lambda_[:, 0]*(1 - lambda_[:, 1]) + (1 - lambda_[:, 0])*lambda_[:, 1]
    prot_frac = chi.mean()
    prot_frac_err = sqrt(prot_frac*(1 - prot_frac) / chi.size)
    print('%5.2f %6.4f %6.4f'%(pH, prot_frac, prot_frac_err))

