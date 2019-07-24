#!/usr/bin/env python

### LIBRARIES ###

import sys       # https://docs.python.org/2/library/sys.html
import argparse  # https://docs.python.org/2/library/argparse.html
import csv       # https://docs.python.org/2/library/csv.html
import math      # https://docs.python.org/2/library/math.html


### PARSE ARGUMENTS ###

# validator for NACA 4-series definition
def valid_naca4(string):
    if (len(string) != 4) or (not string.isdigit()):
        msg = "%r is not a 4 digit NACA aerofoil definition" % string
        raise argparse.ArgumentTypeError(msg)
    return string

# command line parser
parser = argparse.ArgumentParser(
    prog='naca4series',
    description='Create a NACA 4-series aerofoil',
    epilog='See: https://en.wikipedia.org/wiki/NACA_airfoil#Four-digit_series'
)

# parser group for aerofoil definition
group_def = parser.add_argument_group('Aerofoil definition')
group_def.add_argument(
    'MPTT', metavar='MPTT',
    type=valid_naca4,
    help='4 digits specifying MPTT (e.g. 0012)'
)
group_def.add_argument(
    '-r', '--resolution', dest='res', metavar='N',
    type=int, default='100', choices=range(3, 1000),
    help='number of points on camber line (default: %(default)s)'
)
group_def.add_argument(
    '-s', '--spacing', dest='spacing', metavar='TYPE',
    type=str, default='cos', choices=('cos', '2cos', 'lin',),
    help='spacing (default: %(default)s)'
)

# parser group for input/output
group_IO = parser.add_argument_group('Output')
group_IO.add_argument(
    'outfile', metavar='FILE.csv',
    nargs='?', type=argparse.FileType('w'), default=sys.stdout,
    help='file to write to (default: stdout)'
)

# process
args = parser.parse_args()  # reads from sys.argv


### CALCULATE AEROFOIL ###

# auxiliary functions (to avoid using NumPy)
def linspace(a, b, N=100):
    if N == 1:
        return b
    h = (b - a) / (N - 1)
    return [a + h*n  for n in range(N)]
def cos(x):
    # using https://docs.python.org/2/library/functions.html#map
    return map(math.cos, x)
def dealzeros(N, M):
    return [[0]*N for _ in range(M)]

# camberline x-coordinates
if args.spacing == 'cos':
    x = [(1.0-xi)/2.0  for xi in cos(linspace(0.0, math.pi, args.res))]
elif args.spacing == '2cos':
    x = [1.0-xi  for xi in cos(linspace(0.0, math.pi/2.0, args.res))]
elif args.spacing == 'lin':
    x = linspace(0.0, 1.0, args.res)

# extract values from NACA definition
M = float(args.MPTT[0])/100    # maximum camber
P = float(args.MPTT[1])/10     # percentage location of maximum thickness
T = float(args.MPTT[2:3])/100  # maximum thickness

# calculate 
Npts = len(x)
x_C, y_C, x_U, y_U, x_L, y_L = dealzeros(Npts, 6)  # initialise all
for i in range(Npts):
    x_c = x[i]

    # thickness
    y_t = T/0.20 * ( 0.2969 * math.sqrt(x_c)
                   - 0.1260 * x_c
                   - 0.3516 * x_c**2
                   + 0.2843 * x_c**3
                   - 0.1015 * x_c**4 )
    if M==0 and P==0:
        # symmetrical
        y_c = 0
        dycdx = 0
    else:
        # cambered
        if (0 <= x_c) and (x_c <= P):
            y_c = M/(P**2)*(2*P*x_c-x_c**2)
            dycdx = 2*M/(P**2)*(P-x_c)
        else:
            y_c = M/((1-P)**2)*((1-2*P)+2*P*x_c-x_c**2)
            dycdx = 2*M/((1-P)**2)*(P-x_c)
    # coordinates
    theta = math.atan(dycdx)
    # mean camber line
    x_C[i] = x_c
    y_C[i] = y_c
    # upper surface
    x_U[i] = x_c - y_t*math.sin(theta)
    y_U[i] = y_c + y_t*math.cos(theta)
    # lower surface
    x_L[i] = x_c + y_t*math.sin(theta)
    y_L[i] = y_c - y_t*math.cos(theta)


### EXPORT DATA ###

# prepare file
writer = csv.writer(
    args.outfile,
    delimiter=',',
    quotechar='|',
    quoting=csv.QUOTE_MINIMAL
)
# don't write a header row ["x","y","z"], to be universally readable!

for i in range(2*Npts):  # 2x for both lower and upper

    if i < Npts:
        # go backwards
        x_i = x_U[-1-i]
        y_i = y_U[-1-i]
    elif i == Npts:
        continue  # skip this one (duplicate LE)
        #NOTE: TE is not closed by definition
    else:
        # go forward
        x_i = x_L[i-Npts]
        y_i = y_L[i-Npts]
    z_i = 0

    # write data point
    writer.writerow([x_i, y_i, z_i])
