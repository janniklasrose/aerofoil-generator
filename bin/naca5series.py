#!/usr/bin/env python3

"""
Create a NACA 5-series aerofoil.

For usage, run naca5series with the -h/--help option.
Definition: https://en.wikipedia.org/wiki/NACA_airfoil#Five-digit_series
See also: NACA-TR-537 and NACA-TR-610.

Copyright (c) 2026 Jan Niklas Rose
"""

### LIBRARIES ###

import sys  # https://docs.python.org/2/library/sys.html
import argparse  # https://docs.python.org/2/library/argparse.html
import csv  # https://docs.python.org/2/library/csv.html
import math  # https://docs.python.org/2/library/math.html

### CONSTANTS ###

# mean line tables for design lift coefficient cld = 0.3
# table[(P, Q)] = (r, k1)
STANDARD_MEANLINES = {
    ("1", "0"): (0.0580, 361.4),  #  5%
    ("2", "0"): (0.1260, 51.64),  # 10%
    ("4", "0"): (0.2900, 6.643),  # 20%
    ("3", "0"): (0.2025, 15.957),  # 15%
    ("5", "0"): (0.3910, 3.230),  # 25%
}
REFLECT_MEANLINES = {
    # ('1', '1') = NACA 211 was not published in NACA-TR-537.
    ("2", "1"): (0.1300, 51.990),  # 10%
    ("3", "1"): (0.2170, 15.793),  # 15%
    ("4", "1"): (0.3180, 6.520),  # 20%
    ("5", "1"): (0.4410, 3.191),  # 25%
}
MEANLINE_TABLE = {**STANDARD_MEANLINES, **REFLECT_MEANLINES}


### PARSE ARGUMENTS ###


# validator for NACA 5-series definition
def valid_naca5(string):
    """Validate a NACA 5-series definition."""
    if (len(string) != 5) or (not string.isdigit()):
        msg = "%r is not a 5 digit NACA aerofoil definition" % string
        raise argparse.ArgumentTypeError(msg)
    if (string[1], string[2]) not in MEANLINE_TABLE:
        msg = (
            "%r is not a supported NACA 5-series mean line (leading 3 digits)" % string
        )
        raise argparse.ArgumentTypeError(msg)
    if (string[0] == "0") and (string[2] != "0"):
        msg = "%r is symmetric (digit 1) and must not have camber (digit 3)" % string
        raise argparse.ArgumentTypeError(msg)
    return string


# command line parser
parser = argparse.ArgumentParser(
    prog="naca5series",
    description="Create a NACA 5-series aerofoil",
    epilog="See: https://en.wikipedia.org/wiki/NACA_airfoil#Five-digit_series",
)

# parser group for aerofoil definition
group_def = parser.add_argument_group("Aerofoil definition")
group_def.add_argument(
    "LPQTT",
    metavar="LPQTT",
    type=valid_naca5,
    help="5 digits specifying LPQTT (e.g. 23012)",
)
group_def.add_argument(
    "-r",
    "--resolution",
    dest="res",
    metavar="N",
    type=int,
    default="100",
    choices=range(3, 1000),
    help="number of points on camber line (default: %(default)s)",
)
group_def.add_argument(
    "-c",
    "--chordlength",
    dest="chord",
    metavar="C",
    type=float,
    default="1.0",
    help="chord length (default: %(default)s)",
)
group_def.add_argument(
    "-s",
    "--spacing",
    dest="spacing",
    metavar="TYPE",
    type=str,
    default="cos",
    choices=(
        "cos",
        "2cos",
        "lin",
    ),
    help="spacing (default: %(default)s)",
)
group_def.add_argument(
    "-z",
    "--zcoordinate",
    dest="zval",
    metavar="Z",
    type=float,
    default="0.0",
    help="z position of aerofoil section (default: %(default)s)",
)
group_def.add_argument(
    "-p",
    "--plane",
    dest="plane",
    metavar="IJ",
    type=str,
    default="xy",
    choices=(
        "xy",
        "xz",
        "yz",
    ),
    help="aerofoil plane (default: %(default)s)",
)
group_def.add_argument(
    "--meancamberline",
    dest="meancamberline",
    action="store_true",
    help="only write the mean camber line (default: false)",
)

# parser group for input/output
group_IO = parser.add_argument_group("Output")
group_IO.add_argument(
    "outfile",
    metavar="FILE.csv",
    nargs="?",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="file to write to (default: stdout)",
)

# process
args = parser.parse_args()  # reads from sys.argv


### CALCULATE AEROFOIL ###


# auxiliary functions (to avoid using NumPy)
def linspace(a, b, N=100):
    if N == 1:
        return b
    h = (b - a) / (N - 1)
    return [a + h * n for n in range(N)]


def cos(x):
    # using https://docs.python.org/2/library/functions.html#map
    return map(math.cos, x)


def dealzeros(N, M):
    return [[0] * N for _ in range(M)]


# camberline x-coordinates
if args.spacing == "cos":
    x = [(1.0 - xi) / 2.0 for xi in cos(linspace(0.0, math.pi, args.res))]
elif args.spacing == "2cos":
    x = [1.0 - xi for xi in cos(linspace(0.0, math.pi / 2.0, args.res))]
elif args.spacing == "lin":
    x = linspace(0.0, 1.0, args.res)

# extract values from NACA definition
L = int(args.LPQTT[0])
P = int(args.LPQTT[1])
Q = args.LPQTT[2]
T = float(args.LPQTT[3:5]) / 100  # maximum thickness
cld = 0.15 * L  # design lift coefficient
p = 0.05 * P  # position of maximum camber
r, k1_base = MEANLINE_TABLE[(args.LPQTT[1], args.LPQTT[2])]
k1 = k1_base * (cld / 0.3)  # tabulated k1 values are for cld = 0.3

# calculate
Npts = len(x)
x_C, y_C, x_U, y_U, x_L, y_L = dealzeros(Npts, 6)  # initialise all
for i in range(Npts):
    x_c = x[i]

    # thickness (same as naca4series)
    y_t = (
        T
        / 0.20
        * (
            0.2969 * math.sqrt(x_c)
            - 0.1260 * x_c
            - 0.3516 * x_c**2
            + 0.2843 * x_c**3
            - 0.1015 * x_c**4
        )
    )

    # mean camber line
    if cld == 0:
        # no lift means symmetrical and thus no camber
        y_c = 0
        dycdx = 0
    elif Q == "0":
        # standard mean line
        if 0 <= x_c < r:
            y_c = k1 / 6.0 * (x_c**3 - 3 * r * x_c**2 + r**2 * (3 - r) * x_c)
            dycdx = k1 / 6.0 * (3 * x_c**2 - 6 * r * x_c + r**2 * (3 - r))
        else:
            y_c = k1 / 6.0 * r**3 * (1 - x_c)
            dycdx = -k1 / 6.0 * r**3
    else:
        # reflex mean line
        k2_k1 = (3 * (r - p) ** 2 - r**3) / ((1 - r) ** 3)  # from NACA-TR-537
        if 0 <= x_c < r:
            y_c = (
                k1
                / 6.0
                * ((x_c - r) ** 3 - k2_k1 * (1 - r) ** 3 * x_c - r**3 * x_c + r**3)
            )
            dycdx = k1 / 6.0 * (3 * (x_c - r) ** 2 - k2_k1 * (1 - r) ** 3 - r**3)
        else:
            y_c = (
                k1
                / 6.0
                * (
                    k2_k1 * (x_c - r) ** 3
                    - k2_k1 * (1 - r) ** 3 * x_c
                    - r**3 * x_c
                    + r**3
                )
            )
            dycdx = (
                k1 / 6.0 * (3 * k2_k1 * (x_c - r) ** 2 - k2_k1 * (1 - r) ** 3 - r**3)
            )

    # coordinates
    theta = math.atan(dycdx)
    # mean camber line
    x_C[i] = x_c
    y_C[i] = y_c
    # upper surface
    x_U[i] = x_c - y_t * math.sin(theta)
    y_U[i] = y_c + y_t * math.cos(theta)
    # lower surface
    x_L[i] = x_c + y_t * math.sin(theta)
    y_L[i] = y_c - y_t * math.cos(theta)


### EXPORT DATA ###

# prepare file
writer = csv.writer(
    args.outfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
)
# don't write a header row ["x","y","z"], to be universally readable!

for i in range(2 * Npts):  # 2x for both lower and upper
    if args.meancamberline:
        # only print the mean camber line (x_C and y_C)
        if i >= Npts:
            continue
        x_i = x_C[i]
        y_i = y_C[i]
    else:
        if i < Npts:
            # go backwards
            x_i = x_U[-1 - i]
            y_i = y_U[-1 - i]
        elif i == Npts:
            continue  # skip this one (duplicate LE)
            # NOTE: TE is not closed by definition
        else:
            # go forward
            x_i = x_L[i - Npts]
            y_i = y_L[i - Npts]

    # transform
    x_i *= args.chord
    y_i *= args.chord
    z_i = args.zval

    # assign to correct plane and write data point
    if args.plane == "xy":
        row = [x_i, y_i, z_i]
    elif args.plane == "xz":
        row = [x_i, z_i, y_i]
    elif args.plane == "yz":
        row = [z_i, x_i, y_i]
    writer.writerow(row)
