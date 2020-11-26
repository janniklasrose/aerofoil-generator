""" inspiration:
https://github.com/rafael-aero/XFOILinterface/blob/master/%40Airfoil/createNACA5.m
https://github.com/dgorissen/naca/blob/master/naca.py

Documentation:
NACA-TR-537 : https://ntrs.nasa.gov/citations/19930091610
NACA-TR-610 : https://ntrs.nasa.gov/citations/19930091685
"""
def standard(x, k1, r):
    if 0 <= x and x < r:
        y_c = k1/6.0 * (x**3 - 3*r*x**2 + r**2*(3 - r)*x)
    elif m <= x <= 1:
        y_c = k1/6.0 * r**3 * (1 - x)
    else:
        y_c = None  #TODO: or NAN?
    return y_c

def reflex(x, k1, k2_k1, r):
    k2_k1 = (3*(m-p)**2-m**3)/((1-m)**3)  # see NACA Rept No 537
    if 0 <= x and x < r:
        y_c = k1/6 * ((x-r)**3 - k2_k1 * (1-r)**3*x - r**3*x + r**3)
    elif r <= x and x <= 1:
        y_c = k1/6 * (k2_k1*(x-r)**3 - k2_k1*(1-r)**3 - r**3*x + r**3)
    else:
        y_c = None  #TODO: or NAN?
    return y_c

# table[designation] = (r, k1, fn, descr)
#TODO: all for first digit = 2, i.e. Cl=0.3. second digit gives p, not all ranges given (only "useful range")
meanline_table = {}
meanline_table['210'] = (0.0580, 361.4  , standard, "05% standard")
meanline_table['220'] = (0.1260,  51.64 , standard, "10% standard")
meanline_table['230'] = (0.2025,  15.957, standard, "15% standard")
meanline_table['240'] = (0.2900,   6.643, standard, "20% standard")
meanline_table['250'] = (0.3910,   3.230, standard, "25% standard")
#### these reflex foils are almost never used, so not discussed further
# '211' was not reported in the NACA report due to an error, so we have no data on this
meanline_table['221'] = (0.1300,  51.990, reflex  , "10% reflex")
meanline_table['231'] = (0.2170,  15.793, reflex  , "15% reflex")
meanline_table['241'] = (0.3180,   6.520, reflex  , "20% reflex")
meanline_table['251'] = (0.4410,   3.191, reflex  , "25% reflex")

def parse(fivedigits):
    p = 0.1*3/2*int(fivedigits[0])
    x_cmax = 5*int(fivedigits[1])
    definition = meanline_table[fivedigits[0:3]]
    fn = definition[4]

    N = 20
    x = [0+i/(N-1) for i in range(N)]
    y = [0]*N
    for i in range(N):
        x_i = x[i]
        if fivedigits[2] is '0':
            y_i = fn(x_i, definition[1], definition[2])
        elif fivedigits[2] is '1':
            y_i = fn(x_i, definition[1], definition[2], definition[3])
        else:
            print("ERROR")
        y_i *= cld/0.3  # scale wrt default 0.3
        y[i] = y_i  # store
    return x, y

x, y = parse('23100')
print(x)
print(y)
