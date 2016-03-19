''' PSF model functions

Author: Simon Zieleniewski

Last updated: 19-11-15

'''

import numpy as n
from scipy.special import j1


def obsc_airy(x, h, w, q):
    #Search for 0 in array - if not present, pass
    try:
        #loc = n.argwhere(x==0.)[0]
        #x[loc[0],loc[1]] = 1.E-30 #Fix for 0. in array - using 1.E-30 seems to give good result for Airy function
        loc = n.argwhere(x==0.)
        x[loc] = 1.E-30
    except IndexError:
        pass
    return h*((2*j1(n.sin(x)/w) - q*j1(q*(n.sin(x)/w)))/((n.sin(x)/w)*(1 - q**2)))**2


def moffat(x, h, w, q):
    return h/(1 + (x/w)**2)**(q)


def lorentz(x, h, p, w):
    return h/(1 + (((x-p)/w)**2))


def x1(x, a, b):
    return a + b*x


def x2(x, a, b, c):
    return a + b*x + c*x**2


def x6(x, a, b, c, d, e, f, g):
    return a +b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6
