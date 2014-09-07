import numpy
from numpy import exp
from ..constants.physics import h, c, k

def B_nu(nu,T):
    
    numpy.seterr(over="ignore")
    return (2*h*nu**3/c**2)/(exp(h*nu/(k*T))-1.0)
