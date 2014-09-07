import numpy
from numpy import exp
from ..constants.physics import h, c, k

def dB_nu(nu,T):
    
    numpy.seterr(over="ignore")
    return (-2*h**2*nu**4/(c**2*k*T**2))/(exp(h*nu/(k*T))-1)/ \
            (1.-exp(-h*nu/(k*T)))
