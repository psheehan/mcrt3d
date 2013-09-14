#!/usr/bin/env python

from numpy import *
from scipy.integrate import trapz
from misc import B_nu, dB_nu
import matplotlib.pyplot as plt

data = loadtxt("dustkappa_yso.inp",usecols=(0,1,2),skiprows=2)

lam = data[:,0]
kabs = data[:,1]
ksca = data[:,2]

kext = kabs+ksca
albedo = ksca/kext

c_l = 2.99792458e10
sigma = 5.67051e-5

nu = (c_l/(lam*1.0e-4)).astype(float)

T = logspace(-1,4,100).astype(float)

int_Bnu = ones(T.size)
int_Bnukext = ones(T.size)
int_dBnukext = ones(T.size)
for i in range(T.size):
    Bnu = B_nu(nu,T[i])
    dBnu = dB_nu(nu,T[i])

    int_Bnu[i] = trapz(Bnu,x=nu)
    int_Bnukext[i] = trapz(Bnu*kext,x=nu)
    int_dBnukext[i] = trapz(dBnu*kext,x=nu)

planck_opacity = -int_Bnukext/(sigma*T**4/pi)

dplanck_opacity_dT = zeros(T.size)
dint_dBnukext_dT = zeros(T.size)

for i in range(T.size-1):
    dplanck_opacity_dT[i] = (planck_opacity[i+1]-planck_opacity[i])/ \
            (T[i+1]-T[i])
    dint_dBnukext_dT[i] = (int_dBnukext[i+1]-int_dBnukext[i])/ \
            (T[i+1]-T[i])

dkextdnu = zeros(nu.size)
dalbedodnu = zeros(nu.size)
for i in range(nu.size-1):
    dkextdnu[i] = (kext[i+1]-kext[i])/(nu[i+1]-nu[i])
    dalbedodnu[i] = (albedo[i+1]-albedo[i])/(nu[i+1]-nu[i])

f = open("lookup_table1.txt","w")
f.write("{0:4d}\n".format(T.size))
for i in range(T.size):
    f.write("{0:e}   {1:e}   {2:e}   {3:e}   {4:e}\n".format(T[i], \
            planck_opacity[i], int_dBnukext[i], dplanck_opacity_dT[i], \
            dint_dBnukext_dT[i]))
f.close()

f = open("lookup_table2.txt","w")
for i in range(T.size):
    Bnu = B_nu(nu,T[i])
    dBnu = dB_nu(nu,T[i])
    
    if (i < T.size-1):
        dBnudT = (B_nu(nu,T[i+1])-B_nu(nu,T[i]))/(T[i+1]-T[i])
        ddBnudT = (dB_nu(nu,T[i+1])-dB_nu(nu,T[i]))/(T[i+1]-T[i])

    for j in range(nu.size):
        test1 = trapz(Bnu[0:j]*kext[0:j],x=nu[0:j])
        test2 = trapz(dBnu[0:j]*kext[0:j],x=nu[0:j])
        f.write("{0:e}   {1:e}   {2:e}   {3:e}   {4:e}   {5:e}\n".format( \
                Bnu[j], dBnu[j], dBnudT[j], ddBnudT[j], test1, test2))
f.close()

f = open("lookup_table3.txt","w")
for i in range(nu.size):
    f.write("{0:e}   {1:e}\n".format(dkextdnu[i], dalbedodnu[i]))
f.close()
