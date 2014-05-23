#!/usr/bin/env python3

from mcrt3d import *
from time import time
import matplotlib.pyplot as plt

from Cartesian import *
#from Cylindrical import *
#from Spherical import *
#from YSO import *

G = SetupGrid()
#images = SetupImages()
#spectra = SetupSpectra()

# Now give that all to C.

M = MCRT(G)
C = Camera(G)

t1 = time()
M.thermal_mc(1000000,True)
t2 = time()
print(t2-t1)

"""
for i in range(images.size):
    C.make_image(images[i])

for i in range(spectra.size):
    C.make_image(spectra[i])
"""

for i in range(9):
    plt.imshow(G.temp[:,:,i],origin="lower",interpolation="nearest", \
            vmin=G.temp.min(),vmax=G.temp.max())
    plt.colorbar()
    plt.show()

"""
plt.imshow(images[0].intensity[:,:,0],origin="lower",interpolation="nearest")
plt.show()

plt.loglog(c_l / spectra[0].nu, spectra[0].intensity[0,0,:])
plt.axis([1e-5,1e0,1e-19,1e-9])
plt.show()
"""
