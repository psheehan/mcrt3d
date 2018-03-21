#!/usr/bin/env python3

from mcrt3d import MCRT
from mcrt3d.camera import Camera
from time import time
import matplotlib.pyplot as plt
from mcrt3d.constants.physics import c as c_l

#from Cartesian import SetupParams, SetupGrid, SetupImages, SetupSpectra
#from Cylindrical import SetupParams, SetupGrid, SetupImages, SetupSpectra
from Spherical import SetupParams, SetupGrid, SetupImages, SetupSpectra
#from YSO import *

model = MCRT()

SetupParams(model)
SetupGrid(model)
images = SetupImages()
spectra = SetupSpectra()

# Now give that all to C.

t1 = time()
model.run_thermal_mc()
t2 = time()
print(t2-t1)

for i in range(images.size):
    model.run_image(images[i])

for i in range(spectra.size):
    model.run_image(spectra[i])

for i in range(9):
    plt.imshow(model.grid.temperature[0][:,:,i], origin="lower",\
            interpolation="nearest", vmin=model.grid.temperature[0].min(),\
            vmax=model.grid.temperature[0].max())
    plt.colorbar()
    plt.show()

plt.imshow(images[0].intensity[:,:,0],origin="lower",interpolation="none")
plt.show()

plt.loglog(c_l / spectra[0].nu*1.0e4, spectra[0].intensity[0,0,:])
plt.axis([1e-1,1e4,1e-19,1e-9])
plt.show()
