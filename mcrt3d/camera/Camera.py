import os
import ctypes

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../../src/libmcrt3d.so')

lib.new_Camera.restype = ctypes.c_void_p
lib.new_Camera.argtypes = None

lib.set_camera_grid.restype = None
lib.set_camera_grid.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.set_camera_params.restype = None
lib.set_camera_params.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.make_image.restype = None
lib.make_image.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

class Camera:
    def __init__(self, G, Q):
        self.obj = lib.new_Camera()

        lib.set_camera_grid(self.obj, G.obj)
        lib.set_camera_params(self.obj, Q.obj)

    def make_image(self, image):
        lib.make_image(self.obj, image.obj)
