import os
import ctypes

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../../src/libmcrt3d.so')
lib.new_Camera.restype = ctypes.c_void_p

class Camera:
    def __init__(self, G):
        self.obj = lib.new_Camera()

        lib.set_camera_grid(ctypes.c_void_p(self.obj), ctypes.c_void_p(G.obj))

    def make_image(self, image):
        lib.make_image(ctypes.c_void_p(self.obj), ctypes.c_void_p(image.obj))
