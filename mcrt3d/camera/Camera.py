from ..mcrt3d import lib

class Camera:
    def __init__(self, G, Q):
        self.obj = lib.new_Camera()

        lib.set_camera_grid(self.obj, G.obj)
        lib.set_camera_params(self.obj, Q.obj)

    def make_image(self, image):
        lib.make_image(self.obj, image.obj)
