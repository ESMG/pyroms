class Grid_HYCOM(object):
    """
    Grid object for HYCOM
    """

    def __init__(self, lon_t, lat_t, lon_vert, lat_vert, mask_t, z_t, h, angle, name):

        self.name = name

        self.lon_t = lon_t
        self.lat_t = lat_t

        self.lon_vert = lon_vert
        self.lat_vert = lat_vert

        self.mask_t = mask_t

        self.z_t = z_t

        self.h = h

        self.angle = angle
