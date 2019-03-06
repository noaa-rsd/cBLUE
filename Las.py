import logging
import numpy as np
import numexpr as ne
import laspy

logging.basicConfig(format='%(asctime)s:%(message)s', level=logging.INFO)

"""
This class provides the functionality to load las files into cBLUE.  One Las object 
is created for each loaded las file.  
"""


class Las:

    def __init__(self, las):
        self.las = las
        self.las_short_name = las.split('\\')[-1]
        self.las_base_name = self.las_short_name.replace('.las', '')
        self.inFile = laspy.file.File(self.las, mode="r")
        self.num_file_points = self.inFile.__len__()
        self.points_to_process = self.inFile.points['point']

        '''index list that would sort gps_time (to be used to
        later when exporting las data and calculated tpu to a las
        file
        '''
        self.time_sort_indices = None

    def get_flight_line_ids(self):
        """generates a list of unique flight line ids

        This method returns a list of unique flight line ids.

        :return: ndarray
        """
        return np.unique(self.points_to_process['pt_src_id'])

    def get_flight_line_txyz(self):
        """retrieves the x, y, z, and timestamp data from the las data points

        The x, y, and z values in the las file are stored as integers.  The
        scale and offset values in the las file header are used to convert
        the integer values to decimal values with centimeter precision.

        :param ? fl:
        :return: np.array, np.array, np.array, np.array
        """
        scale_x = np.asarray(self.inFile.header.scale[0])
        scale_y = np.asarray(self.inFile.header.scale[1])
        scale_z = np.asarray(self.inFile.header.scale[2])

        offset_x = np.asarray(self.inFile.header.offset[0])
        offset_y = np.asarray(self.inFile.header.offset[1])
        offset_z = np.asarray(self.inFile.header.offset[2])

        t = self.points_to_process['gps_time']
        X = self.points_to_process['X']
        Y = self.points_to_process['Y']
        Z = self.points_to_process['Z']

        x = ne.evaluate("X * scale_x + offset_x")
        y = ne.evaluate("Y * scale_y + offset_y")
        z = ne.evaluate("Z * scale_z + offset_z")

        xyzt = np.vstack([x, y, z, t]).T
        self.time_sort_indices = xyzt[:,3].argsort()

        flight_lines = self.points_to_process['pt_src_id']

        return xyzt, self.time_sort_indices, flight_lines

    @staticmethod
    def get_average_water_surface_ellip_height():
        """returns the average ellipsoid height of the water surface returns

        Currently, this method returns a visually-determined estimate of the average ellipsoid height
        of the water surface returns in the survey area, which is used during tpu
        calculation to calculate the depth of each data point.

        # TODO: define better way to determine average ellipsoid height of surface

        :return: float
        """
        return 23.0


if __name__ == '__main__':
    pass
