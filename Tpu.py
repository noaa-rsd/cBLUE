"""
cBLUE (comprehensive Bathymetric Lidar Uncertainty Estimator)
Copyright (C) 2019 
Oregon State University (OSU)
Center for Coastal and Ocean Mapping/Joint Hydrographic Center, University of New Hampshire (CCOM/JHC, UNH)
NOAA Remote Sensing Division (NOAA RSD)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Contact:
Christopher Parrish, PhD
School of Construction and Civil Engineering
101 Kearney Hall
Oregon State University
Corvallis, OR  97331
(541) 737-5688
christopher.parrish@oregonstate.edu

"""

import logging
from pathlib import Path
from pathos import logger
import pathos.pools as pp
import json
import os
import sys
import laspy
import numpy as np
import numexpr as ne
import pandas as pd
import progressbar
from tqdm import tqdm
from collections import OrderedDict
from Subaerial import Subaerial, SensorModel, Jacobian
from Subaqueous import Subaqueous
from Las import Las
import datetime
import cProfile
import pathos
import PySimpleGUI as sg
from functools import partial

from Merge import Merge


class Tpu:
    """
    TODO:  rework...becasue J & M moved to CBlueApp.py
    This class coordinates the TPU workflow.  Beginning when the user 
    hits *Compute TPU*, the general workflow is summarized below:

    1. Form observation equation (SensorModel class)
    2. Generate Jacobian (Jacobian class)
    3. for each flight line within Las

        * Merge the Las data and trajectory data (Merge class)
        * Calculate subaerial thu and tvu (Subaerial class)
        * Calculate subaqueous thu and tvu (Subaqueous class)
        * Combine subaerial and subaqueous TPU
        * Export TPU (either as Python "pickle' or as Las extra bytes)
    
    """

    def __init__(self, surface_select, surface_ind,
                 wind_selection, wind_val,
                 kd_selection, kd_val, vdatum_region,
                 vdatum_region_mcu, tpu_output, cblue_version, 
                 sensor_model, cpu_process_info, subaqueous_luts,
                 water_surface_ellipsoid_height):

        # TODO: refactor to pass the GUI object, not individual variables (with controller?)
        self.surface_select = surface_select
        self.surface_ind = surface_ind
        self.wind_selection = wind_selection
        self.wind_val = wind_val
        self.kd_selection = kd_selection
        self.kd_val = kd_val
        self.vdatum_region = vdatum_region
        self.vdatum_region_mcu = vdatum_region_mcu
        self.tpu_output = tpu_output
        self.cblue_version = cblue_version
        self.sensor_model = sensor_model
        self.cpu_process_info = cpu_process_info
        self.subaqueous_luts = subaqueous_luts
        self.water_surface_ellipsoid_height = water_surface_ellipsoid_height

        self.subaqu_lookup_params = None
        self.metadata = {}
        self.flight_line_stats = {}

    def update_fl_stats(self, fl, num_fl_points, fl_tpu_data):

        # calc flight line tpu summary stats
        fl_tpu_count = fl_tpu_data.shape[0]
        fl_tpu_min = fl_tpu_data[:, 0:6].min(axis=0).tolist()
        fl_tpu_max = fl_tpu_data[:, 0:6].max(axis=0).tolist()
        fl_tpu_mean = fl_tpu_data[:, 0:6].mean(axis=0).tolist()
        fl_tpu_stddev = fl_tpu_data[:, 0:6].std(axis=0).tolist()

        fl_stat_indx = {
            'total_thu': 0,
            'total_tvu': 1,
            }

        fl_stats_strs = []
        for fl_stat, ind in fl_stat_indx.items():
            fl_stats_vals = (fl_stat,
                             fl_tpu_min[ind], 
                             fl_tpu_max[ind],
                             fl_tpu_mean[ind], 
                             fl_tpu_stddev[ind])
        
            fl_stats_str = '{}: {:6.3f}{:6.3f}{:6.3f}{:6.3f}'.format(*fl_stats_vals)
            fl_stats_strs.append(fl_stats_str)

        fl_header_str = f'{fl} ({fl_tpu_count}/{num_fl_points} points with TPU)'
        self.flight_line_stats.update({fl_header_str: fl_stats_strs})

    def calc_tpu(self, sbetC, sensor_model_name, las_file):
        """

        :param sbetC: sbet object passsed from CBlueApp
        :param sensor_model_name: string name representation of chosen sensor model
        :param las_file: the las file we are calculating tpu for
        :return:
        """
        sbet = self.sbet_get_las_tile(sbetC, las_file)

        merge = Merge()

        # CREATE OBSERVATION EQUATIONS
        sensor_model = SensorModel(sensor_model_name)

        jacobian = Jacobian(sensor_model)

        data_to_output = []

        # CREATE LAS OBJECT TO ACCESS INFORMATION IN LAS FILE
        las = Las(las_file)
        
        diag_dfs = []

        if las.num_file_points:  # i.e., if las had data points
            logging.info('{} ({:,} points)'.format(las.las_short_name, las.num_file_points))
            logging.debug('flight lines {}'.format(las.unq_flight_lines))
            unsorted_las_xyzt, t_argsort, flight_lines = las.get_flight_line_txyz()

            self.flight_line_stats = {}  # reset flight line stats dict
            for fl in las.unq_flight_lines:

                logging.debug('flight line {} {}'.format(fl, '-' * 50))

                # las_xyzt has the same order as points in las (i.e., unordered)
                fl_idx = flight_lines == fl
                fl_unsorted_las_xyzt = unsorted_las_xyzt[fl_idx]
                fl_t_argsort = t_argsort[fl_idx]
                fl_las_idx = t_argsort.argsort()[fl_idx]

                num_fl_points = np.sum(fl_idx)  # count Trues
                logging.debug(f'{las.las_short_name} fl {fl}: {num_fl_points} points')

                # CREATE MERGED-DATA OBJECT M
                logging.debug('({}) merging trajectory and las data...'.format(las.las_short_name))
                merged_data, stddev, unsort_idx, raw_class = merge.merge(las.las_short_name, fl, sbet.values,
                                                                   fl_unsorted_las_xyzt, fl_t_argsort, fl_las_idx)

                if merged_data is not False:  # i.e., las and sbet is merged

                    logging.debug('({}) calculating subaer thu/tvu...'.format(las.las_short_name))
                    subaer_obj = Subaerial(jacobian, merged_data, stddev)
                    subaer_thu, subaer_tvu = subaer_obj.calc_subaerial_tpu()
                    depth = merged_data[4] - self.water_surface_ellipsoid_height

                    logging.debug('({}) calculating subaqueous thu/tvu...'.format(las.las_short_name))
                    subaqu_obj = Subaqueous(self.surface_ind, self.wind_val, 
                                            self.kd_val, depth, self.subaqueous_luts)
                    subaqu_thu, subaqu_tvu = subaqu_obj.fit_lut()
                    self.subaqu_lookup_params = subaqu_obj.get_subaqueous_meta_data()
                    vdatum_mcu = float(self.vdatum_region_mcu) / 100.0  # file is in cm (1-sigma)

                    logging.debug('({}) calculating total thu...'.format(las.las_short_name))
                    total_thu = ne.evaluate('sqrt(subaer_thu**2 + subaqu_thu**2)')

                    logging.debug('({}) calculating total tvu...'.format(las.las_short_name))
                    total_tvu = ne.evaluate('sqrt(subaer_tvu**2 + subaqu_tvu**2 + vdatum_mcu**2)')

                    fl_tpu_data = np.vstack((
                        total_thu,
                        total_tvu,
                        unsort_idx
                        )).T

                    data_to_output.append(fl_tpu_data)

                    self.update_fl_stats(fl, num_fl_points, fl_tpu_data)

                else:
                    logging.warning('SBET and LAS not merged because max delta '
                                    'time exceeded acceptable threshold of {} '
                                    'sec(s).'.format(merge.max_allowable_dt))

                    self.flight_line_stats.update(
                        {'{} (0/{} points with TPU)'.format(fl, num_fl_points): None})

            self.write_metadata(las)  # TODO: include as VLR?
            self.output_tpu_to_las_extra_bytes(las, data_to_output)

        else:
            logging.warning('WARNING: {} has no data points'.format(las.las_short_name))

    def output_tpu_to_las_extra_bytes(self, las, data_to_output):
        """output the calculated tpu to a las file

        This method creates a las file tht contains the contents of the
        original las file and the calculated tpu values as VLR extra bytes.
        The las file is generated using "The laspy way", as documented in
        https://laspy.readthedocs.io/en/latest/tut_part_3.html.

        The following references have additional information describing las
        extra bytes:

        LAS v1.4 specifications:
        https://www.asprs.org/a/society/committees/standards/LAS_1_4_r13.pdf

        The LAS 1.4 Specification (ASPRS PERS article)
        https://www.asprs.org/wp-content/uploads/2010/12/LAS_Specification.pdf

        ASPRS LAS Working Group Github repository
        https://github.com/ASPRSorg/LAS

        The following table lists the information contained as extra bytes:

        .. csv-table:: cBLUE VLR Extra Bytes
            :header: id, dtype, description
            :widths: 14, 20, 20

            total_thu,  unsigned short (2 bytes), total horizontal uncertainty
            total_tvu,  unsigned short (2 bytes), total vertical uncertainty

        :param las:
        :param data_to_output:
        :param output_columns:
        :return:
        """

        out_las_name = os.path.join(self.tpu_output, las.las_base_name) + '_TPU.las'
        logging.debug('logging las and tpu results to {}'.format(out_las_name))
        in_las = laspy.file.File(las.las, mode="r")  # las is Las object
        out_las = laspy.file.File(out_las_name, mode="w", header=in_las.header)

        tpu_data_type = 9  # float

        extra_byte_dimensions = OrderedDict([
            ('total_thu', tpu_data_type),
            ('total_tvu', tpu_data_type)
            ])

        num_extra_bytes = len(extra_byte_dimensions.keys())

        # define new extrabyte dimensions
        for dimension, dtype in extra_byte_dimensions.items():
            logging.debug('creating extra byte dimension for {}...'.format(dimension))
            out_las.define_new_dimension(name=dimension, 
                                         data_type=dtype,
                                         description=dimension)

        if len(data_to_output) != 0:
            tpu_data = np.vstack(data_to_output)
            extra_byte_df = pd.DataFrame(tpu_data[:, 0:num_extra_bytes], 
                                         index=tpu_data[:, num_extra_bytes], 
                                         columns=extra_byte_dimensions.keys())

            print(f'extra_byte_df')
            print(extra_byte_df)

            print(f'las.num_file_points')
            print(las.num_file_points)

            print('las.t_argsort')
            print(las.t_argsort)

            if extra_byte_df.shape[0] == las.num_file_points:
                extra_byte_df = extra_byte_df.sort_index()
                print(f'extra_byte_df --------- sort.index()')
                print(extra_byte_df)
            else:
                logging.info("""
                filling data points for which TPU was not calculated 
                with no_data_value (also sorting by index, or t_idx)
                """)
                no_data_value = -1
                extra_byte_df = extra_byte_df.reindex(las.t_argsort,
                                                      fill_value=no_data_value).sort_index()

            logging.debug('populating extra byte data for total_thu...')
            out_las.total_thu = extra_byte_df['total_thu']
        
            logging.debug('populating extra byte data for total_tvu...')
            out_las.total_tvu = extra_byte_df['total_tvu']
        
        else:
            logging.debug('populating extra byte data for total_thu...')
            out_las.total_thu = np.zeros(las.num_file_points)
        
            logging.debug('populating extra byte data for total_tvu...')
            out_las.total_tvu = np.zeros(las.num_file_points)

        # copy data from in_las
        for field in in_las.point_format:
            logging.debug('writing {} to {} ...'.format(field.name, out_las))
            las_data = in_las.reader.get_dimension(field.name)
            out_las.writer.set_dimension(field.name, las_data[las.t_argsort])

    def write_metadata(self, las):
        """creates a json file with summary statistics and metedata

        This method creates a json file containing summary statistics for each
        tpu field, per flight line, and a record of the environmental and VDatum
        parameters specified by the user.  The file also records certain
        parameters used during the monte carlo simulations used to create
        the lookup tables used in the subaqueous portion of the tpu calculations.

        :param las:
        :return: n/a
        """

        logging.debug('({}) creating TPU meta data file...'.format(las.las_short_name))
        self.metadata.update({
            'subaqueous lookup params': self.subaqu_lookup_params,
            'water surface': self.surface_select,
            'wind': self.wind_selection,
            'kd': self.kd_selection,
            'VDatum region': self.vdatum_region,
            'VDatum region MCU': self.vdatum_region_mcu,
            'flight line stats (min max mean stddev)': self.flight_line_stats,
            'sensor model': self.sensor_model,
            'cBLUE version': self.cblue_version,
            'cpu_processing_info': self.cpu_process_info, 
            'water_surface_ellipsoid_height': self.water_surface_ellipsoid_height,
        })

        try:
            #self.metadata['flight line stats'].update(self.flight_line_stats)  # flight line metadata
            with open(os.path.join(self.tpu_output, '{}.json'.format(las.las_base_name)), 'w') as outfile:
                json.dump(self.metadata, outfile, indent=1, ensure_ascii=False)
        except Exception as e:
            logging.error(e)
            print(e)

    def sbet_las_tiles_generator(self, sbet, las_files, jacobian, merge):
        """This generator is the 2nd argument for the
        run_tpu_multiprocessing method, to avoid
        passing entire sbet or list of tiled
        sbets to the calc_tpu() method
        """

        for las_file in las_files:
            logging.debug('({}) generating SBET tile...'.format(las_file.split('\\')[-1]))

            inFile = laspy.file.File(las_file, mode='r')

            west = inFile.reader.get_header_property('x_min')
            east = inFile.reader.get_header_property('x_max')
            north = inFile.reader.get_header_property('y_max')
            south = inFile.reader.get_header_property('y_min')

            yield sbet.get_tile_data(north, south, east, west), las_file, jacobian, merge

    def sbet_get_las_tile(self, sbet, las_file):
        """Replaces generator function for multiprocessing, exe version
        """

        logging.debug('({}) generating SBET tile...'.format(las_file.split('\\')[-1]))

        inFile = laspy.file.File(las_file, mode='r')

        west = inFile.reader.get_header_property('x_min')
        east = inFile.reader.get_header_property('x_max')
        north = inFile.reader.get_header_property('y_max')
        south = inFile.reader.get_header_property('y_min')

        return sbet.get_tile_data(north, south, east, west)


    def run_tpu_multiprocess(self, num_las, sbet, las_files, sensor_model):
        """runs the tpu calculations using multiprocessing

        This methods initiates the tpu calculations using the pathos
        multiprocessing framework (https://pypi.org/project/pathos/).
        Whether the tpu calculations are done with multiprocessing or not is
        currently determined by which "run_tpu_*" method is manually specified
        in the tpu_process_callback() method of the CBlueApp class.

        TODO: Include user option to select single processing or multiprocessing

        :param num_las:
        :param sbet: sbet object passsed from CBlueApp
        :param las_files: list of las files to process
        :param sensor_model: string name representation of chosen sensor model
        :return:
        """

        print('Calculating TPU (multi-processing)...')
        p = pp.ProcessPool(2)

        #for _ in tqdm(p.imap(self.calc_tpu, self.sbet_las_tiles_generator(sbet, las_files, jacobian, merge)), total=num_las, ascii=True):
        #    pass

        func = partial(self.calc_tpu, sbet, sensor_model)

        for _ in tqdm(p.imap(func, las_files), total=num_las, ascii=True):
            pass

        # Initialize progress bar so that it shows from onset
        # sg.OneLineProgressMeter('Calculating TPU (multi-processing)...', 0, num_las, "pbar")
        # for num in enumerate(p.imap(func, las_files)):
        #    sg.OneLineProgressMeter('Calculating TPU (multi-processing)...', num[0] + 1, num_las, "pbar")

        return p

    def run_tpu_singleprocess(self, num_las, sbet_las_generator):
        """runs the tpu calculations using a single processing

        This methods initiates the tpu calculations using single processing.
        Whether the tpu calculations are done with multiprocessing or not is
        currently determined by which "run_tpu_*" method is manually specified
        the tpu_process_callback() method of the CBlueApp class.  TODO: Include
        a user option to select single processing or multiprocessing

        :param sbet_las_generator:
        :return:
        """

        print('Calculating TPU (single-processing)...')
        with progressbar.ProgressBar(max_value=num_las) as bar:
            for i, sbet_las in enumerate(sbet_las_generator):
                bar.update(i)
                self.calc_tpu(sbet_las)


if __name__ == '__main__':
    # Required for pyinstaller support of multiprocessing
    pathos.helpers.freeze_support()

    pass
