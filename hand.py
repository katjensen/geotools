###################################################################
# Derivation of the HAND (Height Above the Nearest Drainage) index
#
# based on DEM and derived flood accumulation from HydroSHEDS database:
# https://hydrosheds.cr.usgs.gov/webappcontent/HydroSHEDS_TechDoc_v10.pdf
#
# Requirements: GDAL 1.8+, numpy,
#
###################################################################
import argparse
import datetime
import os
import sys

import numpy as np
from osgeo import gdal


class Hand:

    def __init__(self, input_dir, base_name, overwrite, verbose, max_height=2000):
        self.input_dir = input_dir
        self.output_dir = os.path.join(input_dir, "hand")
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.overwrite = overwrite
        self.verbose = verbose
        self.max_height = max_height

        # for a center pixel, going clockwise starting at row, col+1 and looking at all neighbors
        #  we look for a direction value that would drain into center pixel
        self.drain_direction_values = [16, 32, 64, 128, 1, 2, 4, 8]
        self.count = 0  # pixel set counter
        self.wp = 0     # water pixel counter

        # Define input filenames & read in data
        self.flow_dir_file = os.path.join(self.input_dir, base_name + "_flowdir.tif")
        self.flow_dir = self.read_dat(self.flow_dir_file, save_info=True)
        self.dem_file = os.path.join(self.input_dir, base_name + "_dem.tif")
        self.dem = self.read_dat(self.dem_file)
        self.drainage_file = os.path.join(self.dat_dir, base_name + "_drainage.tif")
        self.wat_mask = self.read_dat(self.drainage_file).astype(np.bool)
        self.hand_file = os.path.join(self.output_dir, base_name + "_hand.tif")

    def read_dat(self, filename, save_info=False):
        filehandle = gdal.Open(filename)
        if filehandle is None:
            print('ERROR: dat file no data: ', filename)
            sys.exit(-1)

        self.RasterXSize = filehandle.RasterXSize
        self.RasterYSize = filehandle.RasterYSize
        self.RasterCount = filehandle.RasterCount

        if save_info:
            self.dat_info = filehandle

        band = filehandle.GetRasterBand(1)
        data = band.ReadAsArray(0, 0, filehandle.RasterXSize, filehandle.RasterYSize)

        if self.verbose:
            print 'Size is ', filehandle.RasterXSize, 'x', filehandle.RasterYSize, 'x', filehandle.RasterCount, "Dtype= ", data.dtype
        filehandle = None
        return data

    def progress_bar(self, complete=0.0):
        gdal.TermProgress_nocb(complete)

    def process_hand(self):
        # Initilize hand array
        self.hand = np.full(self.flow_dir.shape, 0, dtype=np.int16)
        cols = self.RasterXSize
        rows = self.RasterYSize
        total = rows * cols
        ncount = 0

        # Number of water pixels
        num_wps = np.sum(self.wat_mask)
        if self.verbose:
            print str(datetime.datetime.now()), "Create Sorted Array of Water pixel data...", num_wps, total

        # Create sorted array of water pixels
        wph = np.empty((num_wps,), dtype=[('row', int), ('col', int), ('height', int)])

        for r in xrange(self.RasterYSize):
            for c in xrange(self.RasterXSize):
                if self.wat_mask[r, c]:
                    elev = self.dem[r, c]
                    self.hand[r, c] = 1     # water pixel is automatically set to 1 --> 1m or less
                    wph[self.wp] = (r, c, elev)
                    self.wp += 1
        wph.sort(order='height')

        if self.verbose:
            print str(datetime.datetime.now()), "Processing water pixel data ... ", num_wps

        if num_wps > 0:
            self.progress_bar(ncount / float(num_wps))

        # Get through all pixels in surface drainage reference
        for i in xrange(num_wps):
            wpx = wph[num_wps - i - 1]  # highest first
            r = wpx["row"]
            c = wpx["col"]
            ht = wpx["height"]

            flow_dir = self.flow_dir[r, c]
            if flow_dir == 247:     # oceans
                self.hand[r, c] = 255
            else:
                self.count += 1
                self.wp += 1
                neighbors = self.get_8_adjacent_neighbors(r, c)
                for idx, n in enumerate(neighbors):
                    if n is not None:
                        pixel = (n[0], n[1], idx)
                        self.process_pixel(pixel, ht)

            if self.verbose:
                ncount += 1
                self.progress_bar(ncount / float(num_wps))
        self.hand[self.hand == 0] = 0

    def get_8_adjacent_neighbors(self, row, col):
        """
        Neighboring elements are numbered like this:
        5  6  7
        4  x  0
        3  2  1

        (matching drainage flow direction value data
        """
        neighbors = [[row, col + 1], [row + 1, col + 1], [row + 1, col], [row + 1, col - 1], [row, col - 1],
                     [row - 1, col - 1], [row - 1, col], [row - 1, col + 1]]
        for i in xrange(8):
            neighbor_row, neighbor_col = neighbors[i]
            if (neighbor_row < 0) or (neighbor_row > self.RasterYSize - 1):
                neighbors[i] = None
            elif (neighbor_col < 0) or (neighbor_col > self.RasterXSize - 1):
                neighbors[i] = None

        return neighbors

    def process_pixel(self, neighbor, wp_height):
        """ Process current pixel and check the 8 adjacent pixels to find which drain into it """
        q = [neighbor]
        while len(q) > 0:
            neighbor = q.pop(0)
            row = neighbor[0]
            col = neighbor[1]
            index = neighbor[2]
            hand = self.hand[row, col]
            drain = self.flow_dir[row, col]

            if drain == 247:
                self.hand[row, col] = 255
                continue

            if hand > 0:
                continue

            if drain == self.drain_direction_values[index]:
                neighbor_height = self.dem[row, col]
                relative_height = neighbor_height - wp_height
                if relative_height <= self.max_height:
                    if relative_height < 0:
                        relative_height = 0
                    self.hand[row, col] = relative_height + 1
                    self.count += 1

                    neighbors = self.get_8_adjacent_neighbors(row, col)
                    for index, n in enumerate(neighbors):
                        if n is not None:
                            pixel = (n[0], n[1], index)
                            q.append(pixel)

    def write_out(self):
        print str(datetime.datetime.now()), "Saving Hand dat..."
        print " water pixels= ", self.wp
        print " processed= ", self.count
        print " total pixels= ", self.RasterXSize * self.RasterYSize
        print " hand pixels= ", np.count_nonzero(self.hand)

        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(self.hand_file, self.RasterXSize, self.RasterYSize, 1, gdal.GDT_Int16)
        band = dst_ds.GetRasterBand(1)
        band.WriteArray(self.hand, 0, 0)

        projection = self.dat_info.GetProjection()
        geotransform = self.dat_info.GetGeoTransform()

        dst_ds.SetGeoTransform(geotransform)
        dst_ds.SetProjection(projection)
        dst_ds.FlushCache()

        dst_ds = None
        self.dat_info = None

        return



# ======
# Main
#
# hand.py -i "my_input_directory" -b "s06w075"


if __name__ == "__main__":
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num < 1800: # because of GetGeoTransform(can_return_null)
        print('ERROR: Python bindings of GDAL 1.8.0 or later required')
        sys.exit(1)

    parser = argparse.ArgumentParser(description='Generate HAND')
    apg_input = parser.add_argument_group('Input')
    apg_input.add_argument("-i", "--inputdir", nargs=1,
                           help="Filepath to directory containing all 3 required input files")
    apg_input.add_argument("-b", "--basename", nargs=1,
                           help="Three required input files (flowdir, dem, drainage) must all have identical basenames "
                                "in each filename preceeding 'X.tif'")
    apg_input.add_argument("-m", "--maxheight", default=20, nargs=1,
                           help="Maximum height above nearest drainage to be considered [meters]")
    apg_input.add_argument("-o", "--overwrite", action='store_true',
                           help="Will overwrite any existing files found in hand directory")
    apg_input.add_argument("-v", "--verbose", action='store_true', help="Verbose on/off")
    options = parser.parse_args()

    H = Hand(input_dir=options.inputdir, base_name=options.base_name, max_height=int(options.maxheight))

    if not os.path.isfile(H.hand_file) or options.overwrite:
        H.process_hand()
        H.write_out()
    else:
        print("hand_img exists", H.hand_file)

    if options.verbose:
        print(str(datetime.now()), "Done.")