import numpy as np
from osgeo import gdal
from osgeo import osr


class GeoTiffFileInfo():

    def __init__(self, reference_file, band_num=1):
        """ Get single raster band information from a reference file """
        filehandle = gdal.Open(reference_file)
        band = filehandle.GetRasterBand(band_num)
        band_data = band.ReadAsArray()
        self.geotransform = filehandle.GetGeoTransform()
        self.geoprojection = filehandle.GetProjection()
        self.shape = band_data.shape
        filehandle = None


def read_singleband_gtiff(filename, band_num=1, verbose=False):
    """ Simple function reading in a single band of values from a geotiff file

    Agrs:
        filename (str):     File path of geotiff
        band_num (int):     Specified band number

    Returns:
        numpy.ndarray: Numpy array of data
    """
    if verbose:
        print('... reading: ', filename)

    filehandle = gdal.Open(filename)
    band = filehandle.GetRasterBand(band_num)
    band_data = band.ReadAsArray()
    filehandle = None

    return np.array(band_data)


def write_singleband_gtiff(outfile, dat, reference_file, dtype=gdal.GDT_Float32):
    """ Writes out a single-band geotiff file, using geospatial information from a reference file

    Args:
        outfile (str):              filepath for output file
        dat (numpy.ndarray):        numpy array of data to write to file
        reference_file (str):       filepath of a reference file to set geospatial information
        dtype (osego.gdalconst):    data type to be written (see: https://naturalatlas.github.io/node-gdal/classes/Constants%20(GDT).html)
    """
    # Extract relevant information from reference file
    file_info = GeoTiffFileInfo(reference_file)

    # Initiate the new file
    output_raster = gdal.GetDriverByName('GTiff').Create(outfile, file_info.shape[1], file_info.shape[0], 1, dtype)
    output_raster.SetGeoTransform(file_info.geotransform)   # Specify coordinates based on reference file info
    output_raster.SetProjection(file_info.geoprojection)    # Set same projection as reference file
    output_raster.GetRasterBand(1).WriteArray(dat)          # Write out array to the raster
    output_raster.FlushCache()
    return
