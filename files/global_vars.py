# Global variables for holding various settings and data.

incl_value = None
'''
Inclination value of the disk in degrees.
Type: float or None
'''

pa_value = None
'''
Position angle of the disk in degrees.
Type: float or None
'''

center_coords_value = [0, 0]
'''
Coordinates of the center of the disk selected by the user.
Type: list of two floats [x, y]
'''

subhd = None
'''
Header of the sub-region cutout of the FITS image.
Type: astropy.io.fits.Header or None
'''

header = None
'''
Header of the original FITS image.
Type: astropy.io.fits.Header or None
'''

f = None
'''
Interpolation function for the cutout image.
Type: scipy.interpolate.interp2d or None
'''

inclr = None
'''
Inclination of the disk in radians.
Type: float or None
'''

xdisk = None
'''
Array of x coordinates in the disk image, in arcseconds.
Type: numpy.ndarray or None
'''

ydisk = None
'''
Array of y coordinates in the disk image, in arcseconds.
Type: numpy.ndarray or None
'''

image_deproj = None
'''
Deprojected image data.
Type: numpy.ndarray or None
'''

disk_delt = None
'''
Angular step size in the disk image.
Type: float or None
'''

PAr = None
'''
Position angle of the disk in radians.
Type: float or None
'''

rout = None
'''
Outer radius of the disk in arcseconds.
Type: float or None
'''