import warnings
# warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u

# From grizli
def get_wcs_pscale(wcs, unit=u.arcsec):
    """Get correct pscale from a `~astropy.wcs.WCS` object
    Parameters
    ----------
    wcs : `~astropy.wcs.WCS` or `~astropy.io.fits.Header`
    unit: `~astropy.units.Unit`
    Returns
    -------
    pscale : float
        Pixel scale
    """
    pscale = proj_plane_pixel_scales(wcs)[0] * 3600 * unit
    return pscale.to(unit).value

def transform_wcs(in_wcs, translation=[0., 0.], rotation=0., scale=1.):
    """Update WCS with shift, rotation, & scale
    Parameters
    ----------
    in_wcs: `~astropy.wcs.WCS`
        Input WCS
    translation: [float, float]
        xshift & yshift in pixels
    rotation: float
        CCW rotation (towards East), radians
    scale: float
        Pixel scale factor
    Returns
    -------
    out_wcs: `~astropy.wcs.WCS`
        Modified WCS
    """
    out_wcs = in_wcs.deepcopy()

    # Compute shift for crval
    crval = in_wcs.all_pix2world([in_wcs.wcs.crpix-np.array(translation)],
                                    1).flatten()

    # Compute shift at image center
    if hasattr(in_wcs, '_naxis1'):
        refpix = np.array([in_wcs._naxis1/2., in_wcs._naxis2/2.])
    else:
        refpix = np.array(in_wcs._naxis)/2.

    c0 = in_wcs.all_pix2world([refpix], 1).flatten()
    c1 = in_wcs.all_pix2world([refpix-np.array(translation)], 1).flatten()

    out_wcs.wcs.crval += c1-c0

    theta = -rotation
    _mat = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])

    out_wcs.wcs.pc = np.dot(out_wcs.wcs.pc, _mat)/scale

    out_wcs.pscale = get_wcs_pscale(out_wcs)
    if hasattr(out_wcs, 'pixel_shape'):
        _naxis1 = int(np.round(out_wcs.pixel_shape[0]*scale))
        _naxis2 = int(np.round(out_wcs.pixel_shape[1]*scale))
        out_wcs._naxis = [_naxis1, _naxis2]
    elif hasattr(out_wcs, '_naxis1'):
        out_wcs._naxis1 = int(np.round(out_wcs._naxis1*scale))
        out_wcs._naxis2 = int(np.round(out_wcs._naxis2*scale))

    return out_wcs
