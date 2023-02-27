# NIST query
# https://astroquery.readthedocs.io/en/latest/nist/nist.html
# from astroquery.nist import Nist
# import astropy.units as u
# table = Nist.query(4000 * u.AA, 7000 * u.AA, linename="H I")


def search_line(wavelength, abs_err=10, element_spectrum=None):
    """Search line in Atomic Line List database.
    Parameters
    ----------
    wavelength : float
        Wavelength in Angstrom.
    abs_err : float
        Absolute error in Angstrom.
    Returns
    -------
    line : str
        Line name.
    """
    from astroquery.atomic import AtomicLineList
    import astropy.units as u
    wavelength_range = (wavelength - abs_err) * u.AA, (wavelength + abs_err) * u.AA
    table = AtomicLineList.query_object(wavelength_range=wavelength_range, wavelength_type='Vacuum', element_spectrum=element_spectrum)
    return table

def air_refractive_index(wavelength, ref='Morton+00', wave_type='vac'):
    """Air refractive index.
    Parameters
    ----------
    wavelength : float or array
        Wavelength in vacuum in units of Angstrom.
    ref : str
        Reference for the refractive index.
        Options:
            - 'Morton+00' (default)
    wave_type : str
        Type of wavelength.
        Options:
            - 'vac' (default)
            - 'air'
    Returns
    -------
    n : float or array
        Refractive index.

    References
    ----------
    https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf
    https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    """
    if wave_type == 'vac':
        if ref == 'Morton+00':
            # Morton (2000) ApJS 130, 403
            # https://ui.adsabs.harvard.edu/abs/2000ApJS..130..403M/abstract
            # Eq. 8
            sigma2 = (1e4 / wavelength) ** 2
            n = 1 + 1e-8 * (8342.13 + 2406030 / (130 - sigma2) + 15997 / (38.9 - sigma2))
        else:
            raise ValueError('Unknown reference for the refractive index.')

    elif wave_type == 'air':
        if ref == 'Morton+00':
            # reverse of Morton (2000)
            # https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
            sigma2 = (1e4 / wavelength) ** 2
            1 + 1e-8 * (8336.624212083 + 2408926.869968 / (130.1065924522 - sigma2) + 15997.40894897 / (38.92568793293 - sigma2))
        else:
            raise ValueError('Unknown reference for the refractive index.')
    return n



def vac_to_air(wavelength):
    """Convert vacuum wavelength to air wavelength.
    Parameters
    ----------
    wavelength : float or array
        Wavelength in vacuum.
    Returns
    -------
    wavelength : float or array
        Wavelength in air.
    """
    refact_index = air_refractive_index(wavelength, ref='Morton+00', wave_type='vac')
    return wavelength / refact_index

def air_to_vac(wavelength):
    """Convert air wavelength to vacuum wavelength.
    Parameters
    ----------
    wavelength : float or array
        Wavelength in air.
    Returns
    -------
    wavelength : float or array
        Wavelength in vacuum.
    """
    refact_index = air_refractive_index(wavelength, ref='Morton+00', wave_type='air')
    return wavelength * refact_index


# # Commonly used lines
class AstroLine:
    def __init__(self, name, wavelength):
        self.name = name
        self.wavelength = wavelength