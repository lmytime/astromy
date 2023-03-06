from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack
import numpy as np


SPEC_EDGE_WAVELENGTH = {"F115W": [1.0, 1.3], "F150W": [1.3, 1.7], "F200W": [1.7, 2.3]} # microns
NIRISS_PIXEL_SCALE = 0.065 # arcsec/pixel




class Grizli1D:
    def __init__(self, path):
        with fits.open(path) as hdul:
            self.hdul = hdul
            self.primary_header = hdul[0].header
            self.id = self.primary_header['ID']
            self.ra = self.primary_header['RA']
            self.dec = self.primary_header['DEC']
            self.target = self.primary_header['TARGET']
            self.data = {}
            for ext_id in range(1, len(hdul)):
                self.data[hdul[ext_id].name] = Table.read(hdul[ext_id])
            self.bands = list(self.data.keys())
    
    def __repr__(self) -> str:
        self.hdul.info().__repr__()
        return ""

    def __str__(self) -> str:
        return self.hdul.__str__()

    @property
    def concat_data(self):
        data = [Table.from_pandas(self.data[band].to_pandas()) for band in self.bands]
        data = vstack(data)
        return data

class GrizliFull:
    """
    This class is used to read the  grizli .full.fits file (line map).
    """
    def __init__(self, path):
        with fits.open(path) as hdul:
            print(hdul.info())
            print(hdul[0].header)

class GrizliRow:
    """
    This class is used to read the  grizli .row.fits file (catalog for individual source).
    """
    def __init__(self, path):
        self.data = Table.read(path)

    def __repr__(self) -> str:
        return self.data.__repr__()
    
    def __str__(self) -> str:
        return self.data.__str__()


import numpy as np



class GrizliStack:
    """
    This class is used to read the  grizli .stack.fits file (two dimensional spectrum).
    """
    def __init__(self, path):
        with fits.open(path) as hdul:
            # print(hdul.info())
            self.hdul = hdul
            self.primary_header = hdul[0].header
            self.id = self.primary_header['ID']
            self.ra = self.primary_header['RA']
            self.dec = self.primary_header['DEC']

            # bands and PAs
            header_keys = list(self.primary_header.keys())
            self.bands = [self.primary_header[key] for key in header_keys if key.startswith('GRISM')]
            self.band_info = {}
            for band in self.bands:
                self.band_info[band] = {'N_PA':self.primary_header[f"N{band}"], 'PA': [] }
                for i in range(self.primary_header[f"N{band}"]):
                    self.band_info[band]['PA'].append(self.primary_header[f"{band}{i+1:02d}"])

            # load composite data
            self.data = {}
            for band in self.bands:
                self.data[band] = self._get_data(band)
            self.header = {}
            for band in self.bands:
                self.header[band] = self._get_header(band)
        
            # load data of different PAs
            self.data_PA = {}
            for band in self.bands:
                for pa in self.band_info[band]['PA']:
                    identifier = f"{band},{pa}"
                    self.data_PA[identifier] = self._get_data(band, pa)
            self.header_PA = {}
            for band in self.bands:
                for pa in self.band_info[band]['PA']:
                    identifier = f"{band},{pa}"
                    self.header_PA[identifier] = self._get_header(band, pa)


    def _get_data(self, band, pa=None):
        """
        Return the data of a given band and pa.
        """
        data = {}
        if pa:
            identifier = f"{band},{pa}"
        else:
            identifier = f"{band}"
        for ext in ['SCI', 'WHT', 'MODEL', 'KERNEL']:
            hdu = self.hdul[ext, identifier]
            data[ext] = hdu.data
        return data

    def _get_header(self, band, pa=None):
        """
        Return the data of a given band and pa.
        """
        header = {}
        if pa:
            identifier = f"{band},{pa}"
        else:
            identifier = f"{band}"
        for ext in ['SCI', 'WHT', 'MODEL', 'KERNEL']:
            hdu = self.hdul[ext, identifier]
            header[ext] = hdu.header
        return header
    

    def spec2d(self, band, bsub=False):
        sci = self.data[band]['SCI']
        kern = self.data[band]['KERNEL']
        model = self.data[band]['MODEL']

        center_line = kern.shape[0]/2
        bottom = -kern.shape[0]/2*NIRISS_PIXEL_SCALE
        top = kern.shape[0]/2*NIRISS_PIXEL_SCALE
        sci_header = self.header[band]['SCI']

        wave = np.linspace(sci_header['WMIN'], sci_header['WMAX'], sci.shape[1])
        space = np.linspace(bottom, top, kern.shape[0])

        if(bsub):
            spec = sci - model
        else:
            spec = sci
        return wave, space, spec


    def spec2d_concat(self, bsub=False):
        # todo: currenly only work for only 2 bands, e.g., F115W and F200W
        """
        Concatenate the 2D spectrum of different bands.
        """
        waves = []
        specs = []
        for band in self.bands:
            w, s, sp = self.spec2d(band, bsub)
            edge = SPEC_EDGE_WAVELENGTH[band]
            mask = (w>=edge[0])&(w<edge[1])
            print(mask)
            waves.append(w[mask])
            specs.append(sp[:, mask])
    
        wave = np.concatenate(waves)
        space = s
        spec = np.concatenate(specs, axis=1)
            
        return wave, space, spec