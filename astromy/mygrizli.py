from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack


class Grizli1D:
    def __init__(self, path):
        with fits.open(path) as hdul:
            self.primary_header = hdul[0].header
            self.id = self.primary_header['ID']
            self.ra = self.primary_header['RA']
            self.dec = self.primary_header['DEC']
            self.target = self.primary_header['TARGET']
            self.data = {}
            for ext_id in range(1, len(hdul)):
                self.data[hdul[ext_id].name] = Table.read(hdul[ext_id])
            self.band = self.data.keys()
            
    def concat_band(self):
        data = [self.data[band] for band in self.band]
        
        data = vstack(data)
        print(data)
        # data = data.sort('wave')
        # return data
