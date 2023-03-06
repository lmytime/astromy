# Develop for twodspec

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astromy import AstroImage, combine_RGB, zscale, gamma_correction
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from reproject import reproject_adaptive, reproject_exact, reproject_interp
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
from scipy.signal import find_peaks
from spectres import spectres

class Astro1dSpec:
    # 1D spectrum class
    def __init__(self, data, wave_axis=None, err=None):
        self.data = data
        self.wave_axis = wave_axis
        self.err = err
    
    def rebin(self, new_wave, fill=None, verbose=True):
        # rebin spectrum to a new wavelength axis
        new_fluxes = spectres(new_wavs=new_wave, spec_wavs=self.wave_axis, spec_fluxes=self.data, spec_errs=self.err, fill=fill, verbose=verbose)
        if(self.err is not None):
            new_fluxes, new_errs = new_fluxes
        else:
            new_errs = None
        return Astro1dSpec(data=new_fluxes, wave_axis=new_wave, err=new_errs)

    def preview(self, gamma=1.0, colorbar=True, **kwargs):
        fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
        ax.plot(self.wave_axis, self.data, drawstyle='steps-mid', color='k')
        if(self.err is not None):
            ax.fill_between(self.wave_axis, self.data-self.err,self.data+self.err, alpha=0.2, color='gray', step='mid')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        return fig, ax

    def __repr__(self, plot=True):
        __info__ = """
        Spectral information:
        ------------------
        Waevelength range: {:.3f} - {:.3f}
        Pixel scale: {:.3f} /pixel
        """.format(self.wave_axis[0], self.wave_axis[-1], self.wave_axis[1]-self.wave_axis[0])
        if(plot):
            self.preview()
        return __info__

    def save(self, url, format='fits', overwrite=True):
        # save file to csv
        if(format=='csv'):
            pd.DataFrame({'wave': self.wave_axis, 'flux': self.data, 'err':self.err}).to_csv(url, index=False)
        # save file to fits
        elif(format=='fits'):
            col1 = fits.Column(name='wavelength', format='D', array=self.wave_axis)
            col2 = fits.Column(name='flux', format='D', array=self.data)
            col3 = fits.Column(name='error', format='D', array=self.err)
            hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
            hdu.writeto(url, overwrite=overwrite)
