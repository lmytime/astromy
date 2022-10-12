import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.stats import sigma_clip
from astropy.convolution import convolve,Gaussian2DKernel
from spectral_cube import SpectralCube
from scipy import ndimage
import regions
from astropy import units as u
from astropy.coordinates import SkyCoord

class AstroCube(SpectralCube):

    def __init__(self, cube, **krg):
        '''
        初始化cube实例
        '''
        # super().__init__(**krg) #因为继承自SpectralCube,这句是初始化父类
        self.cube = cube

    @classmethod#类方法
    def read(cls,filename,**krg):
        '''
        class function, read the datacube from .fits
        '''
        cube = super().read(filename=filename, **krg)
        # cube = SpectralCube(data=cube.hdu.data, wcs=cube.wcs, **krg)
        return AstroCube(cube)

    def aper_spec(self, coord, radius):
        '''
        从cube中提取一个圆形区域的光谱
        coord:中心坐标
        radius:半径
        '''
        regpix = regions.CircleSkyRegion(center=SkyCoord(coord[0], coord[1], unit='deg'), radius=radius*u.arcsec)
        subcube = self.cube.subcube_from_regions([regpix])
        spectrum = subcube.sum(axis=(1,2))
        return spectrum