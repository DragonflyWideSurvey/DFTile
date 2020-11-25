#!/usr/bin/env python
"""
Derive the boundaries of all the dwfs fields and write them into a file 
"""

__author__ = "Shany Danieli"

import os
import pylab as plt
import numpy as np
import astropy.wcs as wcs
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs



datadir = '/Users/shanydanieli/projects/dwfs/data/raw/'


if __name__ == "__main__":

    field_list = os.listdir(datadir)

    lowerleft_coor = []
    upperleft_coor = []
    lowerright_coor = []
    upperright_coor = []

    for i in range(len(field_list)):
        current_field = field_list[i]

        file_name = datadir+str(current_field)+'/CoaddsByFilter/coadd_SloanG_'+str(current_field)+'.fits'
        hdu = fits.open(file_name)
        data = hdu[0].data

        lowerleft_pix = [0, 0]
        upperleft_pix = [0, data.shape[0]]
        lowerright_pix = [data.shape[1], 0]
        upperright_pix = [data.shape[1], data.shape[0]]

        w = wcs.WCS(file_name)

        lowerleft_coor.append(w.wcs_pix2world(lowerleft_pix[0],lowerleft_pix[1], 1))
        upperleft_coor.append(w.wcs_pix2world(upperleft_pix[0],upperleft_pix[1], 1))
        lowerright_coor.append(w.wcs_pix2world(lowerright_pix[0],lowerright_pix[1], 1))
        upperright_coor.append(w.wcs_pix2world(upperright_pix[0],upperright_pix[1], 1))

    'write fields boundary coordinates into a file'
    f = open('../data/fields_boundaries.txt', 'w')
    f.write('# field name\n')
    f.write('# lower left x\n')
    f.write('# lower left y\n')
    f.write('# upper left x\n')
    f.write('# upper left y\n')
    f.write('# lower right x\n')
    f.write('# lower right y\n')
    f.write('# upper right x\n')
    f.write('# upper right y\n')
    for i in range(len(lowerleft_coor)): 
        f.write(str(field_list[i])+'\t'
            +str(lowerleft_coor[i][0])+'\t'+str(lowerleft_coor[i][1])+'\t'
            +str(upperleft_coor[i][0])+'\t'+str(upperleft_coor[i][1])+'\t'
            +str(lowerright_coor[i][0])+'\t'+str(lowerright_coor[i][1])+'\t'
            +str(upperright_coor[i][0])+'\t'+str(upperright_coor[i][1])+'\n')

