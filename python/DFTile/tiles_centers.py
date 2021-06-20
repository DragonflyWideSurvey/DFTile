#!/usr/bin/env python
"""
Create the centers of tiles for the survey 
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

datadir = '/Users/shanydanieli/projects/dwfs/data/'


if __name__ == "__main__":

    # lower left corner: G09_130.5_m1
    # upper left corner G09_130.5_1
    # lower right corner: G15_222_m1
    # upper right corner: G15_222_1

    # get the lower left corner pixel
    image1 = datadir+'df_coadds/G09_130.5_m1/CoaddsByFilter/coadd_SloanG_G09_130.5_m1.fits'
    hdu1 = fits.open(image1)
    data1 = hdu1[0].data
    wcs1 = wcs.WCS(image1)
    # lowerleft_pix = [0, 0]
    lowerright_pix = [data1.shape[1], 0]
    lowerright_coor = wcs1.wcs_pix2world(lowerright_pix[0],lowerright_pix[1], 1)

    # get the upper left corner pixel
    image2 = datadir+'df_coadds/G09_130.5_1/CoaddsByFilter/coadd_SloanG_G09_130.5_1.fits'
    hdu2 = fits.open(image2)
    data2 = hdu2[0].data
    wcs2 = wcs.WCS(image2)
    # upperleft_pix = [0, data2.shape[0]]
    upperright_pix = [data2.shape[1], data2.shape[0]]
    upperright_coor = wcs2.wcs_pix2world(upperright_pix[0],upperright_pix[1], 1)

    # get the lower right corner pixel
    image3 = datadir+'df_coadds/G15_222_m1/CoaddsByFilter/coadd_SloanG_G15_222_m1.fits'
    hdu3 = fits.open(image3)
    data3 = hdu3[0].data
    wcs3 = wcs.WCS(image3)
    # lowerright_pix = [data3.shape[1], 0]
    lowerleft_pix = [0, 0]
    lowerleft_coor = wcs3.wcs_pix2world(lowerleft_pix[0],lowerleft_pix[1], 1)

    # get the upper right corner pixel
    image4 = datadir+'df_coadds/G15_222_1/CoaddsByFilter/coadd_SloanG_G15_222_1.fits'
    hdu4 = fits.open(image4)
    data4 = hdu4[0].data
    wcs4 = wcs.WCS(image4)
    # upperright_pix = [data4.shape[1], data4.shape[0]]
    upperleft_pix = [0, data4.shape[0]]
    upperleft_coor = wcs4.wcs_pix2world(upperleft_pix[0],upperleft_pix[1], 1)

    # GAMA
    lowerleft_coor = np.array([223.9,-2.45])
    upperleft_coor = np.array([223.9,2.45])
    lowerright_coor = np.array([128.5,-2.45])
    upperright_coor = np.array([128.5,2.45])

    current_center_x = lowerleft_coor[0] - 0.35
    current_center_y = lowerleft_coor[1] + 0.35

    brick_centers = []

    while (current_center_x - 0.1 > upperright_coor[0]):
        while  (current_center_y + 0.1 < upperright_coor[1]):

            brick_centers.append([current_center_x,current_center_y])
            current_center_y = current_center_y + 0.7

        current_center_x = current_center_x - 0.7
        current_center_y = lowerleft_coor[1] + 0.35
    

    'write into a region file'
    f1 = open(datadir+'regions/all_bricks_gama.reg', 'w')
    f1.write('# Region file format: DS9 version 4.1\n')
    f1.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f1.write('fk5\n')
    for i in range(len(brick_centers)):
        f1.write('box('+str(brick_centers[i][0])+','+str(brick_centers[i][1])+',2520.000",2520.000",0)\n')


    'write into a file'
    # f2 = open('../data/dwfs_brick_centers.txt', 'w')
    f2 = open(datadir+'dwfs_brick_centers_gama.txt', 'w')
    for i in range(len(brick_centers)):
        f2.write(str(round(brick_centers[i][0],3))+'\t'+str(round(brick_centers[i][1],3))+'\n')

