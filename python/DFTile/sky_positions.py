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

def function():
    pass


if __name__ == "__main__":

    field_list = os.listdir(datadir)

    '''
    G15_222 - has two consecutive fields in the declination 
    '''

    G15_222_m1_G_fits = datadir+'raw/G15_222_m1/CoaddsByFilter/coadd_SloanG_G15_222_m1.fits'
    G15_222_1_G_fits = datadir+'raw/G15_222_1/CoaddsByFilter/coadd_SloanG_G15_222_1.fits'

    hdu_G15_222_m1 = fits.open(G15_222_m1_G_fits)
    data_G15_222_m1 = hdu_G15_222_m1[0].data

    hdu_G15_222_1 = fits.open(G15_222_1_G_fits)
    data_G15_222_1 = hdu_G15_222_1[0].data


    lowerleft_pix = [0, 0]
    upperleft_pix = [0, data_G15_222_1.shape[0]]
    lowerright_pix = [data_G15_222_m1.shape[1], 0]
    upperright_pix = [data_G15_222_1.shape[1], data_G15_222_1.shape[0]]

    w_G15_222_m1 = wcs.WCS(G15_222_m1_G_fits)
    w_G15_222_1 = wcs.WCS(G15_222_1_G_fits)

    lowerleft_coor = w_G15_222_m1.wcs_pix2world(lowerleft_pix[0],lowerleft_pix[1], 1)
    upperleft_coor = w_G15_222_1.wcs_pix2world(upperleft_pix[0],upperleft_pix[1], 1)
    lowerright_coor = w_G15_222_m1.wcs_pix2world(lowerright_pix[0],lowerright_pix[1], 1)
    upperright_coor = w_G15_222_1.wcs_pix2world(upperright_pix[0],upperright_pix[1], 1)


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
    f1 = open('../data/G222_bricks.reg', 'w')
    f1.write('# Region file format: DS9 version 4.1\n')
    f1.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f1.write('fk5\n')
    for i in range(len(brick_centers)):
        f1.write('box('+str(brick_centers[i][0])+','+str(brick_centers[i][1])+',2520.000",2520.000",0)\n')


    'write into a file'
    f2 = open('../data/dwfs_brick_centers.txt', 'w')
    for i in range(len(brick_centers)):
        f2.write(str(round(brick_centers[i][0],3))+'\t'+str(round(brick_centers[i][1],3))+'\n')


