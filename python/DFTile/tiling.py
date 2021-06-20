#!/usr/bin/env python
"""
Create fixed-size tiles for the survey 
"""

__author__ = "Shany Danieli"

import os
from datetime import date

import pylab as plt
import numpy as np
import pandas as pd

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord

from reproject import reproject_interp, reproject_exact
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs

datadir = '/Users/shanydanieli/projects/dwfs/data/'
# tilesdir = '/Users/shanydanieli/projects/dwfs/data/tiles/'
tilesdir = '/Users/shanydanieli/projects/dwfs/data/tiles_test/'

def spot_images(current_tile_center):
    x = current_tile_center[0]
    y = current_tile_center[1]

    fields_boundaries = np.genfromtxt(datadir+'fields_boundaries.txt', delimiter='\t', 
        dtype='S12,f8,f8,f8,f8,f8,f8,f8,f8', names=('field name', 'lower left x', 'lower left y', 
        'upper left x', 'upper left y', 'lower right x','lower right y', 'upper right x', 'upper right y'))

    relevant_fields = []
    for i in range(len(fields_boundaries)):
        lowerleft_x = fields_boundaries[i][1]
        lowerright_x = fields_boundaries[i][5]
        lowerleft_y = fields_boundaries[i][2]
        upperleft_y = fields_boundaries[i][4]

        if (lowerleft_x > x > lowerright_x) and (lowerleft_y < y < upperleft_y):
            relevant_fields.append(fields_boundaries[i][0])

    return(relevant_fields)

def make_cutout(current_tile_center,relevant_fields, filter):
    tile_side_degs = 0.7  # in degrees
    tile_side_pixels = tile_side_degs*60*60/2.5 * u.pixel

    all_cutouts = []
    all_centers_pix = []
    for i in range(len(relevant_fields)):
        field = relevant_fields[i].decode('utf-8')
        if filter == 'SloanG':
            fits_file_name = datadir+'df_coadds/'+str(field)+'/CoaddsByFilter/coadd_SloanG_'+str(field)+'.fits'
        elif filter =='SloanR':
            fits_file_name = datadir+'df_coadds/'+str(field)+'/CoaddsByFilter/coadd_SloanR_'+str(field)+'.fits'                
        'load the image'        
        hdu = fits.open(fits_file_name)[0]
        data = hdu.data
        header = hdu.header
        wcs = WCS(hdu.header)

        'convert tile center to pixel coordinates'
        w = WCS(fits_file_name)
        tile_center_pix = w.wcs_world2pix(current_tile_center[0] , current_tile_center[1] ,1)
        position = (tile_center_pix[0], tile_center_pix[1])
        all_centers_pix.append(position)
        size = tile_side_pixels
        cutout = Cutout2D(data, position, size, wcs=wcs)

        cutout_mbg = subtract_bg(cutout, header)

        # all_cutouts.append(cutout)
        all_cutouts.append(cutout_mbg)

    return(all_cutouts, all_centers_pix)

def make_weightmap_cutout(current_tile_center,relevant_fields, all_centers_pix, filter):
    tile_side_degs = 0.7  # in degrees
    tile_side_pixels = tile_side_degs*60*60/2.5 * u.pixel

    all_weightmap_cutout = []
    for i in range(len(relevant_fields)):
        field = relevant_fields[i].decode('utf-8')
        if filter == 'SloanG':
            fits_file_name = datadir+'df_coadds/'+str(field)+'/CoaddsByFilter/coadd_SloanG_'+str(field)+'_weightmap.fits'
        elif filter == 'SloanR':
            fits_file_name = datadir+'df_coadds/'+str(field)+'/CoaddsByFilter/coadd_SloanR_'+str(field)+'_weightmap.fits'
        'load the image'        
        hdu = fits.open(fits_file_name)[0]
        data = hdu.data
        wcs = WCS(hdu.header)
            
        'convert tile center to pixel coordinates'
        # w = WCS(fits_file_name)
        # w = all_cutouts[i].wcs
        # tile_center_pix = w.wcs_world2pix(current_tile_center[0] , current_tile_center[1] ,1)
        tile_center_pix = all_centers_pix[i]
        # tile_center_pix = wcs.wcs_world2pix(current_tile_center[0] , current_tile_center[1] ,1)
        position = (tile_center_pix[0], tile_center_pix[1])
        size = tile_side_pixels

        # cutout = Cutout2D(data, position, size, wcs=wcs)
        cutout = Cutout2D(data, position, size)

        all_weightmap_cutout.append(cutout)

    return(all_weightmap_cutout)


def make_hdu_list(all_cutouts):
    hdu_list = []
    for i in range(len(all_cutouts)):
        hdu = fits.PrimaryHDU()
        hdu.data = all_cutouts[i].data
        hdu.header.update(all_cutouts[i].wcs.to_header())
        hdu_list.append(hdu)

    return hdu_list


def make_weightmap_hdu_list(all_weightmap_cutout):
    hdu_list = []
    for i in range(len(all_weightmap_cutout)):
        hdu = fits.PrimaryHDU()
        hdu.data = all_weightmap_cutout[i].data
        # hdu.header.update(all_weightmap_cutout[i].wcs.to_header())
        hdu_list.append(hdu)

    return hdu_list


def save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter):
    hdu = fits.PrimaryHDU()
    hdu.data = cutout_data
    hdu.header.update(cutout_wcs.to_header())

    hdu.header = update_header(hdu.header, relevant_fields, filter)


    RA = int(current_tile_center[0]*10)
    Dec = int(abs(current_tile_center[1]*10))
    sign = np.sign(current_tile_center[1])
    if sign>0:
        coor_sign = 'p'
    elif sign<0:
        coor_sign = 'm'
    elif int(sign) == 0:
        coor_sign = 'p'
    if filter=='SloanG':
        filter_name='g'
    elif filter=='SloanR':
        filter_name='r'

    cutout_filename = 'DF'+str(RA)+str(coor_sign)+str(Dec)+str(filter_name)+'.fits'
    # cutout_filename = 'DF'+str(current_tile_center[0])+'-'+str(current_tile_center[1])+'-'+str(filter)+'.fits'
    hdu.writeto(tilesdir+cutout_filename, overwrite=True)

def reproject(hdu_list, input_weights, tile_center):

    wcs_out, shape_out = find_optimal_celestial_wcs(hdu_list, reference=tile_center)
    
    if len(hdu_list) > 1:
        array, footprint = reproject_and_coadd(hdu_list, output_projection=wcs_out,
                                       shape_out=shape_out, input_weights=input_weights, 
                                   reproject_function=reproject_exact, match_background=True)
    elif len(hdu_list) <2:
        array, footprint = reproject_and_coadd(hdu_list, output_projection=wcs_out,
                                       shape_out=shape_out, input_weights=None, 
                                   reproject_function=reproject_exact, match_background=False)

    return(array, wcs_out)

def update_header(header, relevant_fields, filter):

    today = date.today()
    header['DATE'] = (str(today), 'Date this tile was generated')
    
    used_fields = [x.decode('utf-8') for x in relevant_fields]
    used = ','.join(used_fields)
    header['USED'] = (used, 'Survey images used to generate the tile')

    header['FILTER'] = (filter)

    return header

def ready_for_tiling(field):
    fields_for_tiling = []
    with open(datadir+'ready_for_tiles.txt', 'r') as fileobj:
        for row in fileobj:
            fields_for_tiling.append(row.rstrip('\n'))
    fields_for_tiling = fields_for_tiling[1:]

    if field in fields_for_tiling:
        return True
    else:
        return False


def make_tile(tile_center, filter):
    current_tile_center = tile_center
    tc = SkyCoord(tile_center[0], tile_center[1], unit="deg")

    'find all images that overlap with this tile'    
    relevant_fields = spot_images(current_tile_center)

    'check if field is ready for tiling'
    count_ready_fields = 0
    for i in range(len(relevant_fields)):
        field = relevant_fields[i].decode("utf-8")
        ready = ready_for_tiling(field)
        if ready == True:
            count_ready_fields = count_ready_fields + 1

    'if there is at least one relevant field, make a tile'
    if count_ready_fields > 0:
        print('making a new tile, center at: '+str(tile_center[0])+', '+str(tile_center[1]) +', '+str(filter)+' band')
        'get cutouts'
        all_cutouts, all_centers_pix = make_cutout(current_tile_center, relevant_fields, filter)
        all_weightmap_cutout = make_weightmap_cutout(current_tile_center, relevant_fields, all_centers_pix, filter)

        hdu_list = make_hdu_list(all_cutouts)
        input_weights = make_weightmap_hdu_list(all_weightmap_cutout)

        'reproject and combine + save into a fits file'
        cutout_data, cutout_wcs = reproject(hdu_list, input_weights, tc)
        save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter)
        
        # 'reproject and combine + save into a fits file'
        # if len(all_cutouts)<2:
        #     cutout = all_cutouts[0]
        #     save_tile_fits(cutout.data, cutout.wcs, current_tile_center, relevant_fields, filter)
        #     # cutout_data, cutout_wcs = reproject(hdu_list, input_weights, tc)
        #     # save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter)
        # else:
        #     cutout_data, cutout_wcs = reproject(hdu_list, input_weights, tc)
        #     save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter)


def subtract_bg(cutout, header):
    bkgval = float(header['BACKVAL'])
    image = cutout.data
    image -= bkgval
    cutout.data = image

    return cutout


if __name__ == "__main__":

    # read in all tile centers for the GAMA fields
    tile_centers = np.genfromtxt(datadir+'dwfs_brick_centers_gama.txt', comments='#')

    # write tiles into a dataframe
    tiles_cat = pd.DataFrame(tile_centers)
    tiles_cat.columns = ['ra','dec']


    # # read in the names of fields that are ready for tiling
    # fields_for_tiling = []
    # with open(datadir+'ready_for_tiles.txt', 'r') as fileobj:
    #     for row in fileobj:
    #         fields_for_tiling.append(row.rstrip('\n'))
    # fields_for_tiling = fields_for_tiling[1:]


    # loop over all potential tile centers 
    # only make tiles for fields that are indicated as ready 
    # in the ready_for_tiles.txt file
    print('Loop over all potential tiles')
    print('-------------------------------')
    # for i in range(len(tiles_cat)):
    for i in range(1):
        tile_center = [tiles_cat.iloc[i].ra,tiles_cat.iloc[i].dec]
        make_tile(tile_center, filter='SloanG')
        make_tile(tile_center, filter='SloanR')
        
        


