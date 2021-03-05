#!/usr/bin/env python
"""
Create fixed-size tiles for the survey 
"""

__author__ = "Shany Danieli"

import os
from datetime import date

import pylab as plt
import numpy as np

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D

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
    for i in range(len(relevant_fields)):
        field = relevant_fields[i].decode('utf-8')
        if filter == 'SloanG':
            fits_file_name = datadir+'raw/'+str(field)+'/CoaddsByFilter/coadd_SloanG_'+str(field)+'.fits'
        elif filter =='SloanR':
            fits_file_name = datadir+'raw/'+str(field)+'/CoaddsByFilter/coadd_SloanR_'+str(field)+'.fits'                
        'load the image'        
        hdu = fits.open(fits_file_name)[0]
        data = hdu.data
        wcs = WCS(hdu.header)

        'convert tile center to pixel coordinates'
        w = WCS(fits_file_name)
        tile_center_pix = w.wcs_world2pix(current_tile_center[0] , current_tile_center[1] ,1)
        position = (tile_center_pix[0], tile_center_pix[1])
        size = tile_side_pixels
        cutout = Cutout2D(data, position, size, wcs=wcs)
        all_cutouts.append(cutout)

    return(all_cutouts)

def make_weightmap_cutout(current_tile_center,relevant_fields, filter):
    tile_side_degs = 0.7  # in degrees
    tile_side_pixels = tile_side_degs*60*60/2.5 * u.pixel

    all_weightmap_cutout = []
    for i in range(len(relevant_fields)):
        field = relevant_fields[i].decode('utf-8')
        if filter == 'SloanG':
            fits_file_name = datadir+'raw/'+str(field)+'/CoaddsByFilter/coadd_SloanG_'+str(field)+'_weightmap.fits'
        elif filter == 'SloanR':
            fits_file_name = datadir+'raw/'+str(field)+'/CoaddsByFilter/coadd_SloanR_'+str(field)+'_weightmap.fits'
        'load the image'        
        hdu = fits.open(fits_file_name)[0]
        data = hdu.data
        wcs = WCS(hdu.header)
            
        'convert tile center to pixel coordinates'
        w = WCS(fits_file_name)
        tile_center_pix = w.wcs_world2pix(current_tile_center[0] , current_tile_center[1] ,1)
        position = (tile_center_pix[0], tile_center_pix[1])
        size = tile_side_pixels
        cutout = Cutout2D(data, position, size, wcs=wcs)
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
        hdu.data = all_cutouts[i].data
        hdu.header.update(all_cutouts[i].wcs.to_header())
        hdu_list.append(hdu)

    return hdu_list


def save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter):
    hdu = fits.PrimaryHDU()
    hdu.data = cutout_data
    hdu.header.update(cutout_wcs.to_header())

    hdu.header = update_header(hdu.header, relevant_fields, filter)

    cutout_filename = 'tile_'+str(current_tile_center[0])+'_'+str(current_tile_center[1])+'_'+str(filter)+'.fits'
    hdu.writeto(tilesdir+cutout_filename, overwrite=True)

def reproject(hdu_list, input_weights):

    wcs_out, shape_out = find_optimal_celestial_wcs(hdu_list)
    array, footprint = reproject_and_coadd(hdu_list, output_projection=wcs_out,
                                       shape_out=shape_out, input_weights=input_weights, 
                                       reproject_function=reproject_exact, match_background=True)
    return(array, wcs_out)

def update_header(header, relevant_fields, filter):

    today = date.today()
    header['DATE'] = (str(today), 'Date this tile was generated')
    
    used_fields = [x.decode('utf-8') for x in relevant_fields]
    used = ','.join(used_fields)
    header['USED'] = (used, 'Survey images used to generate the tile')

    header['FILTER'] = (filter)


    return header


if __name__ == "__main__":
                
        tile_centers = np.loadtxt(datadir+'dwfs_brick_centers.txt')

        filter = 'SloanG'
        for i in range(len(tile_centers[:5])):
            print('.........Generating a new tile.........')
            print('Tile # '+str(i))
            print('Tile center is at: ' + str(tile_centers[i]))
            print('Filter: ' + str(filter))
            current_tile_center = tile_centers[i]
            'find all images that overlap with this tile'
            relevant_fields = spot_images(current_tile_center)
            'get their cutouts'
            all_cutouts = make_cutout(current_tile_center, relevant_fields, filter)
            all_weightmap_cutout = make_weightmap_cutout(current_tile_center, relevant_fields, filter)

            hdu_list = make_hdu_list(all_cutouts)
            input_weights = make_weightmap_hdu_list(all_weightmap_cutout)

            'reproject and combine + save into a fits file'
            if len(all_cutouts)<2:
                cutout = all_cutouts[0]
                save_tile_fits(cutout.data, cutout.wcs, current_tile_center, relevant_fields, filter)
            else:
                cutout_data, cutout_wcs = reproject(hdu_list, input_weights)
                save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter)


        # filter = 'SloanR'
        # for i in range(len(tile_centers[:6])):
        #     print('.........Generating a new tile.........')
        #     print('Tile # '+str(i))
        #     print('Tile center is at: ' + str(tile_centers[i]))
        #     print('Filter: ' + str(filter))
        #     current_tile_center = tile_centers[i]
        #     'find all images that overlap with this tile'
        #     relevant_fields = spot_images(current_tile_center)
        #     'get their cutouts'
        #     all_cutouts = make_cutout(current_tile_center, relevant_fields, filter)
        #     all_weightmap_cutout = make_weightmap_cutout(current_tile_center, relevant_fields, filter)

        #     hdu_list = make_hdu_list(all_cutouts)
        #     input_weights = make_weightmap_hdu_list(all_weightmap_cutout)

        #     'reproject and combine + save into a fits file'
        #     if len(all_cutouts)<2:
        #         cutout = all_cutouts[0]

        #         # print(cutout.wcs)

        #         save_tile_fits(cutout.data, cutout.wcs, current_tile_center, relevant_fields, filter)
        #     else:
        #         cutout_data, cutout_wcs = reproject(hdu_list, input_weights)                
        #         save_tile_fits(cutout_data, cutout_wcs, current_tile_center, relevant_fields, filter)



