#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:18:38 2018

@author: andrew
"""

import numpy as np
from astropy.convolution import convolve
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image
from psf import psf_construct
import sep
from astropy.io import fits
import sys

def make_sources(x_size, y_size, num_sources, psf=[], flux=50000):
    image_size = (x_size, y_size)
    
    num = num_sources
    fluxes = flux*np.random.random(size=num)
    
    image = np.zeros(image_size)
    
    for i in range(num):
        x_loc, y_loc = np.random.randint(0,x_size), np.random.randint(0,y_size)
        image[x_loc][y_loc] = fluxes[i]
        
        
    for p in psf:
        image = convolve(image, p)
    
    return image, x_loc, y_loc, fluxes

def make_image(size_x, size_y, x_loc=[], y_loc=[], fluxes=[], psf=[], 
               real=False, real_template=False, im='', add_bkg=False, 
               bkg_loc=70, bkg_rms=15, bkg_mode='sep'):
    image = np.zeros((size_x, size_y))
    num_sources = len(fluxes)
    
    if real == True:
        for source in range(num_sources):
            test_image = np.zeros((size_x, size_y))
            test_image[int(x_loc[source])-1][int(y_loc[source])-1] = fluxes[source]
            test_image = convolve(test_image, psf_construct(im, x_loc[source]-1, y_loc[source]-1, template=real_template))
            image += test_image
    else:
        for source in range(num_sources):
            image[x_loc[source]-1][y_loc[source]-1] = fluxes[source]
        for p in psf:
            image = convolve(image, p)
    
    if add_bkg == True:
        if bkg_mode == 'sep':
            im_data = fits.getdata(im)
            im_mask = fits.getdata(im, 1)
            try: weight_check = fits.getval(im, 'WEIGHT')
            except: weight_check = 'N'
            if weight_check == 'Y':
                im_mask = (im_mask-1)*-1
            try: bkg = sep.Background(im_data, mask=im_mask)
            except:
                im_mask = im_mask.astype(bool)
                im_data_sep = im_data.byteswap().newbyteorder()
                bkg = sep.Background(im_data_sep, mask=im_mask)
            back = bkg.back()
            rms = bkg.rms()
            gauss_bkg = np.random.normal(loc=back, scale=rms)
            image += gauss_bkg
        elif bkg_mode == 'gauss':
            gauss_bkg = np.random.normal(loc=bkg_loc, scale=bkg_rms, size=(size_x, size_y))
            image += gauss_bkg
        else:
            print("-> Error: 'bkg_mode' keyword must be either 'sep' or 'gauss'\n-> Exiting...")
            sys.exit()
    
    return image

def get_moffat_gamma(fwhm, alpha=4.765):
    gamma = fwhm/(2*np.sqrt(2**(1/alpha)-1))
    return gamma

def make_gaussian_im(x_size, y_size, fluxes=[100,1000,10000], x_pos=[500,250,750], y_pos=[300,80,460],std=[6,6,6]):
    shape = (x_size, y_size)
    table = Table()
    table['flux'] = fluxes
    table['x_mean'] = x_pos
    table['y_mean'] = y_pos
    table['x_stddev'] = std
    table['y_stddev'] = std
    image = make_gaussian_sources_image(shape, table)
    return image