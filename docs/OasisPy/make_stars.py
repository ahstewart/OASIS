#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:18:38 2018

@author: andrew
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.stats import biweight_location
from astropy.stats import mad_std
import sep
from photutils.datasets import make_gaussian_sources_image
from photutils.psf import create_matching_kernel
from astropy.table import Table
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel, Moffat2DKernel

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

def make_image(size_x, size_y, x_loc=[], y_loc=[], fluxes=[], psf=[]):
    image = np.zeros((size_x, size_y))
    num_sources = len(fluxes)
    for source in range(num_sources):
        image[x_loc[source]-1][y_loc[source]-1] = fluxes[source]
    
    for p in psf:
        image = convolve(image, p)
    
    return image

def add_noise(image, x_size=1000, y_size=1000, norm_loc=90, norm_scale=20):
    image = (np.random.poisson(image)).astype(np.float64)
    noise = np.random.normal(np.random.normal(loc=norm_loc,scale=norm_scale,size=(x_size,y_size)))
    image += noise
    return image

def get_moffat_gamma(fwhm, alpha=4.765):
    gamma = fwhm/(2*np.sqrt(2**(1/alpha)-1))
    return gamma
    
def plot_image(array):
    mean, std = np.mean(array), np.std(array)
    plt.imshow(array, cmap='Greys_r', origin='lower', vmin = mean-std, vmax=mean+std)
    
def make_constant_bkg(x_size, y_size, sky_brightness, sig):
    bkg = np.random.normal(size=(x_size,y_size), loc=sky_brightness, scale=sig)
    return bkg

def make_variable_bkg(x_size, y_size, sky_brightness, iterator, sig):
    bkg = np.random.random(size=(x_size,y_size))*sky_brightness
    for i in range(iterator):
        bkg = gaussian_filter(bkg, sigma=sig)
    return bkg

def sig_clip(array, sig=3.0, iterators=5):
    mean,median,std = sigma_clipped_stats(array, sigma=sig, iters=iterators)
    return mean,median,std
    
def get_back(image):
    med_back = np.median(image)
    biweight_back = biweight_location(image)
    rms = mad_std(image)
    sep_back = sep.Background(image)
    sig_back_mean, sig_back_median, sig_back_rms = sig_clip(image)
    print('\nNumpy Median: %f, %f\nBiweight: %f, %f\nSEP: %f, %f\nSigma clip: %f, %f' 
          % (med_back, rms, biweight_back, rms, sep_back.globalback, sep_back.globalrms, 
             sig_back_median, sig_back_rms))
    
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

def make_other_im(sig):
    h, w = 2000,3000
    pos_x = [1500,1000,1200]
    pos_y = [1600,1400,200]
    array = np.array([ 0.65343465,  0.50675629,  0.84946314])
    fluxes = 2000000.0 + array * 300.0
    img = np.zeros((h,w))
    for x, y, f in zip(pos_x, pos_y, fluxes):
        img[x, y] = f
    
    img = gaussian_filter(img, sigma=sig, mode='constant')
    
    final = fits.PrimaryHDU(img)
    final.writeto('/home/andrew/scripts/g_im2.fits', overwrite=True)
    
    
    

