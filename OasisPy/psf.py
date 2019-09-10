#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 20:09:18 2018

@author: andrew
"""

import glob
import os
import initialize
from astropy.io import fits
import sex
from tqdm import tqdm
import sys
import numpy as np

def psfex(location, images=[]):
    initialize.create_configs(location)
    config_loc = location + '/configs/psfex.config'
    with open(config_loc, 'r') as config:
        data = config.readlines()
        config.close()
    data[83] = "PSF_DIR" + "        " + location + "/psf" + "\n"
    with open(config_loc, 'w') as config:
        config.writelines(data)
        config.close()
    print("-> Calculating PSFs...")
    for im in tqdm(images):
        cat = im.replace('data','psf')
        cat = cat.replace('templates', 'psf')
        cat = cat[:-5] + '.cat'
        os.system("psfex %s > %s.psf -c %s" % (cat, cat[:-4], config_loc))
    if images == []:
        print("-> PSFs already exist")
    
def fwhm(image):
    psf = image.replace('fits', 'psf')
    psf = psf.replace('data', 'psf')
    psf_data = fits.open(psf)
    fwhm = psf_data[1].header['PSF_FWHM']
    return fwhm  

def fwhm_template(image):
    psf = image.replace('fits', 'psf')
    psf = psf.replace('templates', 'psf')
    psf_data = fits.open(psf)
    fwhm = psf_data[1].header['PSF_FWHM']
    return fwhm  

def check_fwhm(location):
    FWHM = []
    temps = glob.glob(location + '/templates/*.fits')
    images = glob.glob(location + '/data/*.fits')
    for i in images:
        FWHM.append(fwhm(i))
    for t in temps:
        FWHM.append(fwhm_template(t))
    return FWHM

def get_first_model(image):
    psf_loc = image.replace('data', 'psf')
    psf_loc = psf_loc.replace('fits', 'psf')
    psf_hdu = fits.open(psf_loc)
    psf_bintable = psf_hdu[1].data
    psf_mask = psf_bintable['PSF_MASK']
    first_model = psf_mask[0][0]
    psf_hdu.close()
    return first_model

def get_all_models(image):
    psf_loc = image.replace('data', 'psf')
    psf_loc = psf_loc.replace('fits', 'psf')
    psf_hdu = fits.open(psf_loc)
    psf_bintable = psf_hdu[1].data
    psf_mask = psf_bintable['PSF_MASK']
    psf_models = []
    for model in psf_mask[0]:
        psf_models.append(model)
    psf_hdu.close()
    return psf_models

def psf_homo(location):
    images = glob.glob("%s/psf_2/*.cat" % (location))
    for i in tqdm(images):
        os.system("psfex %s -c %s/configs/psfex_homo.config" % (i, location))

def psf_construct(image, x, y, normalize=True, verbose=False, template=False):
    if template == False:
        psf = image.replace('/data/', '/psf/')
        psf = psf.replace('.fits', '.psf')
    elif template == True:
        psf = image.replace('/templates/', '/psf/')
        psf = psf.replace('.fits', '.psf')   
    else:
        print("-> Error: 'template' keyword must be boolean\n-> Exiting...")
        sys.exit()
    try: 
        psf_model = (fits.getdata(psf))['PSF_MASK'][0]
        psf_header = fits.getheader(psf, 1)
    except:
        print("-> Error: Can't find PSF file for %s\n-> Exiting..." % (image))
        sys.exit()
    psf_deg = psf_header['POLDEG1']
    x_offset = psf_header['POLZERO1']
    x_scale = psf_header['POLSCAL1']
    X = (x-x_offset)/x_scale
    y_offset = psf_header['POLZERO2']
    y_scale = psf_header['POLSCAL2']
    Y = (y-y_offset)/y_scale
    if verbose == True:
        print("\nPSF degree: %d\nScaled X: %.2f\nScaled Y: %.2f\n" % (psf_deg, X, Y))
    if psf_deg == 2:
        construct = psf_model[0] + (X*psf_model[1]) + (X*X*psf_model[2]) + (Y*psf_model[3]) + (X*Y*psf_model[4]) + (Y*Y*psf_model[5])
    elif psf_deg == 1:
        construct = psf_model[0] + (X*psf_model[1]) + (Y*psf_model[2])
    elif psf_deg == 0:
        construct = psf_model[0]
    if normalize == True:
        return construct/np.sum(construct)
    else:
        return construct
    

def PSF(path):
    '''Computes PSF model of science images. Uses ``PSFex`` (Bertin).
    
    :param str path: Path of data file tree (contains the **configs**, **data**, **psf**, **residuals**, **sources**, **templates** directories). Use a comma-separated list for mapping to multiple datasets.
    :returns: PSF models of each science image are ouputted into the **psf** directory with the *.psf* suffix.
    
    '''
    paths = (path.replace(' ','')).split(',')
    del path
    for path in paths:
        ims = sex.sextractor_psf(path)
        psfex(path, images=ims)   
    
if __name__ == '__main__':
    path = input("-> Enter path to exposure time directory: ")
    PSF(path)