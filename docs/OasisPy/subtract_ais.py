#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 17:24:29 2018

@author: andrew
"""

import glob
import os
import shutil
import initialize
from intensity_match import int_match_to_template
from background import background_match_to_template
from astropy.io import fits
import sys
import multiprocessing
import numpy as np
import psf
#import signal

def run_subtraction(location):
    processors = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=processors)
    r = pool.map(isis_sub, location)
    pool.close()

def isis_sub(location):
    x = 0
    images = glob.glob(location + "/data/*_A_.fits")
    template = glob.glob(location + "/templates/*.fits")
    residuals = glob.glob(location + "/residuals/*residual_.fits")
    images_names = [(i.split('/')[-1])[:-5] for i in images]
    res_names = [(r.split('/')[-1])[:-14] for r in residuals]
    resids = [res for res in images_names if res not in res_names]
    ims = []
    for rs in resids:
        ims.append(location+'/data/'+rs+'.fits')
    if ims != []:
        if len(template) == 1:
            ais_loc = os.path.dirname(initialize.__file__) + "/AIS/package/bin/./mrj_phot"
            initialize.create_configs(location)
            ais_config_loc = location + '/configs/default_config'
            cwd = os.getcwd()
            psf_data = glob.glob(location + '/psf/*')
            template_mask = fits.getdata(template[0], 1)
            if len(psf_data) == 2*(len(images)+1):
#                fwhm_template = psf.fwhm_template(template[0])
                try:
                    os.mkdir(cwd + "/AIS_temp")
                except FileExistsError:
                    pass
                os.chdir(cwd + "/AIS_temp")
                length = len(location) + 5
                print("\n-> Subtracting images...")
                for i in ims:
                    int_match_to_template(location, i, template[0])
#                    background_match_to_template(i, template[0], sigma=20)
#                    fwhm_image = psf.fwhm(i)
#                    if fwhm_template < fwhm_image:
#                        os.system(ais_loc + " " + template[0] + " " + i + " -c " + ais_config_loc)
#                        os.system("mv -f %s/AIS_temp/conv.fits %s/residuals/%sresidual_.fits" % (cwd, location, i[length:-5]))
#                        invert_image(location + '/residuals/' + i[length:-5] + 'residual_.fits')
#                    else:
                    os.system(ais_loc + " " + i + " " + template[0] + " -c " + ais_config_loc)
                    os.system("mv -f %s/AIS_temp/conv.fits %s/residuals/%sresidual_.fits" % (cwd, location, i[length:-5]))
                    hdu = fits.open(location + '/residuals/' + i[length:-5] + 'residual_.fits', mode='update')
                    hdr = hdu[0].header
                    hdr.set('OPTIMIZE', 'N')
                    hdu.close()
                    image_hdu = fits.open((i.replace('residual_', '')).replace('residuals', 'data'))
                    image_hduMask = np.logical_or(np.logical_not(image_hdu[1].data), np.logical_not(template_mask)).astype(int)
                    image_hdu.close()
                    hdu = fits.open(location + '/residuals/' + i[length:-5] + 'residual_.fits')
                    data = hdu[0].data
                    hdr = hdu[0].header
                    hdu.close()
                    hduData = fits.PrimaryHDU(data, header=hdr)
                    hduMask = fits.ImageHDU(image_hduMask)
                    hduList = fits.HDUList([hduData, hduMask])
                    hduList.writeto(location + '/residuals/' + i[length:-5] + 'residual_.fits', overwrite=True)
                    x += 1
                    per = float(x)/float(len(ims)) * 100
                    print("\t %.1f%% subtracted..." % (per))
    #                    else:
    #                        os.system(ais_loc + " " + i + " " + template[0] + " -c " + ais_config_loc)
    #                        os.system("mv -f %s/AIS_temp/conv.fits %s/residuals/%sresidual_.fits" % (cwd, location, i[length:-5]))
    #                        x += 1
    #                        per = float(x)/float(len(images)) * 100
    #                        print("\n-> %.1f%% subtracted..." % (per))
            else:
                print("-> Error: Need PSFs before running subtraction\n-> Run psf.py first")
                print("-> If any images have been manually removed from the data directory, delete all contents of the psf directory and run OasisPy again\n")
                sys.exit()
        else:
            print("-> Subtraction failure: Template missing")
            sys.exit()
        os.chdir(cwd)
        shutil.rmtree(cwd + "/AIS_temp")
    else:
        print("-> Images have already been subtracted")
        
def isis_sub_sim(location, image):
    template = glob.glob(location + "/templates/*.fits")
    if len(template) == 1:
        ais_loc = os.path.dirname(initialize.__file__) + "/AIS/package/bin/./mrj_phot"
        initialize.create_configs(location)
        ais_config_loc = location + '/configs/default_config'
        cwd = os.getcwd()
        try:
            os.mkdir(cwd + "/AIS_temp")
        except FileExistsError:
            pass
        os.chdir(cwd + "/AIS_temp")
        length = len(location) + 5
        print("\n-> Subtracting template from fake image...")
        int_match_to_template(location, image, template[0])
#       background_match_to_template(i, template[0])
        os.system(ais_loc + " " + template[0] + " " + image + " -c " + ais_config_loc)
        os.system("mv -f %s/AIS_temp/conv.fits %s/residuals/%sresidual_.fits" % (cwd, location, image[length:-5]))
        invert_image(location + '/residuals/' + image[length:-5] + 'residual_.fits')
        hdu = fits.open(location + '/residuals/' + image[length:-5] + 'residual_.fits', mode='update')
        hdr = hdu[0].header
        hdr.set('OPTIMIZE', 'N')
        hdu.close()
        image_hdu = fits.open((image.replace('residual_', '')).replace('residuals', 'data'))
        image_hduMask = (image_hdu[1].data).astype(int)
        image_hdu.close()
        hdu = fits.open(location + '/residuals/' + image[length:-5] + 'residual_.fits')
        data = hdu[0].data
        hdr = hdu[0].header
        MSK = image_hduMask
        hdu.close()
        hduData = fits.PrimaryHDU(data, header=hdr)
        hduMask = fits.ImageHDU(MSK)
        hduList = fits.HDUList([hduData, hduMask])
        hduList.writeto(location + '/residuals/' + image[length:-5] + 'residual_.fits', overwrite=True)
    else:
        print("\n-> Subtraction failure: Template missing")
        
def invert_residuals(location):
    residuals = glob.glob(location + '/residuals/*.fits')
#    length = len(location) + 11
    for r in residuals:
#        name = r[length:]
        r_hdu = fits.open(r)
        data = r_hdu[0].data
        data *= -1
        hdu = fits.PrimaryHDU(data, header=r_hdu[0].header)
        hdu.writeto(r, overwrite=True)
#        hdu.writeto(location+'/residuals/'+'INVERT_'+name)
#        os.system("rm %s" % (r))
    
def invert_image(image):
    r_hdu = fits.open(image)
    data = r_hdu[0].data
    data *= -1
    hdu = fits.PrimaryHDU(data, header=r_hdu[0].header)
    hdu.writeto(image, overwrite=True)