#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 20:12:31 2018

@author: andrew
"""

import glob
from image_registration import chi2_shift
from scipy.ndimage import shift
from astropy.io import fits
import numpy as np
import os

def chi2(location):
    print("-> Aligning images with the chi2 method...")
    x = 0
    ref = glob.glob(location + '/*_ref_A_.fits')
    images = glob.glob(location + '/*_N_.fits')
    if len(ref) == 1:
        ref_data = fits.getdata(ref[0])
        ref_data = np.array(ref_data, dtype='float64')
        for i in images:
            if i != ref[0]:
                hdu = fits.open(i)
                data = hdu[0].data
                hdr = hdu[0].header
                mask = hdu[1].data
                data = np.array(data, dtype='float64')
                dx,dy,edx,edy = chi2_shift(ref_data, data, upsample_factor='auto')
                corrected_image = shift(data, [-1*dy, -1*dx])
                corrected_mask = shift(mask, [-1*dy, -1*dx])
                #delete any masked boundaries
                shape = corrected_image.shape
                rowNum_top = 0
                colNum_left = 0
                rowNum_bottom = 0
                colNum_right = 0
                for rows in range(shape[0]):
                    if (corrected_image[:(rows+1)]==0).all():
                        rowNum_top += 1
                    else:
                        break
                for cols in range(shape[1]):
                    if (corrected_image[:,:(cols+1)]==0).all():
                        colNum_left += 1
                    else:
                        break
                for rows2 in range(shape[0]):
                    if (corrected_image[((rows2+1)*-1):]==0).all():
                        rowNum_bottom += 1
                    else:
                        break
                for cols2 in range(shape[1]):
                    if (corrected_image[:,((cols2+1)*-1):]==0).all():
                        colNum_right += 1
                    else:
                        break
                rowNum_bottom *= -1
                colNum_right *= -1
                if rowNum_top == 0:
                    rowNum_top = None
                if colNum_left == 0:
                    colNum_left = None
                if rowNum_bottom == 0:
                    rowNum_bottom = None
                if colNum_right == 0:
                    colNum_right = None
                corrected_image = corrected_image[rowNum_top:rowNum_bottom,colNum_left:colNum_right]
                corrected_mask = corrected_mask[rowNum_top:rowNum_bottom,colNum_left:colNum_right]
                #write aligned array and mask to original image location
                hduData = fits.PrimaryHDU(corrected_image, header=hdr)
                hduMask = fits.ImageHDU(corrected_mask)
                hduList = fits.HDUList([hduData, hduMask])
                hduList.writeto(i[:-8] + '_A_.fits')
                os.remove(i)
                x += 1
                print("-> %.1f%% aligned..." % (float(x)/float(len(images)) * 100))
                hdu.close()
    else:
        print("-> Error: Alignment failed- reference image missing")
