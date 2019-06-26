#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 15:55:42 2018

@author: andrew
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from astropy.io import fits

data = fits.getdata('/home/andrew/OASIS/targets/M31_S13_C_R1/00:44:03.710_+41:14:46.19/air/10/temp/03:14:29.335_A_.fits')

img_shape = data.shape

n_stars = 3
pos_x = [1500,2000,1200]
pos_y = [1600,1400,2200]
array = np.array([ 0.65343465,  0.50675629,  0.84946314])
fluxes = 200000.0 + array * 300.0
img = np.zeros(img_shape)
for x, y, f in zip(pos_x, pos_y, fluxes):
    img[x, y] = f

img = gaussian_filter(img, sigma=15.0, mode='constant')

data += img

hdu = fits.PrimaryHDU(data)
hdu.writeto('/home/andrew/OASIS/targets/M31_S13_C_R1/00:44:03.710_+41:14:46.19/air/10/temp/03:14:29.335_A_.fits', overwrite=True)

plt.imshow(img)