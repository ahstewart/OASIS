#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 19:20:13 2018

@author: andrew
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.mlab as mlab
import math

def hist(image):
    hdu = fits.open(image)
    data = hdu[0].data
#    N = poiss("/home/andrew/ISIS_images/images/07:41:27.555_A_.fits", "/home/andrew/ISIS_images/images/08:00:03.041_A_.fits")
#    data = data/N
    data = data/np.sqrt(np.var(data)+data)
    plt.hist(data.flatten(), bins=50, range=(-5,5))
    mu = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    plt.plot(x,mlab.normpdf(x, mu, sigma)*800000)
    plt.show()

def dates(location):
    images = glob.glob(location + "/*.fits")
    image_names = []
    for i in images:
        image_names.append(i[(len(location)+1):])
    dates = location[:-4] + "/register/dates"
    with open(dates, 'w') as date:
        for im in image_names:
            date.write(im + "\n")
        
def var(image):
    data = fits.getdata(image)
    var = np.var(data)
    return var

def poiss(image, reference_image):
    V = var(reference_image)
    data = fits.getdata(image)
    norm = np.ones(data.shape)
    for i in range(np.size(data, axis=0)):
        for j in range(np.size(data, axis=1)):
            norm[i][j] = np.sqrt(V+data[i][j])
    return norm