#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 19:48:05 2018

@author: andrew
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy.special import iv
from astropy.io import fits
import glob
import sys
#import gc
from scipy.stats import chisquare
from scipy.stats import chi2
from scipy.stats import skellam as skell
from scipy.special import erf
from scipy.integrate import quad
from astropy.stats import sigma_clip
import sep

def stdev_dist(data, n_x=3, n_y=3):
    shape_x, shape_y = data.shape[0], data.shape[1]
#    while shape_x % n_x != 0:
#        n_x += 1
#    while shape_y % n_y != 0:
#        n_y += 1
#    stdevs = []
#    stamps_vert = np.split(data, n_x)
#    stamps_horiz = np.split(data, n_y, axis=1)
#    for x in stamps_vert:
#        stdevs.append(np.ma.median(x))
#    for y in stamps_horiz:
#        stdevs.append(np.ma.median(y))
#    return np.std(stdevs)
    width_x, width_y = np.floor(shape_x/n_x), np.floor(shape_y/n_y)
    d_x, d_y = np.floor(width_x/2), np.floor(width_y/2)
    x = []
    y = []
    stdevs = []
#    stamps = []
    for i in range(n_x):
        x.append(d_x + (i*width_x))
    for j in range(n_y):
        y.append(d_y + (j*width_y))
    for h in x:
        for k in y:
            stamp = data.data[int((k-d_y)):int((k+d_y+1)), int((h-d_x)):int((h+d_x+1))]
#            stamps.append(stamp)
            stdevs.append(np.ma.std(stamp))
    return np.std(stdevs)

def split_im(data, n_x=2, n_y=2):
#    shape_x, shape_y = data.shape[0], data.shape[1]
#    while shape_x % n_x != 0:
#        n_x += 1
#    while shape_y % n_y != 0:
#        n_y += 1
#    if n_x > 5 or n_y > 5:
#        print("-> Warning: odd image dimensions-- optimization may take a long time")
    stamps = []
    stamps_vert = np.array_split(data, n_x)
    stamps_horiz = np.array_split(data, n_y, axis=1)
    for v in stamps_vert:
        stamps_horiz = np.split(v, n_y, axis=1)
        for h in stamps_horiz:
            stamps.append(h)
    return stamps

def make_im(x=1000,y=1000,numStars=50,seeing=2,maxFlux=10000,bkg=50):
    Xs = np.random.randint(0,high=x,size=numStars)
    Ys = np.random.randint(0,high=y,size=numStars)
    fluxes = np.random.random(size=numStars)*maxFlux
    image = np.zeros((x,y))
    for i in range(numStars):
            image[Xs[i],Ys[i]] = fluxes[i]
    image = gaussian_filter(image,seeing)
    image += bkg
    return image

def add_poisson(image):
    im = image + np.random.poisson(image)
    return im

def plot_res(im1,im2):
    res = im1 - im2
    median1 = np.median(im1)
    median2 = np.median(im2)
    medianAvg = round(np.mean([median1,median2]))
    bns = int(round(2*(len(res.flatten()))**(1/3)))
    return plt.hist(res.flatten(), bins=bns, range=[-1*medianAvg,medianAvg], normed=True), medianAvg, bns

def get_res_data(image, N_x=2, N_y=2, opt=True):
    name = image.split('/')[-1]
    location = image.replace('/data/' + name,'')
    residual = location + '/' + 'residuals/' + name[:-5] + 'residual_.fits'
    template = glob.glob(location + '/templates/*.fits')
    if opt == True:
        res = '%s/temp/conv.fits' % (location)
        try: weight_check = fits.getval(res, 'WEIGHT')
        except: weight_check = 'N'
        if weight_check == 'N':
            res_data = np.ma.masked_array((fits.getdata(res))*-1, mask=fits.getdata(residual, 1))
        else:
            res_data = np.ma.masked_array((fits.getdata(res))*-1, mask=((fits.getdata(residual, 1)-1)*-1))
    else:
        try: weight_check = fits.getval(residual, 'WEIGHT')
        except: weight_check = 'N'
        if weight_check == 'N':
            res_data = np.ma.masked_array(fits.getdata(residual), mask=fits.getdata(residual, 1))
        else:
            res_data = np.ma.masked_array(fits.getdata(residual), mask=((fits.getdata(residual, 1)-1)*-1))
    res_data_clipped = sigma_clip(res_data.data, sigma=3)
    res_data_master_mask = np.logical_or(res_data.mask, res_data_clipped.mask)
    res_data = np.ma.MaskedArray(res_data.data, mask=res_data_master_mask)
    if len(template) == 1:
        im_data = fits.getdata(image)
        temp_data = fits.getdata(template[0])
        im_data = (im_data).byteswap().newbyteorder()
        temp_data = (temp_data).byteswap().newbyteorder()
        image_data = np.ma.masked_array(im_data, mask=np.logical_not(fits.getdata(image, 1)))
        template_data = np.ma.masked_array(temp_data, mask=np.logical_not(fits.getdata(template[0], 1)))
#        res_stamps = split_im(res_data, n_x=N_x, n_y=N_y)
        res_stamps = np.array_split(res_data, N_x)
        image_stamps = split_im(image_data, n_x=N_x, n_y=N_y)
        template_stamps = split_im(template_data, n_x=N_x, n_y=N_y)
#        try: image_median = float(fits.getval(image, 'MEDIAN'))
#        except: image_median = np.median(np.ma.MaskedArray(fits.getdata(image), mask=np.logical_not(fits.getdata(image, 1))))
#        try: template_median = float(fits.getval(template[0], 'MEDIAN'))
#        except: template_median = np.median(np.ma.MaskedArray(fits.getdata(template[0]), mask=np.logical_not(fits.getdata(template[0], 1))))
#        image_median = np.ma.median(image_data)
#        template_median = np.ma.median(template_data)
#        image_median = [np.ma.median(i) for i in image_stamps]
#        template_median = [np.ma.median(i) for i in template_stamps]
        image_backs = [sep.Background(i.data, mask=i.mask) for i in image_stamps]
        template_backs = [sep.Background(j.data, mask=j.mask) for j in template_stamps]
        image_median = [i.globalback for i in image_backs]
        template_median = [j.globalback for j in template_backs]
        del res_data, name, location, residual, template
        return res_stamps, image_median, template_median   
    elif len(template) > 1:
        print("\n-> Error: Too many template images in 'templates' directory\n")
        sys.exit()
    else:
        print("\n-> Error: Problem with number of template images\n")
        sys.exit()

#@profile
def get_res_hist(res_data, image_median_list, template_median_list):
    image_median = np.mean(image_median_list)
    template_median = np.mean(template_median_list)
    medianAvg = np.mean([image_median, template_median])
    res_median = np.ma.median(res_data)
    res_median_theoretical = image_median - template_median
    if np.ma.is_masked(np.ma.median(res_data)):
        res_median = 0
    rangeRes = abs((-1*5*medianAvg)+res_median_theoretical) + abs(res_median_theoretical+(5*medianAvg))
    range_low = res_median_theoretical - (5*abs(medianAvg))
    range_high = res_median_theoretical + (5*abs(medianAvg))
    bns = int(round((rangeRes/400)*2*(len((res_data.data).flatten()))**(1/3))/2)
#    if medianAvg > 0:
    return np.histogram(res_data.data, weights=(np.logical_not(res_data.mask)).astype(int), bins=bns, range=[range_low, range_high], density=True), res_median, np.ma.std(res_data), bns
#    elif medianAvg < 0:
#        return np.histogram(res_data.data, weights=(np.logical_not(res_data.mask)).astype(int), bins=bns, range=[2*medianAvg+res_median, res_median+(-1*2*medianAvg)], density=True), res_median, np.ma.std(res_data), bns
#    else:
#        return np.histogram(res_data.data, weights=(np.logical_not(res_data.mask)).astype(int), bins=bns, range=[-100,100], density=True), res_median, np.ma.std(res_data), bns

def plot_im(im,rnge=[0,500],bns=100):
    return plt.hist(im.flatten(), bins=bns, range=rnge)
   
def normal(x):
    return (1/(np.sqrt(2*np.pi)))*np.exp(-0.5*(x**2))

def int_gauss(x1, x2):
    I = quad(normal, x1, x2)
    return I[0]

def chi2_fit(res_data, median1, median2, conf_level=0.995, dx=100, dy=100, delta_df=2):
#    data_mid_horiz = data.shape[0] / 2
#    data_mid_vert = data.shape[1] / 2
#    res_data = data[(round(data_mid_vert)-dy):(round(data_mid_vert)+dy+1), (round(data_mid_horiz)-dx):(round(data_mid_horiz)+dx+1)]
    expected_iter = len(res_data.flatten())
    res_noise = np.sqrt(median1+median2)
    res_data = (res_data - (median1-median2)) / res_noise
#    expected_interval = skell.interval(0.99999999, median1, median2)
    medianAvg = np.mean([median1, median2])
#    data_range= abs((-1*2*medianAvg)+np.ma.median(res_data)) + abs(np.ma.median(res_data)+(2*medianAvg))
#    data_range = expected_interval[1] - expected_interval[0]
#    bns = int(round((data_range/400)*2*(len((res_data.data).flatten()))**(1/3))/2)
    bns = int((np.log2(expected_iter) + 1))
    expected_interval = (-5, 5)
#    expected = skellam(np.linspace(expected_interval[0], expected_interval[1], num=expected_iter), median1, median2)
#    expected = [skell.rvs(median1, median2) for i in range(expected_iter)]
#    expected = np.random.normal(size=expected_iter)
    if medianAvg >= 0:
        data_hist = np.histogram(res_data.data, weights=(np.logical_not(res_data.mask)).astype(int), bins=bns, range=expected_interval, density=False)
#        expected = [skellam(((5000)*np.random.random())-2.5, median1, median2) for i in range(expected_iter)]
#        expected_hist = np.histogram(expected, bins=bns, range=expected_interval, density=False)
    elif medianAvg < 0:
        data_hist = np.histogram(res_data.data, weights=(np.logical_not(res_data.mask)).astype(int), bins=bns, range=expected_interval, density=False)
#        expected = skellam(np.arange(hist_range[0], hist_range[1], expected_step_size), median1, median2)
#        expected_hist = np.histogram(expected, bins=bns, range=[2*medianAvg+np.ma.median(res_data), np.ma.mean(res_data)+(-1*2*medianAvg)], density=False)
    data_freq = data_hist[0]
    expected_hist = np.array([abs(int_gauss(data_hist[1][i], data_hist[1][i+1])) for i in range(len(data_hist[1])-1)])
    expected_freq = expected_hist * np.sum(data_freq)
#    expected_freq = expected_hist[0]
    data_freq = [data_freq[i] for i in range(len(data_freq)) if expected_freq[i] != 0]
    expected_freq = [expected_freq[i] for i in range(len(expected_freq)) if expected_freq[i] != 0]
#    expected_freq[expected_freq == 0] = 1 * 10**-16
    df = len(res_data.flatten()) - 1 - delta_df
    chi2_test = chisquare(data_freq, f_exp=expected_freq, ddof=delta_df)
    chi2_real = chi2.ppf(conf_level, df)
    Q = 1 / (1 + (chi2_test[0] / chi2_real))
    b_mid = 0.5*(data_hist[1][1:]+data_hist[1][:-1])
    b_mid = b_mid.astype(np.float64)
    plt.plot(b_mid, data_freq, 'b-')
    plt.plot(b_mid, expected_freq, 'r--')
    return Q, medianAvg, bns, data_hist, expected_hist, expected_freq, chi2_test, chi2_real
    
#@profile
def fit(image, n_x=1, n_y=1, dist='skellam', optimize=True, qThresh=0.35):
#    signal.alarm(10)
    Qs = []
    res_stamps, median1, median2 = get_res_data(image, N_x=n_x, N_y=n_y, opt=optimize)
    for i in range(len(res_stamps)):
        try:
            hist, dataMean, stdev, bns = get_res_hist(res_stamps[i], median1, median2)
            if stdev >= 600 or stdev == 0:
                Q = 0
            else:
                a, b = hist[0], hist[1]
                b_mid = 0.5*(b[1:]+b[:-1])
                b_mid = b_mid.astype(np.float64)
        #        name = image.split('/')[-1]
        #        location = image.replace('/data/' + name,'')
        #        res = location + '/' + 'residuals/' + name[:-5] + 'residual_.fits'
        #        res_data = fits.getdata(res)
        #        if optimize == True:
        #            res_data = fits.getdata(location + '/temp/' + 'conv.fits')
            #    bZeros = [b_mid[i] for i in range(len(a)) if 10 > (a[i]*len(im1.flatten()))]
            #    a = [i for i in a if 10 < (i*len(im1.flatten()))]
            #    for z in bZeros:
            #        b_mid = list(filter((z).__ne__, b_mid))
            #    a = np.array(a, dtype='float64')
            #    b_mid = np.array(b_mid, dtype='float64')
#                if medianAvg > 0:
#                xdata = np.linspace((-1*2*np.mean([median1, median2]))+np.median(res_stamps[i]),np.median(res_stamps[i])+(2*np.mean([median1, median2])),10*bns)
#                elif medianAvg < 0:
#                    xdata = np.linspace((2*medianAvg)+np.median(res_data),np.median(res_data)+(-1*2*medianAvg),10*bns)
#                else:
#                    xdata = np.linspace(-100,100,10*bns)
                if dist == 'skellam':
                    params, covars = curve_fit(skellam,b_mid,a,[median1[i], median2[i]])
#                    plt.plot(xdata, skellam(xdata,*params),'r-')
#                    plt.plot(b_mid,a,'b-')
                    residuals = a - skellam(b_mid,*params)
                if dist == 'gauss':
                    params, covars = curve_fit(gauss,b_mid,a)
        #            plt.plot(xdata, gauss(xdata,*params),'r-')
                    residuals = a - gauss(b_mid,*params)
                residuals = residuals[~np.isnan(residuals)]
                residuals = residuals[~np.isinf(residuals)]
                residuals = np.array(residuals)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((a-np.mean(a))**2)
                rSquared = (1-(ss_res/ss_tot))
                res_fitMean = abs(params[0] - params[1])
                res_realMean = abs(median1[i] - median2[i])
                Q = rSquared/(1+abs(res_fitMean-res_realMean))
                print(rSquared, res_fitMean, res_realMean)
                if list(residuals) == []:
                    Q = 0
        #        plt.plot(b_mid,a)
                del params, residuals, res_fitMean
        except Exception as e:
                Q = 0
                print(e)
#    signal.alarm(0)
        Qs.append(Q)
    if Qs == []:
        Q = 0
    if np.min(Qs) < qThresh:
        Q_final = np.min(Qs)
    else:
        Q_final = np.mean(Qs)
    del Qs, res_stamps, median1, median2
    try:
        del hist,  dataMean, stdev, residuals, a, b, b_mid
    except:
        pass
    return Q_final 


def skellam(k,u1,u2):
    return np.exp(-1*(u1+u2))*(u1/u2)**(k/2)*iv(k, 2*np.sqrt(u1*u2))

def gauss(x,mean,sig):
    return (1/np.sqrt((2*np.pi*sig**2)))*np.exp(-1*(x-mean)**2/(2*sig**2))
