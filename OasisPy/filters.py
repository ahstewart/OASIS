#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 18:28:42 2018

@author: andrew
"""

from astropy.io import fits
import numpy as np
import glob
import copy
import sys
import psf
import datetime
import initialize
import os
import sex
from tqdm import tqdm
import math
from scipy.stats import skellam
import sep

def get_sources(image, filtered=True, MR=False, phose=False):
    '''
    gets all the point sources and fluxes deteced with SExtractor on a certain residual image
    the format of the outputted data is
    first column: flux
    second: x pixel position 
    third: y pixel position
    '''
    image_name = image.split("/")[-1]
    image_name_res = image_name[:-5] + 'residual_.fits'
    length = (len(image_name)+6)*-1
    location = image[:length]
    filt_source = location + "/sources/filtered_sources.txt"
    if phose == True:
        filt_source = "%s/data/phose_sources.txt" % (location)
        temp_source = "%s/templates/phose_template.cat" % (location)
        image_name_res = image_name
    elif phose == False:
        if filtered == False:
            filt_source = location + "/sources/sources.txt"
        if MR == True and filtered == False:
            location = image.split('/')[:-1]
            location = '/'.join(location)
            location = location[:-9]
            filt_source = location + "sources/MR_sources.txt"
            image_name_res = image.split('/')[-1]
        if MR == True and filtered == True:
            location = image.split('/')[:-1]
            location = '/'.join(location)
            location = location[:-9]
            filt_source = location + "sources/MR_sources_filtered.txt"
            image_name_res = image.split('/')[-1]
    else:
        print("-> Error: 'phose' keyword must be boolean\n-> Exiting...")
        sys.exit()
    data = []
    inds = []
    with open(filt_source, 'r') as filt:
        sources = filt.readlines()
        filt.close()
    for s in sources:
        if s == (image_name_res+'\n'):
            check = True
            x = sources.index(s)
            while check == True:
                x += 1
                try:
                    float(sources[x].split()[0])
                    data.append(sources[x].split())
                    inds.append(x)
                except:
                    try: sources[x].split()[0]
                    except: check = False
            break
    for d in range(len(data)):
        if phose == True:
            data[d] = data[d][1:]
        else:
            data[d] = data[d][1:-1]
        for i in range(len(data[d])):
            data[d][i] = float(data[d][i])
    if phose == True:
        temp_data = []
        with open(temp_source, 'r') as temp:
            templines = temp.readlines()
        for t in templines:
            check = True
            x2 = templines.index(t)
            while check == True:
                x2 += 1
                try:
                    float(templines[x2].split()[0])
                    temp_data.append(templines[x2].split())
                except:
                    try: templines[x2].split()[0]
                    except: check = False
            break
        for td in range(len(temp_data)):
            temp_data[td] = temp_data[td][1:]
            for i in range(len(temp_data[td])):
                temp_data[td][i] = float(temp_data[td][i])
        return data, temp_data
    else:
        return data, inds

def get_all_sources(location, filt=True):
    images = glob.glob(location + '/data/*.fits')
    sources = []
    indices = []
    for i in images:
        d,ind = get_sources(i, filtered=filt)
        sources.append(d)
        indices.append(ind)
    return sources, indices

def source_count(location, filtered=True):
    sources, indices = get_all_sources(location, filt=filtered)
    xyCoords = []
    finalSources = []
    for s in sources:
        for i in range(len(s)):
            xyCoords.append((round(s[i][1]), round(s[i][2])))
    xyCoords_copy = copy.deepcopy(xyCoords)
    for xy in xyCoords_copy:
        xyCoords.remove(xy)
        if xy not in xyCoords and (xy[0]+1,xy[1]) not in xyCoords and (xy[0]-1,xy[1]) not in xyCoords and (xy[0],xy[1]+1) not in xyCoords and (xy[0],xy[1]-1) not in xyCoords:
            finalSources.append(xy)
    return finalSources, len(xyCoords_copy)

def get_image_names(location, filtered=False):
    source_loc = "%s/sources/sources.txt" % (location)
    if filtered == True:
        source_loc = "%s/sources/filtered_sources.txt" % (location)
    try:
        with open(source_loc, 'r') as src:
            src_lines = src.readlines()
            image_names = []
        for s in src_lines:
            if len(s.split()) == 1:
                image_names.append(s)
        return image_names
    except FileNotFoundError:
        return []

def write_total_sources(location):
    Q_min = float(initialize.get_config_value('qFloor'))
    print("\n-> Calculating detection statistics...\n")
    uniqueSources, numFilteredSources = source_count(location)
    originalSources, numSources = source_count(location, filtered=False)
    MR_sources, MR_inds = get_sources("%s/residuals/MR.fits" % (location), MR=True, filtered=False)
    MR_sources_filt, MR_inds_filt = get_sources("%s/residuals/MR.fits" % (location), MR=True, filtered=True)
    MR_sources = len(MR_sources)
    MR_sources_filt = len(MR_sources_filt)
    total_source_loc = location + '/sources/total_sources.txt'
    residuals = glob.glob(location + '/residuals/*_residual_.fits')
    bad_subtractions = 0
    date = datetime.datetime.now()
    for r in residuals:
        if (fits.getdata(r, 1) == 0).all():
            bad_subtractions += 1
    with open(total_source_loc, 'a+') as total:
        total.write('Date Run: %d/%d/%d %d:%d:%d | Number of Images Subtracted = %d\n' % (date.month, date.day, date.year, date.hour, date.minute, date.second, len(residuals)))
        total.write('Total Initial Sources: %d\n' % (numSources))
        print('\nTotal Initial Sources: %d\n' % (numSources))
        total.write('Total Filtered Sources: %d\n' % (numFilteredSources))
        print('Total Filtered Sources: %d\n' % (numFilteredSources))
        total.write('Total Unique Detections: %d\n' % (len(uniqueSources)))
        print('Total Unique Detections: %d\n' % (len(uniqueSources)))
        total.write('\nTotal Master Residual Sources: %d\n' % (MR_sources))
        print('\nTotal Master Residual Sources: %d\n' % (MR_sources))
        total.write('Total Filtered Master Residual Sources (best representation of real # of sources): %d\n' % (MR_sources_filt))
        print('Total Filtered Master Residual Sources (best representation of real # of sources): %d\n' % (MR_sources_filt))
        total.write('\nBad Subtractions (Q-Value < 0.50): %d/%d' % (bad_subtractions, len(residuals)))
        print('\nBad Subtractions (Q-Value < %.2f): %d/%d\n' % (Q_min, bad_subtractions, len(residuals)))
        total.write('\nAverage Number of Sources Per Image: %.2f\n\n\n' % (float(numFilteredSources/len(residuals))))
        print('\nAverage Number of Sources Per Image: %.2f\n\n\n' % (float(numFilteredSources/len(residuals))))
        total.close()
    print("\n-> Complete!\n")
        
def reoccuring(location, pix_dist=1.5, use_config_file=True):
    if use_config_file == True:
        pix_dist = initialize.get_config_value('pix_dist', file_loc=location+'/configs')
    sources,indices = get_all_sources(location)
    for i in range(len(sources)):
        new_sources,new_indices = get_all_sources(location)
        del_inds = []
        for j in range(len(new_sources[i])):
            inds = []
            x_low = round(new_sources[i][j][1]) - pix_dist
            x_high = round(new_sources[i][j][1]) + pix_dist
            y_low = round(new_sources[i][j][2]) - pix_dist
            y_high = round(new_sources[i][j][2]) + pix_dist
            check = 0
            for h in range(len(new_sources)):
                if h != i:
                    for k in range(len(new_sources[h])):
                        x = new_sources[h][k][1]
                        y = new_sources[h][k][2]
                        if x_low < x < x_high and y_low < y < y_high:
                            if check == 0:
                                inds.append(new_indices[i][j])
                            inds.append(new_indices[h][k])
                            check += 1
            if len(inds) >= (len(sources)/2):
                for index in inds:
                    del_inds.append(index)
        update_filtered_sources(location, del_inds)
            
def divot(image, nx=30, ny=30, MR=False):
    #this checks for 'diveted' detections usually indicative of errors in the
    #AIS psf convolution
    ind = []
    if MR == True:
        location = image.split('/')[:-1]
        location = '/'.join(location)
        location = location.replace('/residuals', '')
        image_res = image
        image_sources, indices = get_sources("%s/residuals/MR.fits" % (location), MR=True, filtered=True)
    else:
        image_res = image.replace('data','residuals')
        image_res = image_res[:-5] + 'residual_.fits'
        image_sources, indices = get_sources(image)
    image_data = fits.getdata(image_res)
    sigma = np.std(image_data)
    for s in image_sources:
        x, y = s[1], s[2]
        stamp = image_data[(round(y)-ny):(round(y)+ny+1), (round(x)-nx):(round(x)+nx+1)]
        f1 = (stamp<(sigma*-2)).sum()
        f2 = (stamp<(sigma*-4)).sum()
        f3 = (stamp<(sigma*-6)).sum()
        if f1>200 or f2>20 or f3>2:
            ind.append(indices[image_sources.index(s)])
    return ind

def update_filtered_sources(location, inds, MR=False):
    filt_source = location + "/sources/filtered_sources.txt"
    if MR == True:
        filt_source = location + "/sources/MR_sources_filtered.txt"
    upd_sources = []
    with open(filt_source, 'r') as filt:
        sources = filt.readlines()
        filt.close()
    for i in inds:
        upd_sources.append(sources[i])
    for s in upd_sources:
        sources = list(filter((s).__ne__, sources))
    with open(filt_source, 'w+') as filt:
        filt.writelines(sources)
        filt.close()    

def spread_model_filter(location, spread_model_min=-0.025, spread_model_max=0.1, MR=False, use_config_file=True):
    #filter source file by spread_model and puts the results in filtered_sources.txt
    source_loc = location + '/sources'
    source_txt_loc = source_loc + '/sources.txt'
    source_txt_filtered_loc = source_loc + '/filtered_sources.txt'
    if MR == True:
        source_txt_loc = source_loc + '/MR_sources.txt'
        source_txt_filtered_loc = source_loc + '/MR_sources_filtered.txt'
    del_lin = []
    with open(source_txt_loc, 'r') as src:
        lines = src.readlines()
        src.close()
    if use_config_file == True:
        spread_model_min = initialize.get_config_value('spread_model_min', file_loc=location+'/configs')
        spread_model_max = initialize.get_config_value('spread_model_max', file_loc=location+'/configs')
    for lin in lines:
        parse = lin.split()
        if parse != []:
            try:
                int(parse[0])
                if float(parse[-1]) < spread_model_min or float(parse[-1]) > spread_model_max:
                    del_lin.append(lin)
            except ValueError or IndexError:
                pass
    lines = [a for a in lines if a not in del_lin]
    with open(source_txt_filtered_loc, 'w+') as fil_src:
        fil_src.writelines(lines)
        fil_src.close()
        
def mask_sources_image(res_image, aperture_diam=1.5, use_config_file=True):
    if use_config_file == True:
        location = res_image.split('/')[:-2]
        location = '/'.join(location)
        aperture_diam = initialize.get_config_value('aperture_diam', file_loc=location+'/configs')
    res_data = fits.getdata(res_image)
    res_mask = fits.getdata(res_image, 1)
    weight_check = False
    if fits.getval(res_image, 'WEIGHT') == 'Y':
        weight_check = True
        res_mask = (res_mask-1) * -1
    image = res_image.replace('_residual', '')
    image = image.replace('residuals', 'data')
    im_fwhm = psf.fwhm(image)
    unfiltered_sources, unfiltered_inds = get_sources(image, filtered=False)
    filtered_sources, filtered_inds = get_sources(image, filtered=True)
    for unf in unfiltered_sources:
        if unf not in filtered_sources:
            new_mask = mask_source(res_data.shape[0], res_data.shape[1], (unf[1], unf[2]), aperture_diam*im_fwhm)
            res_mask = np.logical_or(res_mask, new_mask)
    data_hdu = fits.PrimaryHDU(res_data, header=fits.getheader(res_image))
    if weight_check == True:
        mask_hdu = fits.ImageHDU((res_mask-1) * -1)
    else:
        mask_hdu = fits.ImageHDU(res_mask)
    list_hdu = fits.HDUList([data_hdu, mask_hdu])
    list_hdu.writeto(res_image, overwrite=True)
        
def mask_source(x_dim, y_dim, centroid, radius, method='disk'):
    y,x = np.ogrid[-centroid[0]:x_dim-centroid[0], -centroid[1]:y_dim-centroid[1]]
    if method == 'disk':
        mask = x*x + y*y <= radius*radius
    elif method == 'ring':
        mask = x*x + y*y == (radius*radius)
    else:
        print("-> Error: Unkown method entered\n-> Exiting...")
        sys.exit()
    return mask.astype(int)

def stamp(data, x, y, x_width, y_width):
    return data[y-y_width:y+y_width+1, x-x_width:x+x_width+1]
            
def phose_sex(location):
    """
    create source catalogs for each science image using SExtractor, combine them
    into a master catalog called 'phose_sources.txt' in the data directory
    
    also creates a template catalog used for phose thresholding
    """
    
    initialize.create_configs(location)
    try: os.mkdir("%s/data/cats" % (location))
    except: pass
    images = glob.glob("%s/data/*.fits" % (location))
    config = "%s/configs/phose_default.sex" % (location)
    param = "%s/configs/phose_default.param" % (location)
    for i in tqdm(images):
        output_cat = i.replace('data', 'data/cats')
        output_cat = output_cat.replace('fits', 'txt')
        with open(config, 'r') as conf:
            conf_lines = conf.readlines()
        conf_lines[9] = "PARAMETERS_NAME" + "        " + param + "\n"
        conf_lines[32] = "WEIGHT_IMAGE" + "        " + "%s[1]" % (i) + "\n"
        with open(config, 'w') as conf:
            conf.writelines(conf_lines)
        os.system("sextractor %s[0] > %s -c %s" % (i, output_cat, config))
    sex.src_join(location, phose=True)
    template = glob.glob("%s/templates/*.fits" % (location))
    if template == []:
        print("-> Error: template not found\n-> Run combine.py first\n-> Exiting...")
        sys.exit()
    template = template[0]
    output_template_cat = "%s/templates/phose_template.cat" % (location)
    with open(config, 'r') as conf:
        conf_lines = conf.readlines()
    conf_lines[9] = "PARAMETERS_NAME" + "        " + param + "\n"
    conf_lines[32] = "WEIGHT_IMAGE" + "        " + "%s[1]" % (template) + "\n"
    with open(config, 'w') as conf:
        conf.writelines(conf_lines)
    os.system("sextractor %s[0] > %s -c %s" % (template, output_template_cat, config))
    
def phose(image, dpixel=1, thresh=3.5, kron_min=0.01, fill_method='gauss', negative=False, fit='moffat'):
    data_image = image.replace('residual_', '')
    data_image = data_image.replace('residuals', 'data')
    res_data = fits.getdata(image)
    if negative == True:
        res_data *= -1
    res_mask = fits.getdata(image, 1)
    try: weight_check = fits.getval(image, 'WEIGHT')
    except: weight_check = 'N'
    if weight_check == 'Y':
        res_mask = (res_mask-1)*-1
    location = image.split('/')[:-2]
    location = '/'.join(location)
    template = glob.glob("%s/templates/*.fits" % (location))[0]
    try: template_median = float(fits.getval(template, 'MEDIAN'))
    except: template_median = np.median(fits.getdata(template))
    try: science_median = float(fits.getval(data_image, 'MEDIAN'))
    except: science_median = np.median(fits.getdata(data_image))
    res_data_sep = res_data.byteswap().newbyteorder()
    try: res_bkg = sep.Background(res_data_sep, mask=res_mask)
    except ValueError: res_bkg = res_bkg = sep.Background(res_data, mask=res_mask)
    res_rms = res_bkg.globalrms
    res_back = res_bkg.globalback
    if fill_method == 'gauss':
        fill_bkg = np.random.normal(loc=res_bkg.globalback, scale=res_bkg.globalrms, size=res_data.shape)
    elif fill_method == 'skellam':
        fill_bkg = skellam.rvs(float(science_median), float(template_median), size=(res_data.shape))
        fill_bkg = fill_bkg.astype(np.float64)
    else:
        print("-> Error: Invalid value for 'fill_method' keyword\n-> Exiting...")
        sys.exit()
    FWHM = psf.fwhm(data_image)
    sigma = FWHM/2.355
    if fit == 'moffat':
        fit_param = moffat_fwhm_to_a(FWHM)
    elif fit == 'gauss':
        fit_param = FWHM/2.355
    else:
        print("-> Error: Invalid value for 'fit' parameter\n-> Exiting...")
        sys.exit()
    unfiltered_sources, temp_sources = get_sources(data_image, filtered=False, phose=True)
    filtered_sources, filtered_inds = get_sources(data_image, filtered=True)
    bad_mask = np.zeros(res_data.shape)
    good_mask = np.zeros(res_data.shape)
    for s in unfiltered_sources:
        og_flux = s[0]
        x = s[2]
        y = s[3]
        kron_radius = s[1]
        a_image = s[4]
        b_image = s[5]
        min_flux = (np.pi*np.mean((kron_radius*a_image, kron_radius*b_image))**2) * res_back
        if (og_flux < min_flux) or (kron_radius == 0):
            continue
        theta_image = s[6]
        if kron_radius < kron_min:
            indices = [0]
        else:
            indices = [v[1] for i, v in enumerate(temp_sources) if (v[2]-dpixel<=x<=v[2]+dpixel and v[3]-dpixel<=y<=v[3]+dpixel)
                        and phose_check(res_data, x, y, kron_radius, a_image, b_image, theta_image, thresh, og_flux, v[0])]
#        phoseCheck = phose_check(res_data, x, y, kron_radius, a_image, b_image, theta_image, thresh, og_flux, og_flux)
#        if phoseCheck == True or (s not in filtered_sources):
        if indices != []:
            if np.mean(indices) > kron_radius:
                kron_radius = np.mean(indices)
            kron_radius *= 1.5
            size_check = aperture_resize(res_rms, fit_param, sigma, og_flux, kron_radius, x, y, a_image, b_image, dist=fit)
            while size_check == True:
                kron_radius *= 1.5
                size_check = aperture_resize(res_back, fit_param, sigma, og_flux, kron_radius, x, y, a_image, b_image, dist=fit)
            ellipse_mask(bad_mask, (y-1), (x-1), kron_radius, a_image, b_image, theta_image, fill_value=1)
        else:
            ellipse_mask(good_mask, (y-1), (x-1), kron_radius, a_image, b_image, theta_image, fill_value=1)
#            res_data = phose_fill(res_data, phose_check(res_data, x, y, kron_radius, a_image, b_image, theta_image, mask_return=True), fill_bkg)
#        for t in temp_sources:
#            t_flux = t[0]
#            t_x = t[2]
#            t_y = t[3]
#            if (t_x-dpixel<=x<=t_x+dpixel) and (t_y-dpixel<=y<=t_y+dpixel):
#                temp_mask = np.ones(res_data.shape)
#                ellipse_mask(temp_mask, (y-1), (x-1), kron_radius, a_image, b_image, theta_image, fill_value=0)
#                res_data_phot = np.ma.MaskedArray(res_data, mask=temp_mask, fill_value=0)
#                s_flux = np.sum(res_data_phot.filled())
#                if s_flux < (thresh * np.sqrt(og_flux + t_flux)):
#                    temp_mask_fill = (temp_mask - 1) * -1
#                    res_data_final = np.ma.MaskedArray(res_data, mask = temp_mask_fill, fill_value=0)
#                    res_data_final = res_data_final.filled()
#                    phose_patch = np.multiply(temp_mask_fill, fill_bkg)
#                    res_data = res_data_final + phose_patch
#                    del temp_mask, res_data_phot, temp_mask_fill, res_data_final
#    bad_hdu = fits.PrimaryHDU(bad_mask.astype(np.float64))
#    good_hdu = fits.PrimaryHDU(good_mask.astype(np.float64))
#    bad_hdu.writeto("%s/%s_bad.fits" % (location, (image.split('/'))[-1]))
#    good_hdu.writeto("%s/%s_good.fits" % (location, (image.split('/'))[-1]))
    and_mask = np.logical_and(bad_mask, good_mask)
    bad_mask -= and_mask
    bad_mask = np.logical_or(bad_mask, res_mask)
    res_data = phose_fill(res_data, bad_mask, fill_bkg)
#    res_hdu = fits.PrimaryHDU(res_data)
#    res_hdu.writeto("%s/%s_res.fits" % (location, (image.split('/'))[-1]))
#    bad_hdu = fits.PrimaryHDU(bad_mask.astype(np.float64))
#    bad_hdu.writeto("%s/%s.fits" % (location, (image.split('/'))[-1]))
    return (res_data-(science_median-template_median))/np.sqrt(science_median + template_median)

def phose_check(res_data, x, y, kron, a, b, theta, thresh=0, og_flux=0, t_flux=0, mask_return=False):
    temp_mask = np.ones(res_data.shape)
    ellipse_mask(temp_mask, (y-1), (x-1), kron, a, b, theta, fill_value=0)
    res_data_phot = np.ma.MaskedArray(res_data, mask=temp_mask, fill_value=0)
    s_flux = np.sum(res_data_phot.filled())
    del res_data_phot
    if mask_return == False:
        if abs(s_flux) < (thresh * np.sqrt(abs(og_flux + t_flux))):
            return True
        else:
            return False
    else:
        return temp_mask

def phose_fill(res_data, temp_mask, bkg):
    res_data_final = np.ma.MaskedArray(res_data, mask=temp_mask, fill_value=0)
    res_data_final = res_data_final.filled()
    phose_patch = np.multiply(temp_mask, bkg)
    res_data = res_data_final + phose_patch
    del temp_mask, res_data_final
    return res_data

def ellipse_mask(array, x_centroid, y_centroid, kron_radius, a, b, theta, fill_value=1):
    X = array.shape[0]
    Y = array.shape[1]
    theta = math.radians(theta-90)
    a = a * kron_radius
    b = b * kron_radius
    if a != 0 and b != 0:
        x,y = np.ogrid[-x_centroid:X-x_centroid, -y_centroid:Y-y_centroid]
        mask = ((x*np.cos(theta)+y*np.sin(theta))**2/a**2) + ((x*np.sin(theta)-y*np.cos(theta))**2/b**2) < 1
        array[mask] = fill_value 
        return mask
    else:
        pass

def aperture_resize(bkg, param, sigma, flux, kron_r, x, y, a, b, dist='moffat'):
    amplitude = get_gaussian_amplitude(flux, sigma)
    mean_radius = np.mean((kron_r*a, kron_r*b))
    if dist == 'moffat':
        check = abs(moffat(mean_radius, amplitude, param)) > abs(bkg)
    elif dist == 'gauss':
        check = abs(gauss(mean_radius, amplitude, param)) > abs(bkg)
    else:
        print("-> Error: Invalid value for 'dist' parameter\n-> Exiting...")
        sys.exit()
    if check == True:
        return True
    else:
        return False
    
def gauss(x, amplitude, sigma, mean=0):
    return (amplitude * (np.exp(-((x-mean)**2)/(2*(sigma**2)))))

def moffat(x, amplitude, a, b=4.765):
    return amplitude * (1 + (x/a)**2)**(-b)

def moffat_fwhm_to_a(fwhm, b=4.765):
    return fwhm / (2 * np.sqrt(2**(1/b) - 1))    

def get_gaussian_amplitude(flux, sigma):
    return (flux / (2 * np.pi * (sigma**2)))

def check(flux1, flux2):
    res = abs(flux1 - flux2)
    poiss1 = np.sqrt(flux1 + flux2)
    poiss2 = np.sqrt(flux1 * 2)
    print(res, poiss1, poiss2)
    if res > poiss1:
        print("passes test 1")
    if res > poiss2:
        print("passes test 2")