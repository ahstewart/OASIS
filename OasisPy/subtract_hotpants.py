#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 13:04:11 2018

@author: andrew
"""

import glob
import os
import psf
import numpy as np
import sys
from astropy.io import fits
import initialize

def hotpants(location):
    '''
    subtract using hotpants
    '''
    images = glob.glob(location + "/data/*_A_.fits")
    template = glob.glob(location + "/templates/*.fits")
    residuals = glob.glob(location + "/residuals/*residual_.fits")
    psf_data = glob.glob(location + '/psf/*')
    images_names = [(i.split('/')[-1])[:-5] for i in images]
    res_names = [(r.split('/')[-1])[:-14] for r in residuals]
    resids = [res for res in images_names if res not in res_names]
    ims = []
    for rs in resids:
        ims.append(location+'/data/'+rs+'.fits')
    if ims != []:
        if len(psf_data) == 2*(len(images)+1):
            if len(template) == 1:
                for i in ims:
                    hotpants_subtract(i, default=True)
            else:
                print("-> Error with number of templates")
                sys.exit()
        else:
            print("-> Error: Need PSFs before running subtraction\n-> Run psf.py first")
            print("-> If any images have been manually removed from the data directory, delete all contents of the psf directory and run OasisPy again\n")
            sys.exit()
    else:
        print("-> Images have already been subtracted")
        
def hotpants_subtract(image, default=True, opt_ind=0, spatial_deg=0, params=None):
    location = image.split('/')[:-2]
    location = '/'.join(location)
    name = (image.split('/')[-1])[:-5]
    template = glob.glob("%s/templates/*.fits" % location)
    template = template[0]
    image_header = fits.getheader(image)
    sig_image = psf.fwhm(image)/2.355
    template_header = fits.getheader(template)
    sig_template = psf.fwhm_template(template)/2.355
    tu = template_header['SATURATE']
    tg = template_header['GAIN']
    tmi = "%s/hotpants_mask.fits" % location
    if os.path.exists(tmi) == False:
        try: weight_check = fits.getval(template, 'WEIGHT')
        except: weight_check = 'N'
        temp_mask = fits.getdata(template, 1)
        if weight_check == 'Y':
            temp_mask = (temp_mask - 1) * -1
        hotpants_hdu = fits.PrimaryHDU(temp_mask)
        hotpants_hdu.writeto("%s/hotpants_mask.fits" % (location), overwrite=True)
    iu = image_header['SATURATE']
    ig = image_header['GAIN']
    try: weight_check = fits.getval(image, 'WEIGHT')
    except: weight_check = 'N'
    if weight_check == 'Y':
        im_mask = fits.getdata(image, 1)
        im_mask = (im_mask-1)*-1
        im_mask_hdu = fits.PrimaryHDU(im_mask)
        im_mask_hdu.writeto("%s/configs/temp_htpts_mask.fits" % (location), overwrite=True)
        imi = "%s/configs/temp_htpts_mask.fits" % (location)
    else:
        imi = "%s[1]" % (image)
    cwd = os.getcwd()
    os.chdir(os.path.dirname(initialize.__file__))
    if default == True:
        omi = "%s/residuals/BPM.fits" % location
        outim = "%s/residuals/%sresidual_.fits" % (location, name)
        tl = -1000
        il = -1000
        if sig_template < sig_image:
            sigma_match = np.sqrt((sig_image)**2-(sig_template)**2)
            s1 = .5*sigma_match
            s2 = sigma_match
            s3 = 2*sigma_match
            os.system("./hotpants -inim %s -tmplim %s -outim %s -tl %d -il %d -tu %d -tg %d -tmi %s -iu %d -ig %d -imi %s -omi %s -ng 3 6 %.5f 4 %.5f 2 %.5f -bgo 0" % (image, template, outim, tl, il, tu, tg, tmi, iu, ig, imi, omi, s1, s2, s3))
        elif sig_template >= sig_image:
            sigma_match = np.sqrt((sig_template)**2-(sig_image)**2)
            s1 = .5*sigma_match
            s2 = sigma_match
            s3 = 2*sigma_match
            os.system("./hotpants -inim %s -tmplim %s -outim %s -tl %d -il %d -tu %d -tg %d -tmi %s -iu %d -ig %d -imi %s -omi %s -ng 3 6 %.5f 4 %.5f 2 %.5f -bgo 0" % (template, image, outim, il, tl, iu, ig, imi, tu, tg, tmi, omi, s1, s2, s3))
            invert(outim)
        res_data = fits.getdata(outim)
        res_header = fits.getheader(outim)
        BPM = fits.getdata(omi)
        BPM_median = np.median(BPM)
        BPM_std = np.std(BPM)
        BPM[BPM > (BPM_median + BPM_std)] = 0
        BPM[BPM <= (BPM_median + BPM_std)] = 1
        resHDU = fits.PrimaryHDU(res_data, header=res_header)
        resMask = fits.ImageHDU(BPM)
        resList = fits.HDUList([resHDU, resMask])
        resList.writeto(outim, overwrite=True)
        os.remove(omi)
        os.remove("%s/configs/temp_htpts_mask.fits" % location)
        hdu = fits.open("%s/residuals/%sresidual_.fits" % (location, name), mode='update')
        hdr = hdu[0].header
        hdr.set('OPTIMIZE', 'N')
        hdu.close()
    elif default == False:
        omi = "%s/temp/BPM.fits" % location
        outim = "%s/temp/conv.fits" % location
        if 0 <= opt_ind <= 1:
            s1, s2, s3 = 0.7, 1, 2.5
        elif opt_ind == 2:
            s1, s2, s3 = 0.7, 1.5, 3
        elif opt_ind >= 3:
            s1, s2, s3 = 0.7, 2, 4
        if sig_template < sig_image:
            os.system("./hotpants -inim %s -tmplim %s -outim %s -tl %d -il %d -tu %d -tg %d -tmi %s -iu %d -ig %d -imi %s -omi %s -ng 3 6 %.5f 4 %.5f 2 %.5f -bgo %d -r %d -nsx %d -nsy %d -rss %d -ko %d" % (image, template, outim, params['tl'][0], params['il'][0], tu, tg, tmi, iu, ig, imi, omi, s1, s2, s3, params['bgo'][0], (params['r'])[opt_ind], (params['nsx'])[spatial_deg], (params['nsy'])[spatial_deg], (params['rss'])[opt_ind], (params['ko'])[spatial_deg]))
        else:
            os.system("./hotpants -inim %s -tmplim %s -outim %s -tl %d -il %d -tu %d -tg %d -tmi %s -iu %d -ig %d -imi %s -omi %s -ng 3 6 %.5f 4 %.5f 2 %.5f -bgo %d -r %d -nsx %d -nsy %d -rss %d -ko %d" % (template, image, outim, params['il'][0], params['tl'][0], iu, ig, imi, tu, tg, tmi, omi, s1, s2, s3, params['bgo'][0], (params['r'])[opt_ind], (params['nsx'])[spatial_deg], (params['nsy'])[spatial_deg], (params['rss'])[opt_ind], (params['ko'])[spatial_deg]))
            invert(outim)
    else:
        print("\n-> Error: Default must be boolean value\n-> Exiting...")
        os.chdir(cwd)
        sys.exit()
    os.chdir(cwd)

def invert(image):
    data = fits.getdata(image)
    try: mask = fits.getdata(image, 1)
    except: mask = None
    hdr = fits.getheader(image)
    if mask == None:
        hdu = fits.PrimaryHDU(data*-1, header=hdr)
    else:
        hduData = fits.PrimaryHDU(data*-1, header=hdr)
        hduMask = fits.ImageHDU(mask)
        hdu = fits.HDUList([hduData, hduMask])
    hdu.writeto(image, overwrite=True)