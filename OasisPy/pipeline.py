#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 17:24:29 2018

@author: andrew
"""

def PIPELINE(path):
    '''The **OASIS Pipeline**. Runs all **OASIS** functions in succession.
    
    :param str path: Path of data file tree (contains the **configs**, **data**, **psf**, **residuals**, **sources**, **templates** directories). Use a comma-separated list for mapping to multiple datasets.
    :returns: All-in-one difference imaging pipeline. Raw science images are placed in the **data** directory, and residual images and source catalogs are outputted into the **residuals** and **sources** directories, respectively.
    
    '''
    import align_astroalign
    from ref_image import ref_image
    import combine_swarp
    import extract
    import check_saturation
    import subtract
    import mask
    import sys
    import psf
    import MR
    import os
    import time
    import initialize
    paths = (path.replace(' ','')).split(',')
    del path
    for path in paths:
        startTime = time.time()
        if os.path.exists(path) == True:
            initialize.create(path)
            location = path + '/data'
#            sat = check_saturation.check_saturate(location)
#            if sat == 0:
            ref_image(location)
            align_astroalign.align2(location)
#            else:
#                check = input("-> Saturated images found\n-> Move saturated images to OASIS archive? (y/n): ")
#                if check == 'y':
#                    check_saturation.move_arch(sat)
#                    ref_image(location)
#                    align_astroalign.align2(location)
#                elif check == 'n':
#                    ref_image(location)
#                    align_astroalign.align2(location)
#                else:
#                    print("-> Error: Unknown input")
#                    sys.exit()
            mask.MASK(path)
            print("-> Combining images using swarp method...")
            psf.PSF(path)
            combine_swarp.swarp(location)
            psf.PSF(path)
            print("\n-> Subtracting images...")
            subtract.SUBTRACT(path)
            print("-> Running SExtractor on residuals...")
            extract.EXTRACT(path, method='indiv')
            MR.MR(path)
            extract.EXTRACT(path, method='MR')
            endTime = time.time()
            print("-> Finished!\n-> Total runtime: %.2f seconds\n" % (endTime-startTime))
        else:
            print("\n-> Error: Unknown path entered\n-> Please enter the path to an existing exposure time directory\n-> Exiting...\n")
            sys.exit()

def pipeline_run_sim(path, sim_x=1, sim_y=1, sim_width=30, sim=True):
    import align_astroalign
    from ref_image import ref_image
    import combine_swarp
    import extract
    import subtract
    import mask
    import sys
    import psf
    import MR
    import os
    import time
    startTime = time.time()
    if os.path.exists(path) == True:
        location = path + '/data'
        ref_image(location)
        align_astroalign.align2(location)
        mask.MASK(path)
        print("-> Combining images using swarp method...")
        psf.PSF(path)
        combine_swarp.swarp(location)
        psf.PSF(path)
        print("\n-> Subtracting images...")
        subtract.SUBTRACT(path)
#        print("\n-> Subtracting images using AIS method...")
#        subtract_ais.isis_sub(path)
#        optimize.perform_optimization(path, s_x=sim_x, s_y=sim_y, s_width=sim_width, simulate=sim, qFloor=0.80, qValue=0.90, qInitial=0.95, use_config_file=False)
        print("-> Running SExtractor on residuals...")
        extract.EXTRACT(path, method='indiv')
        MR.MR(path)
        extract.EXTRACT(path, method='MR')
        endTime = time.time()
        print("-> Finished!\n-> Total runtime: %.2f seconds\n" % (endTime-startTime))
    else:
        print("\n-> Error: Unknown path entered\n-> Please enter the path to an existing exposure time directory\n-> Exiting...\n")
        sys.exit()
        
if __name__ == '__main__':
    path = input("\n-> Enter path to exposure time directory: ")
    PIPELINE(path)