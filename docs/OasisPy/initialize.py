#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 17:24:29 2018

@author: andrew
"""

def INITIALIZE():
    '''
    Sets up the **OASIS** environment on a new machine. Creates the **OASIS** file tree and installs Alard's ``ISIS`` program (see documentation for details).
    '''
    alert = input("-> Build OASIS file tree in %s? (y/n): " % (loc))
    if alert == 'y':
        initialize(loc)
    elif alert == 'n':
        print("-> Change loc variable in initialize.py to desired OASIS directory path, then run script again")
        sys.exit()
    else:
        print("-> Error: unknown input")
        sys.exit()
    ais_install = input("-> Install ISIS image subtraction software on this machine (requires C shell)? (y/n): ")
    if ais_install == 'y':
        ais_run = os.path.dirname(make_stars.__file__) + '/AIS/package/./install.csh'
        os.system(ais_run)
    elif ais_install == 'n':
        pass
    else:
        print("-> Error: unknown input")
        sys.exit()
