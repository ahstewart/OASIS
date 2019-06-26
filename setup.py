from distutils.core import setup
from setuptools.command.install import install
import os
import sys

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        loc = input("\n-> Location to build OASIS file tree? (default=cwd) ")
        if loc == '':
            loc = os.getcwd()
        if os.path.exists(loc) == True:
            if os.path.exists(loc + "/OASIS") == False:
                os.system("mkdir %s/OASIS" % (loc))
                os.system("mkdir %s/OASIS/targets" % (loc))
                os.system("mkdir %s/OASIS/temp" % (loc))
                os.system("mkdir %s/OASIS/archive" % (loc))
                os.system("mkdir %s/OASIS/config" % (loc))
                os.system("mkdir %s/OASIS/simulations" % (loc))
                os.system("mkdir %s/OASIS/archive/data" % (loc))
                os.system("mkdir %s/OASIS/archive/templates" % (loc))
                os.system("mkdir %s/OASIS/archive/residuals" % (loc))
                print("-> OASIS file system created in %s\n" % (loc))
            else:
                print("-> OASIS architecure already exists on this computer")
            ais_run = os.path.dirname(__file__) + '/OASIS/AIS/package/./install.csh'
            os.system(ais_run)
            with open(os.path.dirname(__file__) + '/OASIS/config/OASIS.config', 'a') as conf:
                conf.write("loc \t %s \t # location of OASIS file tree. DO NOT CHANGE." % loc)
            os.system("cp %s/config/OASIS.config %s/OASIS/config" % (os.path.dirname(make_stars.__file__), loc))
        else:
            print("-> Error: Location does not exist\n-> Exiting...\n")
            sys.exit()
        install.run(self)

with open("README", "r") as fh:
    long_description = fh.read()

setup(name='OasisPy',
      version='1.0',
      description='Difference Imaging Engine for Optical SETI Applications',
      long_description=long_description,
      author='Andrew Stewart',
      url='https://github.com/andrewhstewart/OasisPy.git',
      author_email='ah.stewart@outlook.com',
      packages=['OasisPy'],
      package_dir={'OasisPy': 'OasisPy'},
      package_data={'OasisPy': ['AIS/package/*', 'AIS/package/abs/*', 'AIS/package/bin/*', 'AIS/package/Bphot/*', 'AIS/package/Cphot/*', 'AIS/package/cross/*', 'AIS/package/czerny/*', 'AIS/package/detect/*', 'AIS/package/extract/*', 'AIS/package/fit2d/*', 'AIS/package/images/*', 'AIS/package/interp/*', 'AIS/package/images2/*', 'AIS/package/phot_ref/*', 'AIS/package/register/*', 'AIS/package/stack/*', 'AIS/package/subtract/*', 'AIS/package/utils/*', 'config/*', 'test_config/*']},
      include_package_data = True,
      cmdclass={'install': PostInstallCommand},
      install_requires=[
          'astropy',
          'numpy',
          'glob',
          'copy',
          'sys',
          'datetime',
          'os',
          'shutil',
          'tabulate',
          'scipy',
          'requests',
          'time',
          'MontagePy',
          'tqdm',
          'math',
          'astroscrappy',
          'astroalign',
          'image_registration',
          'operator',
          'calendar',
          'zipfile',
          'gzip',
          'photutils']
      )