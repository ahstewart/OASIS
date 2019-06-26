.. OasisPy documentation master file, created by
   sphinx-quickstart on Wed May  8 12:21:09 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**O**\ptic\ **A**\l **S**\eti **I**\mage **S**\ubtraction-- ``OasisPy``
==========================================================================

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   install
   gettingstarted
   howitworks
   tutorial
   troubleshooting
   API


**OASIS** is a toolkit for detecting variable objects in astronomical images by means of difference imaging. Includes the **OASIS** Pipeline, an all-in-one difference imaging utility that takes a set of input images and performs all necessary difference imaging steps, outputting a set of sources catalogs upon completion. Difference imaging is a notoriously cumbersome task, especially for widely varying data. The **OASIS Pipeline** was built as a way to largely automate many of the menial processing steps involved in a difference imaging project.

The code is designed to perform quality difference imaging on data that vary widely in pointing, background, seeing, etc. Originally used in processing images of large galaxies, OASIS should work well for both extended objects and simple star fields.

It was developed for use in UC Santa Barbara's Optical SETI program (`project homepage <https://www.deepspace.ucsb.edu/projects/implications-of-directed-energy-for-seti>`_), but can be deployed in any application involving anomaly detection in astronomical data (see section **something**). 

The **OasisPy** package is a set of Python modules that facilitate access to **OASIS**'s main functionalities.


Features
--------
* **Masking** -- Masks cosmic rays, hot pixels, saturated stars, CCD defects, etc. Supports the use of weight maps often used in AstrOmatic programs.
* **PSF Modeling** -- Computes PSF models of all input images using the AstrOmatic software `PSFex`.
* **Quality Control** -- Ignores images below a user-defined *S/N* threshold and/or above a seeing threshold.
* **Registration** -- Registers images to a chosen reference frame to subpixel precision. 
* **Photometric Alignment** -- Linearly rescales each image image's intensity scale to match that of the reference image. 
* **Stacking** -- Performs a weighted coaddition of the input images to construct a deep, high *S/N* template image for use in the image subtraction step.
* **Background Matching** -- Matches the background of the input images to the template image, using an image subtraction method that works well for extended objects with complicated backgrounds.
* **Image Subtraction** -- Computes a PSF-matching convolution kernel to convolve with the template image, then subtractsthe template from the input image. Uses the Optimal Image Subtraction (OIS) algorithm from Christophe Alard (alard paper link).
* **Parameter Optimization** -- Iterates over a range of OIS parameter configurations looking for the one that yields the best residual image. If a residual fails to meet a certain quality threshold for all parameter configurations, it is masked. Allows for a more robust subtraction that guarantees all residuals in the dataset will be of optimal quality.
* **Source Extraction** -- Uses the AstrOmatic program `SExtractor` to extract variable objects from residual frames.
* **Source Filtering** -- Filters out subtraction artifacts and other phony variable sources.

