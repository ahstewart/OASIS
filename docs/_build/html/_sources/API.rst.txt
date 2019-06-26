API
===


**OASIS** Methods
-----------------

.. autofunction:: OasisPy.initialize.INITIALIZE()

.. autofunction:: OasisPy.get.GET()

.. autofunction:: OasisPy.mask.MASK(path)

.. autofunction:: OasisPy.align.ALIGN(path, align_method='standard')

.. autofunction:: OasisPy.psf.PSF(path)

.. autofunction:: OasisPy.combine.COMBINE(path)

.. autofunction:: OasisPy.subtract.SUBTRACT(path, method='ois', use_config_file=True)

.. autofunction:: OasisPy.MR.MR(path, method='swarp', sig_thresh=4, gauss_sig=3, gauss_filt=False, use_config_file=True)

.. autofunction:: OasisPy.extract.EXTRACT(path, method='both')


Convenience Functions
---------------------

.. autofunction:: OasisPy.pipeline.PIPELINE(path)

.. autofunction:: OasisPy.run.RUN()


Auxillary Functions
-------------------


.. autofunction:: OasisPy.simulation.sim_fakes(location, n_fakes, iterations, input_mode='flux', PSF='moffat', subtract_method='ois', f_min=0, f_max=40000)

.. autofunction:: OasisPy.simulation.sim_sameField(location, numIms=100, bkg_mag=22.5, fwhm_min=3, fwhm_max=6, rot_min=-2.5, rot_max=2.5, shift_min=-2, shift_max=2, scale_mult=(0,1.5), scale_add=(-20,50), zero_point=25, mode='gauss')

.. autofunction:: OasisPy.simulation.SIM()

.. autofunction:: OasisPy.test.TEST()

.. autofunction:: OasisPy.montage.MOSAIC()
