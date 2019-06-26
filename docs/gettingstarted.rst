Getting Started
===============


**OasisPy** Basics
------------------

**OasisPy** is a collection of python modules that facilitate the OASIS difference imaging process. Each module is in charge of one step in the process. There are a total of 13 modules, and each can either be called individually or in succession, depending on the project at hand. We call these modules "\ **OASIS** methods."

* ``initialize``- Sets up the **OASIS** environment. Is automatically run during install and will likely never have to be run again unless you want to duplicate or move the **OASIS** file tree.
* ``get``- Downloads images from an online archive or finds images on the local hardrive and moves them into the **OASIS** file tree to be processed.
* ``mask``- Masks cosmic rays, saturated stars, CCD defects, and any other artifacts that will inhibit the difference imaging process.
* ``psf``- Computes a PSF model for each image in the data set. This model will be used in many of the subsequent steps.
* ``align``- Chooses the highest S/N frame to be the "reference image," then registers all other images to this reference pointing (to subpixel precision).
* ``combine``- Stacks images into a high S/N template frame to be used for subtraction.
* ``subtract``- Subtracts the template from each image in the data set to create a set of residual frames. This is the most complex and computationally expensive step, and will take the longest to execute.
* ``mr``- Stands for **m**\ aster **r**\ esidual. The residual images created with ``subtract`` are stacked to for a master residual frame.
* ``extract``- Searches the residual images for sources of significant flux. Filters out non point source-like objects and other false positives. Outputs a complete source catalog with all detected variable objects, their position, and their total flux.
* ``test``- Tests the installation of **OasisPy** by fetching a set of images of exoplanet WASP-50b and running them through the **OASIS Pipeline**.
* ``simulation``- Allows further testing of the **OASIS Pipeline** through two different simulated data sets.
* ``run``- Master method that facilitates access to all other methods.
* ``pipeline``- The **OASIS Pipeline**, an all-in-one method that executes the entire difference imaging process, from ``get`` to ``extract``.

For more information on each method and how OasisPy works see :doc:`howitworks`.


Using **OasisPy**
-----------------

**OASIS** methods can be executed in two ways.

As A Python Module
^^^^^^^^^^^^^^^^^^

The simplest way to use **OasisPy** is to treat each method like a simple python module and import them into directly into your code. If comfortable with scripting this is the recommended method of operation as it provides the user with the maximum amount of control over the difference imaging process.

Methods can be called within your python script with the following import::

	from OasisPy import methodname

See :doc:`API` for the details on each method.


Shell Execution
^^^^^^^^^^^^^^^

Methods can also be run from the shell (ANSI) using the formula

.. code-block:: bash

	$ oasis-methodname

For example, to align a set of images one would type

.. code-block:: bash

	$ oasis-align

Each method has a certain number of input parameters the user must provide. Often this just consists of the location of the images being differenced. These parameters are provided by the user through input prompts executed after the initial calling of the method.

This shell execution of **OasisPy** was included to allow those not comfortable with scripting in python, specifcally undergradute or high school students with little or no coding experience, a relatively easy way to access the package's main functionalities.


Convenience Methods
-------------------

Understanding what each of the 13 methods does is critical for complex projects and will allow you to get the most out of the software, but it is not mandatory. For many projects, all you may want to do is call a high level function to do all of the difference imaging steps for you. For this reason we have included several convenience functions for those that do not want to deal with the inner workings of **OASIS**.

Run
^^^
This is the most important convenience function to know. It can be thought of as the "main page" of **OASIS**, from which you can execute any other **OASIS** method. To call ``run`` simply type

.. code-block:: bash

	$ oasis-run

in your terminal. A list of all possible **OASIS** methods is displayed, and from this list you can pick a command and type it into the prompt. It is reccommended to users new to linux and python that ``run`` be the only method directly used.

Pipeline
^^^^^^^^
The ``pipeline`` convenience function makes up what is called the '\ **OASIS Pipeline** \.' This is simply a conglomerate of every **OASIS** method into a single master method. Input data are fed into each method one-by-one and then piped to the next. Using ``pipeline``, a user can send a set of images through the entire difference imaging process with a single high-level command, without worrying about what is actually being done in the intermediary steps. To execute it, type ``oasis-pipeline`` in the terminal or select 'pipeline' if using ``run``.
