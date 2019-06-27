Installing **OasisPy**
======================

Getting **OasisPy** is easy.

For those with ``pip``, simply fetch the package with::

	$ pip install OasisPy

Upon installation you will be prompted with a few questions. The first is the path where **OasisPy** will place the **OASIS** directory, which contains the directory tree in which all of **OasisPy**'s image proccessing will take place. You can make this path anything you want, but to make it easier to find later it usually is best to stick with simple location like your **home** directory. The second prompt you will recieve will be asking if you would like to install the ``ISIS`` image subtraction program along with your **OasisPy** installation. The default to this is yes, and should always be so unless you know for a fact you have this software already installed on your machine.

System Requirements
^^^^^^^^^^^^^^^^^^^

**OasisPy** was developed on a Linux machine and must be run on a POSIX-compliant system.

Software Requirements
^^^^^^^^^^^^^^^^^^^^^

There are a number of outside programs **OasisPy** calls on that will need to be installed as well. Most of these are image processessing software written by Emmanuel Bertin and can be fetched from `<https://www.astromatic.net/>`_.

* ``SExtractor`` (*astromatic*)
* ``PSFex`` (*astromatic*)
* ``SWarp`` (*astromatic*)
* ``SkyMaker`` (optional, used for simulations) (*astromatic*)
* ``PyRAF`` (recommended install through Anaconda)
* ``ISIS`` (automatically installed during **OasisPy** setup)

