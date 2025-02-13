
Prerequisites and Installation
=================================

Prerequisites
----------------

To run TROVA, you need:

- `Python 3.8+ <https://www.python.org/downloads/release/python-380/>`__ 
- `Git <https://git-scm.com/>`__ 
- `Anaconda 3 <https://www.anaconda.com/>`__ 
- `Linux <https://www.linux.org/>`__ 
- `Fortran <https://fortran-lang.org/>`__ 


1- Create a python environment with conda, for example:

    conda create -n py38 python=3.8
    conda activate py38

2. The main Python packages that must be installed are the following (consider using the proposed options):

- `numpy` (use `conda install numpy`)
- `mpi4py` (use `pip install mpi4py` or `conda install mpi4py`)
- `time` (use `conda install -c conda-forge time`)
- `netCDF4` (use `conda install -c conda-forge netcdf4`)
- `scipy` (use `conda install scipy`)
- `importlib` (use `conda install -c conda-forge importlib`)
- `cartopy` (use `conda install -c conda-forge cartopy`)
- `setuptools` (use `pip install setuptools==58.2.0`)
- `hdf5` (use `pip install hdf5`)


Installation
------------------


1- First option 
  
You must check that all the packages are installed and that there is no error message when they are imported.

- Clone the repository:

.. code-block:: bash

    git clone https://github.com/tramo-ephyslab/TROVA-master.git


- Enter the TROVA-master/src/ directory and execute the *install_trova.sh* code.

.. code-block:: bash

    sh install_trova.sh

2- Second option:

With this option conda will install the necessary TROVA dependencies.

.. code-block:: bash
 
    conda install -c tramo-ephyslab trova

Once installed to check and create the fortran functions the first time it is used, open an ipython and run the following command:

.. code-block:: python

    import trova

NOTE: From now on it must have been installed in the python environment and can be used directly like any library.

Possible problems with python packages:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- If you have a problem with the mpi4py library, please check that all the necessary executables are in the created environment (e.g. **"libmpi.so.12"**).
If they do not exist, create a symbolic link to the environment you are using (e.g **ln -s /home/jose/WRF/Library/lib/libmpi.so.12 /home/jose/anaconda3/envs/test_env/lib**)

- If there is a problem with netcdf4, do not use *conda install -c conda-forge netcdf4* but the *pip install netcdf4* option.
- If the problem is not resolved, contact the developers: jose.carlos.fernandez.alvarez@uvigo.es or jcfernandez@cesga.es.


