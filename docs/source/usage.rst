TROVA running
=================================

How do I run TROVA?
----------------------

To run TROVA version 1.1, create a file with the following code (could be **run_TROVA.py**)

.. code-block:: python

   #!/usr/bin/env python3
   import trova as tv
   import sys
   input_file=sys.argv[1]
   tv.TROVA_main(input_file)

Once this code is created, TROVA can be used in the following way:

1- On a Linux computer:

.. code-block:: bash

   mpirun -np num_CPU python TROVA.py input_file_path

   e.g:  mpirun -np 4 python TROVA.py input_back_TC.cfg

2- On a HPC with Linux (See  https://github.com/tramo-ephyslab/TROVA-master/tree/main/run_example_HPC):

.. code-block:: bash

   #!/bin/bash -l
   #SBATCH --mem=120GB
   #SBATCH -N 1
   #SBATCH -n 4
   #SBATCH -t 00:10:00

   module --purge
   module load cesga/2020
   module load miniconda3/4.9.2
   conda activate envname

   srun -n $SLURM_NTASKS --mpi=pmi2 python  run_TROVA.py input_file_path >> py_${SLURM_JOB_ID}.log


- num_CPU: CPU numbers to use (preferably divisible by 4).
- input_file_path: Input file path for TROVA run.

Below is an image of what a run with TROVA would look like, from when it initializes 
showing information about the input data, mask, etc. used, such as the time it takes 
to process the number of files or particles.


.. image:: _static/run_trova.png
   :alt: run trova
   :align: center
   :width: 400px

NOTE: This code is not bug-free. Please report any bugs through 'Issues': https://github.com/tramo-ephyslab/TROVA-master/issues

Input file details
------------------

Below are the parameters of the TROVA configuration file. This is just an example, 
so the user should modify it based on the input data.

.. code-block:: bash

   #*****************************************************************************************
   #*                    EPhysLab (Environmental Physics Laboratory), Spain                 *
   #*                        Galician Supercomputing Center, Spain                          *
   #*                        TRansport Of water VApor (TROVA)                               *
   #*                             version 1.1.1 (15-02-2025)                                *
   #*                        _____ __    ____                                               *
   #*                          |  |  |  /    \ \        //\                                 *
   #*                          |  |__| /      \ \      //__\                                *
   #*                          |  |  \ \      /  \    //    \                               *
   #*                          |  |   \ \____/    \__//      \                              *
   #*                                                                                       *
   #*                       Edificio Campus da Auga/Edificio CESGA                          *
   #*                            University of Vigo/CESGA                                   *
   #*                          www.ephyslab.uvigo.es/www.cesga.es                           *
   #*      contact: jose.carlos.fernandez.alvarez@uvigo.es (jcfernandez@cesga.es),          * 
   #*                         albenis.perez.alarcon@uvigo.es                                *
   #*****************************************************************************************
   #------------------------------------------------------------------------------------------
   #Path to FLEXPART or FLEXPART-WRF partposit binary files [str]
   path_data = "/home/jose/Documentos/TROVArun/Data/"

   #Path for TROVA outputs [str]
   path_output = "output/"

   #Lagrangian tracking mode: ('backward' / 'forward'/ 'wvrt' / 'partposit') [str]
   #backward: moisture sources, forward: moisture sinks, wvrt: water vapor residence time, partposit: particle variables over target region
   mode = "backward"

   #Atmospheric mass [float]
   mass = 1.165725e+18

   #Total number of atmospheric parcels in model simulation [int]
   numP = 2045128

   #Type of file: Set 1 for FELXPART-WRF and FLEXPART newler than version 9. Set 2 for FLEXPART older than version 9.  [int]
   type_file = 1

   #Spatial resolution for TROVA outputs [float]
   resolution = 0.25 

   #Number of point in x-direction for TROVA outputs [int]
   numPdX = 600

   #Number of point in y-direction for TROVA outputs [int]
   numPdY = 325

   #Lower longitude for TROVA output domain [float]
   x_lower_left = -110

   #Lower latitude for TROVA output domain [float]
   y_lower_left = -15

   #Time step for parcel tracking (minutes) [int]
   dtime = 360

   #Total time for parcel tracking (minutes) [int]
   totaltime = 1440

   #Start date for tracking [int]
   year = 2014
   month = 10
   day = 17
   hour = 00
   min = 00

   #Number of days to perform parcel tracking from start day [int]
   ndays = 1

   #path to mask fil (netcdf)
   file_mask = "Masks/CAN.nc"

   #Mask name variable in the mask filee [str]
   maskname = "mask"     

   #Latitude variable name  in the mask file [str]
   maskvar_lat = "lat"

   #Longitude variable name in the mask file [str]
   maskvar_lon = "lon"

   #Mask value for filterirng parcels [int]
   mask_value = 1

   #Subdomain limits for regional models [float]
   #x_left_lower_corner: longitude min, y_left_lower_corner: latitude min, x_right_upper_corner: longitude max, y_right_upper_corner: latitude max
   x_left_lower_corner = -100.0
   y_left_lower_corner = -15.0
   x_right_upper_corner = 39.86
   y_right_upper_corner = 57.0

   #model type: ['FLEXPART' / 'FLEXPART-WRF'] [str]
   model = "FLEXPART-WRF"

   #Set method = 1 for Stohl and James (2004,2005). Set method = 2 for Sodemann et al. (2008) [int]
   method = 1

   #To filter precipitating parcels ["True" / "False"]  [str]
   filter_parcels_dqdt = False

   #Threshold for filtering precipitating parcels [float]. It is only necessary if filter_parcels_dqdt = True.
   dqdt_threshold = -0.0001

   #To filter parcels by heigh ["True" / "False"]  [str]
   filter_parcels_height = False

   #Vertical layer for filtering parcels by height [lower_layer, upper_layer] [meters]. It is only necessary if filter_parcels_height = True.
   filter_vertical_layers = [0,25000]

   #To compute the moisture uptake in vertical layers ["True" / "False"]  [str]
   use_vertical_layers = False

   #Vertical layers to compute moisture uptake
   vertical_layers = [0, 750, 900, 1500, 2250, 3000, 4000, 6000, 9000, 12000, 15000, 20000]

   #File output format. Set 1 to activate output format and 0 to deactivate [int]
   output_txt = 0
   output_npy = 0
   output_nc = 1

   #Target region name [str]
   name_target_region = "CAN"

   #Set file_gz=1 if partposit files are compressed in gz format, else file_gz=0 [int]
   file_gz = 0

   #To save particle positions for each time step [str]
   save_position_part = True

   #To save dqdt positions for each dt [str]
   save_position_dqdt = True

   #Plotting identified parcels within the target region at time t0 (year_month_day_hour_min) [True /  False] [str]
   plotting_parcels_t0 = False

   #Ploting identified parcels trajectories on a map [True /  False] [str]
   plotting_parcels_tracks_on_map = False

   #Map limits for plotting [latmin, lonmin, latmax, lonmax, mapcenter, dlat, dlon] [float]
   #map center must be 0 or 180. If center=180, provide lonmin and lonmax in 0-360 format
   maps_limits = [0, -110, 75, 15, 0, 5, 25]

   #Plotting 3D parcels trajectories [True /  False]
   plotting_3Dparcels_tracks = False

   #Calendar leap/noleap ["True" / "False"]  [str]
   noleap = False

   #Parameter to limit the particles to the domain limits. Consider only in regional models ["True" / "False"]  [str]
   limit_domain = True

   #Parameter to activate the calculation of the water vapor residence time together with the calculation of the humidity sources in backward mode ["True" / "False"]  [str]
   method_wvrt = False

   #Plotting the pattern of moisture sources or sinks [True /  False] [str]
   plotting_moisture_sink_source = True

   #Color pallete limits for plotting [min, max, step] [float]
   limits_plot = [0, 3.2, 0.2]


