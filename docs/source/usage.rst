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

NOTE: This code is not bug-free. Please report any bugs through 'Issues': https://github.com/tramo-ephyslab/TROVA-master/issues

Input file details
------------------

.. code-block:: bash

   #---------------------
   # Input data features
   #---------------------
 
   # Path to FLEXPART or FLEXPART-WRF partposit binary files [str]
   path_data = "/mnt/lustre/hsm/nlsas/notape/home/uvi/fi/tramo/FLEXPART_DATA/"

   # Path for TROVA outputs [str]
   path_output = "/mnt/lustre/scratch/nlsas/home/uvi/fi/mst/JoseC/TROVA_NEW/output/"

   # Lagrangian tracking mode: ('backward' / 'forward') [str]
   mode = "backward"

   # Atmospheric mass [float]
   mass = 5.148e+18

   # Total number of atmospheric parcels in model simulation [int]
   numP = 1997651

   # Type of file: Set 1 for FELXPART-WRF. Set 2 for FLEXPART older or newer than version 9.  [int]
   type_file = 2

   #--------------------------------------------------------
   # TROVA output domain configuration and simulation period
   #--------------------------------------------------------

   # Spatial resolution for TROVA outputs [float]
   resolution = 1 

   # Number of points in the x-direction for TROVA outputs [int]
   numPdX = 360

   # Number of points in the y-direction for TROVA outputs [int]
   numPdY = 180

   # Lower longitude for TROVA output domain [float]
   x_lower_left = -180

   # Lower latitude for TROVA output domain [float]
   y_lower_left = -90

   # Time step for parcel tracking (minutes) [int]
   dtime = 360

   # Total time for parcel tracking (minutes) [int]
   totaltime = 14400

   # Start date for tracking [int]
   year = 2014
   month = 10
   day = 17
   hour = 00
   min = 00

   # Number of days to perform parcel tracking from start day [int]
   ndays = 1

   #------------------
   # Mask data
   #------------------

   # path to mask file (netcdf)
   file_mask = "/mnt/lustre/scratch/nlsas/home/uvi/fi/mst/JoseC/TROVA_NEW/Masks/CAN.nc"

   # Mask name variable in the mask file [str]
   maskname = "mask"     

   # Latitude variable name in the mask file [str]
   maskvar_lat = "lat"

   # Longitude variable name in the mask file [str]
   maskvar_lon = "lon"

   # Mask value for filtering parcels [int]
   mask_value = 1

   #-----------------------------------
   # Configuration for particle tracking
   #-----------------------------------

   # Subdomain limits for regional models [float]
   # x_left_lower_corner: longitude min, y_left_lower_corner: latitude min, x_right_upper_corner: longitude max, y_right_upper_corner: latitude max
   x_left_lower_corner = -180.0
   y_left_lower_corner = -90.0
   x_right_upper_corner = 180
   y_right_upper_corner = 90.0

   # Model type: ['FLEXPART' / 'FLEXPART-WRF'] [str]
   model = "FLEXPART"

   # Set method = 1 for Stohl and James (2005). Set method = 2 for Sodemann et al. (2008) [int]
   method = 1

   # To filter precipitating parcels ["True" / "False"]  [str]
   filter_parcels_dqdt = False

   # Threshold for filtering precipitating parcels [float]. It is only necessary if filter_parcels_dqdt = True.
   dqdt_threshold = -0.0001

   # To filter parcels by height ["True" / "False"]  [str]
   filter_parcels_height = False

   # Vertical layer for filtering parcels by height [lower_layer, upper_layer] [meters]. It is only necessary if filter_parcels_height = True.
   filter_vertical_layers = [0, 25000]

   # To compute the moisture uptake in vertical layers ["True" / "False"]  [str]
   use_vertical_layers = False
   
   #Vertical layers to compute moisture uptake
   vertical_layers = [0, 750, 1500, 2250, 3000, 4000, 6000, 9000, 12000, 15000, 20000]

   #File output format. Set 1 to activate output format and 0 to deactivate [int]
   output_txt = 0
   output_npy = 0
   output_nc = 1

   #-----------------
   # Other parameters
   #-----------------

   #Target region name [str]
   name_target_region = "CAN"

   #Set file_gz=1 if partposit files are compressed in gz format, else file_gz=0 [int]
   file_gz = 0

   #---------------
   #Auxiliar tools
   #---------------

   #To save particle positions for each time step [str]
   save_position_part = False

   #To save dqdt positions for each dt [str]
   save_position_dqdt = False

   #Plotting identified parcels within the target region at time t0 (year_month_day_hour_min) [True /  False] [str]
   plotting_parcels_t0 = False

   #Ploting identified parcels trajectories on a map [True /  False] [str]
   plotting_parcels_tracks_on_map = False

   #Map limits for plotting [latmin, lonmin, latmax, lonmax, mapcenter, dlat, dlon] [float]
   #map center must be 0 or 180. If center=180, provide lonmin and lonmax in 0-360 format
   maps_limits = [0, -110, 75, 15, 0, 5, 25]

   #Plotting 3D parcels trajectories [True /  False]
   plotting_3Dparcels_tracks = False

   #Calendar leap/noleap ["True" / "False"] [only when certain simulations do not use leap calendar] [str]
   noleap = False

   #Parameter to limit the particles to the domain limits. Consider only in regional models ["True" / "False"]  [str] 
   limit_domain = False
