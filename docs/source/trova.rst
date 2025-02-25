API reference
=====================

TROVA module 1 
--------------

.. module:: trova

Below is the description of the functions within the main TROVA module, 
specifically in the *trova.py* script.

.. function:: Kdif_python(matrix1, matrix2, paso)

   Computes the difference between two matrices based on a given step.

   **Parameters** 

   - **matrix1** (*numpy.ndarray*): The first input matrix.

   - **matrix2** (*numpy.ndarray*): The second input matrix.

   - **paso** (*float*): The step value, can be either -1 (backward) or 1 (forward).  

   **Returns**

   - *numpy.ndarray*: The output matrix with computed differences and selected values from the input matrices. 

.. function:: search_row_python(matrix, lista)

   Searches for rows in a matrix that match values in a given list.

   This function takes a matrix and a list of values, and returns a new matrix
   where each row corresponds to a row in the input matrix whose first element
   matches a value in the list. If a value from the list is not found in the matrix,
   the corresponding row in the output matrix is filled with -999.9.

   **Parameters**

   - **matrix** (*numpy.ndarray*): The input matrix to search within.

   - **lista** (*list* or *numpy.ndarray*): The list of values to search for in the first column of the matrix.

   **Returns**

   - *numpy.ndarray*: A new matrix with rows from the input matrix that match the values in the list.
                      Rows that do not match are filled with -999.9.

.. function:: determined_id_python(value_mascara, value_mask)
  
    Determines the indices of elements in a mask that match a given value.

    This function takes a mask array and a value to match, and returns an array
    where each element corresponds to the index of the matching value in the mask.
    If a value from the mask is not found, the corresponding element in the output
    array is filled with -999.

   **Parameters**

   - **value_mascara** (*numpy.ndarray*): The mask array to search within.

   - **value_mask** (*int*): The value to search for in the mask array.

   **Returns**

   - *numpy.ndarray*: An array with indices of the matching values in the mask.
                   Indices that do not match are filled with -999.

.. function:: search_row_python(matrix, lista)

   Searches for rows in a matrix that match values in a given list.

   This function takes a matrix and a list of values, and returns a new matrix
   where each row corresponds to a row in the input matrix whose first element
   matches a value in the list. If a value from the list is not found in the matrix,
   the corresponding row in the output matrix is filled with -999.9.

   **Parameters**

   - **matrix** (*numpy.ndarray*): The input matrix to search within.

   - **lista** (*list* or *numpy.ndarray*): The list of values to search for in the first column of the matrix.

   **Returns**

   - *numpy.ndarray*: A new matrix with rows from the input matrix that match the values in the list.
                      Rows that do not match are filled with -999.9.

.. function:: determined_id_python(value_mascara, value_mask)
  
    Determines the indices of elements in a mask that match a given value.

    This function takes a mask array and a value to match, and returns an array
    where each element corresponds to the index of the matching value in the mask.
    If a value from the mask is not found, the corresponding element in the output
    array is filled with -999.

   **Parameters**

   - **value_mascara** (*numpy.ndarray*): The mask array to search within.

   - **value_mask** (*int*): The value to search for in the mask array.

   **Returns**

   - *numpy.ndarray*: An array with indices of the matching values in the mask.
                   Indices that do not match are filled with -999.

.. function:: check_paths(pfile, path)

   Checks if a given path attribute exists in the provided file object.

   This function attempts to retrieve the value of a specified path attribute
   from a given file object. If the attribute does not exist, it returns an
   empty string.

   **Parameters**

   - **pfile** (*object*): The file object to check for the path attribute.

   - **path** (*str*): The name of the path attribute to retrieve.

   **Returns**

   - *str*: The value of the path attribute if it exists, otherwise an empty string.

.. function:: str2boolean(arg)

   Converts a string representation of truth to a boolean value.

   This function takes a string argument and returns its corresponding boolean value.
   It recognizes several common string representations of true and false values.

   **Parameters**

   - **arg** (*str*): The string to convert to a boolean value. Recognized true values are
                      "yes", "true", "t", "y", "1". Recognized false values are "no", "false",
                      "f", "n", "0". The comparison is case-insensitive.

   **Returns**

   - *bool*: The boolean value corresponding to the input string.

   **Raises**

   - *argparse.ArgumentTypeError*: If the input string does not match any recognized true or false values.

.. function:: ProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='»', printEnd="\r")

   Displays a progress bar in the terminal.

   This function prints a progress bar to the terminal to indicate the progress of a task.
   The progress bar updates with each iteration and shows the percentage of completion.

   **Parameters**

   - **iteration** (*int*): Current iteration (must be between 0 and total).

   - **total** (*int*): Total number of iterations.

   - **prefix** (*str*): Prefix string (optional).

   - **suffix** (*str*): Suffix string (optional).

   - **decimals** (*int*): Positive number of decimals in percent complete (optional).

   - **length** (*int*): Character length of the bar (optional).

   - **fill** (*str*): Bar fill character (optional).

   - **printEnd** (*str*): End character (e.g. "\r", "\r\n") (optional).

   **Returns**

   - *None*

.. function:: get_currentversion()

   Retrieves the current version of the TROVA software.

   This function reads the version information from the VERSION file located
   in the same directory as the script and returns it as a string.

   **Returns**

   - *str*: The current version of the TROVA software.

.. function:: get_lastupdate()

   Retrieves the last update date of the TROVA software.

   This function reads the last update date from the LAST_UPDATE file located
   in the same directory as the script and returns it as a string.

   **Returns**

   - *str*: The last update date of the TROVA software.

.. function:: plotting_tracks_3d(particle_positions, fname)

   Plots 3D tracks of parcels.

   This function creates a 3D plot of parcel tracks using their positions and saves the plot to a file.

   **Parameters**

   - **particle_positions** (*numpy.ndarray*): Array containing the positions of the parcels.

   - **fname** (*str*): The filename to save the plot.

   **Returns**

   - *None*

.. function:: ploting_parcels_tracks_map(particle_positions, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, fname)

   Plots parcel tracks on a map.

   This function creates a 2D map plot of parcel tracks using their positions and saves the plot to a file.

   **Parameters**

   - **particle_positions** (*numpy.ndarray*): Array containing the positions of the parcels.

   - **maps_limits** (*list*): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].

   - **paso** (*int*): Step value indicating the direction of the plot (-1 for backward, 1 for forward).

   - **lat_masked** (*numpy.ndarray*): Array containing the masked latitudes.

   - **lon_masked** (*numpy.ndarray*): Array containing the masked longitudes.

   - **mascara** (*numpy.ndarray*): Array containing the mask values.

   - **value_mask** (*int*): The value to use for the mask.

   - **fname** (*str*): The filename to save the plot.

   **Returns**

   - *None*

.. function:: create_map(maps_limits)

   Creates a map with specified limits.

   This function creates a map with the given latitude and longitude limits and returns the map and its coordinate reference system (CRS).

   **Parameters**

   - **maps_limits** (*list*): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].

   **Returns**

   - *tuple*: A tuple containing the map and its CRS.

.. function:: plotting_parcels_within_target_region(particle_positions, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, fname)

   Plots parcels within the target region on a map.

   This function creates a 2D map plot of parcels within the target region using their positions and saves the plot to a file.

   **Parameters**

   - **particle_positions** (*numpy.ndarray*): Array containing the positions of the parcels.

   - **maps_limits** (*list*): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].

   - **paso** (*int*): Step value indicating the direction of the plot (-1 for backward, 1 for forward).

   - **lat_masked** (*numpy.ndarray*): Array containing the masked latitudes.

   - **lon_masked** (*numpy.ndarray*): Array containing the masked longitudes.

   - **mascara** (*numpy.ndarray*): Array containing the mask values.

   - **value_mask** (*int*): The value to use for the mask.

   - **fname** (*str*): The filename to save the plot.

   **Returns**

   - *None*

.. function:: generate_fecha_simulation(ndias, cyear, cmonth, cday, chours, cminutes)

   Generates a list of simulation dates.

   This function generates a list of dates for the simulation based on the number of days and the initial date and time components provided.

   **Parameters**

   - **ndias** (*int*): Number of days for the simulation.

   - **cyear** (*int* or *list*): Initial year(s) of the simulation.

   - **cmonth** (*int* or *list*): Initial month(s) of the simulation.

   - **cday** (*int* or *list*): Initial day(s) of the simulation.

   - **chours** (*int* or *list*): Initial hour(s) of the simulation.

   - **cminutes** (*int* or *list*): Initial minute(s) of the simulation.

   **Returns**

   - *tuple*: A tuple containing lists of years, months, days, hours, and minutes for the simulation dates.

.. function:: function(latitude, longitude, var, var_layers, use_vlayers, vlayers, method, varpor, filename, path, name_var, unit_var, date_save)

   Creates a NetCDF file with the given data.

   This function creates a NetCDF file with the specified latitude, longitude, variable data, and other attributes.

   **Parameters**

   - **latitude** (*numpy.ndarray*): Array of latitude values.

   - **longitude** (*numpy.ndarray*): Array of longitude values.

   - **var** (*numpy.ndarray*): Array of variable data.

   - **var_layers** (*numpy.ndarray*): Array of variable data for layers.

   - **use_vlayers** (*bool*): Whether to use vertical layers.

   - **vlayers** (*list*): List of vertical layers.

   - **method** (*int*): Method used for processing.

   - **varpor** (*numpy.ndarray*): Array of variable data for sources contribution.

   - **filename** (*str*): Name of the output file.

   - **path** (*str*): Path to save the output file.

   - **name_var** (*str*): Name of the variable.

   - **unit_var** (*str*): Unit of the variable.

   - **date_save** (*numpy.ndarray*): Array of dates for the time dimension.

   **Returns**

   - *None*

.. function:: write_nc(dates, tensor, vartype, filename="output")

   Writes data to a NetCDF file.

   This function writes the given tensor data to a NetCDF file with the specified filename and variable type.

   **Parameters**

   - **dates** (*numpy.ndarray*): Array of dates for the time dimension.

   - **tensor** (*numpy.ndarray*): Tensor data to be written to the NetCDF file.

   - **vartype** (*str*): Type of variable data (e.g., "partpos" or "dqdt").

   - **filename** (*str*): Name of the output file (default is "output").

   **Returns**

   - *None*

.. function:: create_directory(path)

   Creates a directory if it does not exist.

   This function checks if a directory exists at the specified path, and if not, it creates the directory.

   **Parameters**

   - **path** (*str*): The path of the directory to create.

   **Returns**

   - *None*

.. function:: read_binaryFile_fortran(filename, type_file, x_left_lower_corner, y_left_lower_corner, x_right_upper_corner, y_right_upper_corner, limit_domain)

   Reads a binary file using Fortran routines.

   This function reads a binary file based on the specified type and domain limits, and returns the data.

   **Parameters**

   - **filename** (*str*): The name of the binary file to read.

   - **type_file** (*int*): The type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).

   - **x_left_lower_corner** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_left_lower_corner** (*float*): Y-coordinate of the lower left corner of the domain.

   - **x_right_upper_corner** (*float*): X-coordinate of the upper right corner of the domain.

   - **y_right_upper_corner** (*float*): Y-coordinate of the upper right corner of the domain.

   - **limit_domain** (*int*): Whether to limit the domain (1 for yes, 0 for no).

   **Returns**

   - *numpy.ndarray*: The data read from the binary file.

.. function:: load_mask_grid_NR(filename, name_mascara, name_variable_lon, name_variable_lat)

   Loads a mask grid from a NetCDF file.

   This function loads the latitude, longitude, and mask variables from a NetCDF file.

   **Parameters**

   - **filename** (*str*): The name of the NetCDF file to read.

   - **name_mascara** (*str*): The name of the mask variable in the NetCDF file.

   - **name_variable_lon** (*str*): The name of the longitude variable in the NetCDF file.

   - **name_variable_lat** (*str*): The name of the latitude variable in the NetCDF file.

   **Returns**

   - *tuple*: A tuple containing the latitude, longitude, and mask arrays.

.. function:: funtion_interpol_mascara(lat_mascara, lon_mascara, mascara, data)

   Interpolates a mask onto data points.

   This function interpolates the values of a mask onto the given data points using nearest neighbor interpolation.

   **Parameters**

   - **lat_mascara** (*numpy.ndarray*): Array of latitudes for the mask.

   - **lon_mascara** (*numpy.ndarray*): Array of longitudes for the mask.

   - **mascara** (*numpy.ndarray*): Array of mask values.

   - **data** (*numpy.ndarray*): Array of data points to interpolate the mask onto.

   **Returns**

   - *numpy.ndarray*: The interpolated mask values at the data points.

.. function:: plot_point_(lat, lon, mascara)

   Plots points on a map using a mask.

   This function creates a scatter plot of points on a map using the given latitude, longitude, and mask values.

   **Parameters**

   - **lat** (*numpy.ndarray*): Array of latitude values.

   - **lon** (*numpy.ndarray*): Array of longitude values.

   - **mascara** (*numpy.ndarray*): Array of mask values.

   **Returns**

   - *None*

.. function:: funtion_interpol_mascara_2(lat_mascara, lon_mascara, mascara, data)

   Interpolates a mask onto data points.

   This function interpolates the values of a mask onto the given data points using nearest neighbor interpolation.

   **Parameters**

   - **lat_mascara** (*numpy.ndarray*): Array of latitudes for the mask.

   - **lon_mascara** (*numpy.ndarray*): Array of longitudes for the mask.

   - **mascara** (*numpy.ndarray*): Array of mask values.

   - **data** (*numpy.ndarray*): Array of data points to interpolate the mask onto.

   **Returns**

   - *numpy.ndarray*: The interpolated mask values at the data points.

.. function:: determine_id_binary_grid_NR_fortran(data, lat_mascara, lon_mascara, value_mascara, value_mask)

   Determines the indices of elements in a binary grid that match a given value using Fortran routines.

   This function interpolates the mask values onto the data points and determines the indices of elements that match the given value.

   **Parameters**

   - **data** (*numpy.ndarray*): Array of data points.

   - **lat_mascara** (*numpy.ndarray*): Array of latitudes for the mask.

   - **lon_mascara** (*numpy.ndarray*): Array of longitudes for the mask.

   - **value_mascara** (*numpy.ndarray*): Array of mask values.

   - **value_mask** (*int*): The value to search for in the mask array.

   **Returns**

   - *numpy.ndarray*: A submatrix of data points that match the given value.

.. function:: search_row_fortran(lista, matrix)

   Searches for rows in a matrix that match values in a given list using Fortran routines.

   This function takes a matrix and a list of values, and returns a new matrix where each row corresponds to a row in the input matrix whose first element matches a value in the list.

   **Parameters**

   - **lista** (*list* or *numpy.ndarray*): The list of values to search for in the first column of the matrix.

   - **matrix** (*numpy.ndarray*): The input matrix to search within.

   **Returns**

   - *numpy.ndarray*: A new matrix with rows from the input matrix that match the values in the list.

.. function:: calc_A(resolution, lat, lon)

   Calculates the area of grid cells based on latitude and longitude.

   This function calculates the area of each grid cell defined by the given latitude and longitude arrays and the specified resolution.

   **Parameters**

   - **resolution** (*float*): The resolution of the grid cells.

   - **lat** (*numpy.ndarray*): Array of latitude values.

   - **lon** (*numpy.ndarray*): Array of longitude values.

   **Returns**

   - *numpy.ndarray*: An array of the same shape as the input latitude and longitude arrays, containing the area of each grid cell.

.. function:: grid_point(resolution, numPdX, numPdY, x_lower_left, y_lower_left)

   Generates a grid of points based on the specified resolution and domain limits.

   This function generates a grid of latitude and longitude points based on the specified resolution and the coordinates of the lower left corner of the domain.

   **Parameters**

   - **resolution** (*float*): The resolution of the grid cells.

   - **numPdX** (*int*): Number of grid points in the X direction.

   - **numPdY** (*int*): Number of grid points in the Y direction.

   - **x_lower_left** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_lower_left** (*float*): Y-coordinate of the lower left corner of the domain.

   **Returns**

   - *tuple*: A tuple containing two numpy arrays: the latitude and longitude points of the grid.

.. function:: grid_plot_final(lat, lon)

   Generates a grid of points for plotting based on the input latitude and longitude arrays.

   This function generates a grid of latitude and longitude points for plotting, based on the input latitude and longitude arrays.

   **Parameters**

   - **lat** (*numpy.ndarray*): Array of latitude values.

   - **lon** (*numpy.ndarray*): Array of longitude values.

   **Returns**

   - *tuple*: A tuple containing two numpy arrays: the latitude and longitude points for plotting.

.. function:: time_calc(init_time, h_diff)

   Calculates a new time based on the initial time and a time difference in hours.

   This function calculates a new time by adding the specified time difference in hours to the initial time.

   **Parameters**

   - **init_time** (*str*): The initial time in the format "YYYY-MM-DD HH:MM:SS".

   - **h_diff** (*float*): The time difference in hours.

   **Returns**

   - *datetime*: The calculated time.

.. function:: time_calcminutes(init_time, h_diff)

   Calculates a new time based on the initial time and a time difference in minutes.

   This function calculates a new time by adding the specified time difference in minutes to the initial time.

   **Parameters**

   - **init_time** (*str*): The initial time in the format "YYYY-MM-DD HH:MM:SS".

   - **h_diff** (*float*): The time difference in minutes.

   **Returns**

   - *datetime*: The calculated time.

.. function:: generate_file(paso, dtime, totaltime, fecha, path, key_gz, noleap)

   Generates a list of file names and dates for the simulation.

   This function generates a list of file names and corresponding dates for the simulation based on the specified parameters.

   **Parameters**

   - **paso** (*int*): Step value indicating the direction of the simulation (-1 for backward, 1 for forward).

   - **dtime** (*int*): Time step in minutes.

   - **totaltime** (*int*): Total simulation time in minutes.

   - **fecha** (*str*): Initial date and time in the format "YYYY-MM-DD HH:MM:SS".

   - **path** (*str*): Path to save the output files.

   - **key_gz** (*int*): Whether to use gzip compression (1 for yes, 0 for no).

   - **noleap** (*int*): Whether to exclude leap years (1 for yes, 0 for no).

   **Returns**

   - *tuple*: A tuple containing two lists: the list of file names and the list of corresponding dates.

.. function:: read_proccesor(lista_partposi, submatrix, rank, x_left_lower_corner, y_left_lower_corner, x_right_upper_corner, y_right_upper_corner, model, key_gz, type_file, limit_domain)

   Reads and processes binary files in parallel.

   This function reads binary files in parallel using MPI, processes the data, and returns a tensor of the processed data.

   **Parameters**

   - **lista_partposi** (*list*): List of file paths to read.

   - **submatrix** (*numpy.ndarray*): Submatrix of data points to process.

   - **rank** (*int*): Rank of the current MPI process.

   - **x_left_lower_corner** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_left_lower_corner** (*float*): Y-coordinate of the lower left corner of the domain.

   - **x_right_upper_corner** (*float*): X-coordinate of the upper right corner of the domain.

   - **y_right_upper_corner** (*float*): Y-coordinate of the upper right corner of the domain.

   - **model** (*str*): Model type (e.g., "FLEXPART").

   - **key_gz** (*int*): Whether to use gzip compression (1 for yes, 0 for no).

   - **type_file** (*int*): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).

   - **limit_domain** (*int*): Whether to limit the domain (1 for yes, 0 for no).

   **Returns**

   - *numpy.ndarray*: A tensor of the processed data.

.. function:: remove_rows_with_value(tensor, tensor_por, idPart, qIni, ref_index=0, value=-999.9)

   Removes rows from a tensor that contain a specified value.

   This function removes rows from the input tensor, tensor_por, idPart, and qIni arrays where any element in the specified reference index row contains the given value.

   **Parameters**

   - **tensor** (*numpy.ndarray*): The input tensor to filter.

   - **tensor_por** (*numpy.ndarray*): The tensor containing percentage values to filter.

   - **idPart** (*numpy.ndarray*): Array of parcel IDs to filter.

   - **qIni** (*numpy.ndarray*): Array of initial specific humidity values to filter.

   - **ref_index** (*int*): The reference index of the row to check for the specified value (default is 0).

   - **value** (*float*): The value to check for in the reference row (default is -999.9).

   **Returns**

   - *tuple*: A tuple containing the filtered tensor, tensor_por, idPart, and qIni arrays.

.. function:: _backward_dq(lista_partposi, file_mask, name_mascara, name_variable_lon, name_variable_lat, lat_f, lon_f, rank, size, comm, type_file, x_left_lower_corner, y_left_lower_corner, x_right_upper_corner, y_right_upper_corner, model, method, threshold, filter_value, value_mask, key_gz, path_output, use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers, limit_domain, dates)

   Processes backward parcel tracking data.

   This function processes backward parcel tracking data, filters the data based on specified criteria, and returns the results.

   **Parameters**

   - **lista_partposi** (*list*): List of file paths to read.

   - **file_mask** (*str*): Path to the mask file.

   - **name_mascara** (*str*): Name of the mask variable.

   - **name_variable_lon** (*str*): Name of the longitude variable in the mask file.

   - **name_variable_lat** (*str*): Name of the latitude variable in the mask file.

   - **lat_f** (*numpy.ndarray*): Array of latitude values.

   - **lon_f** (*numpy.ndarray*): Array of longitude values.

   - **rank** (*int*): Rank of the current MPI process.

   - **size** (*int*): Total number of MPI processes.

   - **comm** (*MPI.Comm*): MPI communicator.

   - **type_file** (*int*): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).

   - **x_left_lower_corner** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_left_lower_corner** (*float*): Y-coordinate of the lower left corner of the domain.

   - **x_right_upper_corner** (*float*): X-coordinate of the upper right corner of the domain.

   - **y_right_upper_corner** (*float*): Y-coordinate of the upper right corner of the domain.

   - **model** (*str*): Model type (e.g., "FLEXPART").

   - **method** (*int*): Method used for processing.

   - **threshold** (*float*): Threshold value for filtering.

   - **filter_value** (*int*): Value used for filtering.

   - **value_mask** (*int*): Value to search for in the mask array.

   - **key_gz** (*int*): Whether to use gzip compression (1 for yes, 0 for no).

   - **path_output** (*str*): Path to save the output files.

   - **use_vertical_layers** (*bool*): Whether to use vertical layers.

   - **vertical_layers** (*list*): List of vertical layers.

   - **filter_parcels_height** (*bool*): Whether to filter parcels by height.

   - **filter_vertical_layers** (*list*): List of vertical layers for filtering.

   - **limit_domain** (*int*): Whether to limit the domain (1 for yes, 0 for no).

   - **dates** (*list*): List of dates for the simulation.

   **Returns**

   - *tuple*: A tuple containing the processed data and additional information.

.. function:: _forward_dq(lista_partposi, file_mask, name_mascara, name_variable_lon, name_variable_lat, lat_f, lon_f, rank, size, comm, type_file, x_left_lower_corner, y_left_lower_corner, x_right_upper_corner, y_right_upper_corner, model, value_mask, key_gz, path_output, use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers, limit_domain)

   Processes forward parcel tracking data.

   This function processes forward parcel tracking data, filters the data based on specified criteria, and returns the results.

   **Parameters**

   - **lista_partposi** (*list*): List of file paths to read.

   - **file_mask** (*str*): Path to the mask file.

   - **name_mascara** (*str*): Name of the mask variable.

   - **name_variable_lon** (*str*): Name of the longitude variable in the mask file.

   - **name_variable_lat** (*str*): Name of the latitude variable in the mask file.

   - **lat_f** (*numpy.ndarray*): Array of latitude values.

   - **lon_f** (*numpy.ndarray*): Array of longitude values.

   - **rank** (*int*): Rank of the current MPI process.

   - **size** (*int*): Total number of MPI processes.

   - **comm** (*MPI.Comm*): MPI communicator.

   - **type_file** (*int*): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).

   - **x_left_lower_corner** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_left_lower_corner** (*float*): Y-coordinate of the lower left corner of the domain.

   - **x_right_upper_corner** (*float*): X-coordinate of the upper right corner of the domain.

   - **y_right_upper_corner** (*float*): Y-coordinate of the upper right corner of the domain.

   - **model** (*str*): Model type (e.g., "FLEXPART").

   - **value_mask** (*int*): Value to search for in the mask array.

   - **key_gz** (*int*): Whether to use gzip compression (1 for yes, 0 for no).

   - **path_output** (*str*): Path to save the output files.

   - **use_vertical_layers** (*bool*): Whether to use vertical layers.

   - **vertical_layers** (*list*): List of vertical layers.

   - **filter_parcels_height** (*bool*): Whether to filter parcels by height.

   - **filter_vertical_layers** (*list*): List of vertical layers for filtering.

   - **limit_domain** (*int*): Whether to limit the domain (1 for yes, 0 for no).

   **Returns**

   - *tuple*: A tuple containing the processed data and additional information.

.. function:: time_calc_day(init_time, day_diff)

   Calculates a new date based on the initial date and a day difference.

   This function calculates a new date by adding the specified day difference to the initial date.

   **Parameters**

   - **init_time** (*str*): The initial date in the format "YYYY-MM-DD HH:MM:SS".

   - **day_diff** (*int*): The day difference to add.

   **Returns**

   - *datetime*: The calculated date.

.. function:: convert_date_to_ordinal(year, month, day, hour, minute, second)

   Converts a date to an ordinal number.

   This function converts the specified date and time components to an ordinal number based on the NetCDF convention.

   **Parameters**

   - **year** (*int*): The year component of the date.

   - **month** (*int*): The month component of the date.

   - **day** (*int*): The day component of the date.

   - **hour** (*int*): The hour component of the date.

   - **minute** (*int*): The minute component of the date.

   - **second** (*int*): The second component of the date.

   **Returns**

   - *float*: The ordinal number representing the date.

.. function:: decompose_date(value)

   Decomposes an ordinal date value into its components.

   This function decomposes an ordinal date value into its year, month, day, hour, minute, and second components.

   **Parameters**

   - **value** (*float*): The ordinal date value to decompose.

   **Returns**

   - *tuple*: A tuple containing the year, month, day, hour, minute, and second components.

.. class:: InputNotInRangeError

   Exception raised for errors in the input parameters.

   This exception is raised when an input parameter is not within the expected range.

   **Attributes**

   - **message** (*str*): Explanation of the error.

.. function:: to_check_params(paso, type_file, numPdx, numPdY, method, resolution, cant_plazo, file_mask)

   Checks the validity of input parameters.

   This function checks if the input parameters are within the expected range and raises an exception if they are not.

   **Parameters**

   - **paso** (*int*): Step value indicating the direction of the simulation (-1 for backward, 1 for forward).

   - **type_file** (*int*): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).

   - **numPdx** (*int*): Number of grid points in the X direction.

   - **numPdY** (*int*): Number of grid points in the Y direction.

   - **method** (*int*): Method used for processing (1 for Stohl and James, 2 for Sodemann).

   - **resolution** (*float*): The resolution of the grid cells.

   - **cant_plazo** (*int*): Time step in minutes.

   - **file_mask** (*str*): Path to the mask file.

   **Returns**

   - *None*

   **Raises**

   - *InputNotInRangeError*: If any input parameter is not within the expected range.

.. function:: function_proof(lat, lon)

   Checks if latitude and longitude points are within the valid range.

   This function checks if the latitude and longitude points are within the valid range and raises an exception if they are not.

   **Parameters**

   - **lat** (*numpy.ndarray*): Array of latitude values.

   - **lon** (*numpy.ndarray*): Array of longitude values.

   **Returns**

   - *None*

   **Raises**

   - *InputNotInRangeError*: If any latitude or longitude point is not within the valid range.

.. function:: desc_gz(name_file)

   Decompresses a gzip file.

   This function decompresses the specified gzip file and saves the decompressed content to a new file.

   **Parameters**

   - **name_file** (*str*): The name of the gzip file to decompress.

   **Returns**

   - *None*

.. function:: TROVA_LOGO()

   Prints the TROVA logo.

   This function prints the TROVA logo to the terminal.

   **Returns**

   - *None*

.. function:: main_process(path, paso, comm, size, rank, resolution, numPdX, numPdY, dtime, totaltime, year, month, day, hour, minn, time, path_output, file_mask, name_mascara, name_variable_lon, name_variable_lat, x_lower_left, y_lower_left, type_file, masa, numP, x_left_lower_corner, y_left_lower_corner, x_right_upper_corner, y_right_upper_corner, model, method, threshold, filter_value, output_txt, output_npy, output_nc, value_mask, key_gz, save_position_part, use_vertical_layers, vertical_layers, save_position_dqdt, filter_parcels_height, filter_vertical_layers, plotting_parcels_t0, plotting_parcels_tracks_on_map, plotting_3Dparcels_tracks, maps_limits, noleap, limit_domain)

   Main processing function for TROVA.

   This function performs the main processing tasks for TROVA, including reading input files, processing data,
   and generating output files and plots.

   **Parameters**

   - **path** (*str*): Path to the input data.

   - **paso** (*int*): Step value indicating the direction of the simulation (-1 for backward, 1 for forward).

   - **comm** (*MPI.Comm*): MPI communicator.

   - **size** (*int*): Total number of MPI processes.

   - **rank** (*int*): Rank of the current MPI process.

   - **resolution** (*float*): The resolution of the grid cells.

   - **numPdX** (*int*): Number of grid points in the X direction.

   - **numPdY** (*int*): Number of grid points in the Y direction.

   - **dtime** (*int*): Time step in minutes.

   - **totaltime** (*int*): Total simulation time in minutes.

   - **year** (*str*): Initial year of the simulation.

   - **month** (*str*): Initial month of the simulation.

   - **day** (*str*): Initial day of the simulation.

   - **hour** (*str*): Initial hour of the simulation.

   - **minn** (*str*): Initial minute of the simulation.

   - **time** (*float*): Start time of the simulation.

   - **path_output** (*str*): Path to save the output files.

   - **file_mask** (*str*): Path to the mask file.

   - **name_mascara** (*str*): Name of the mask variable.

   - **name_variable_lon** (*str*): Name of the longitude variable in the mask file.

   - **name_variable_lat** (*str*): Name of the latitude variable in the mask file.

   - **x_lower_left** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_lower_left** (*float*): Y-coordinate of the lower left corner of the domain.

   - **type_file** (*int*): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).

   - **masa** (*float*): Mass of the parcels.

   - **numP** (*int*): Number of parcels.

   - **x_left_lower_corner** (*float*): X-coordinate of the lower left corner of the domain.

   - **y_left_lower_corner** (*float*): Y-coordinate of the lower left corner of the domain.

   - **x_right_upper_corner** (*float*): X-coordinate of the upper right corner of the domain.

   - **y_right_upper_corner** (*float*): Y-coordinate of the upper right corner of the domain.

   - **model** (*str*): Model type (e.g., "FLEXPART").

   - **method** (*int*): Method used for processing.

   - **threshold** (*float*): Threshold value for filtering.

   - **filter_value** (*int*): Value used for filtering.

   - **output_txt** (*int*): Whether to output results in TXT format (1 for yes, 0 for no).

   - **output_npy** (*int*): Whether to output results in NPY format (1 for yes, 0 for no).

   - **output_nc** (*int*): Whether to output results in NetCDF format (1 for yes, 0 for no).

   - **value_mask** (*int*): Value to search for in the mask array.

   - **key_gz** (*int*): Whether to use gzip compression (1 for yes, 0 for no).

   - **save_position_part** (*bool*): Whether to save parcel positions at each time step.

   - **use_vertical_layers** (*bool*): Whether to use vertical layers.

   - **vertical_layers** (*list*): List of vertical layers.

   - **save_position_dqdt** (*bool*): Whether to save dq/dt at each time step.

   - **filter_parcels_height** (*bool*): Whether to filter parcels by height.

   - **filter_vertical_layers** (*list*): List of vertical layers for filtering.

   - **plotting_parcels_t0** (*bool*): Whether to plot identified parcels within the target region at time t0.

   - **plotting_parcels_tracks_on_map** (*bool*): Whether to plot identified parcels' trajectories on a map.

   - **plotting_3Dparcels_tracks** (*bool*): Whether to plot 3D parcels' trajectories.

   - **maps_limits** (*list*): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].

   - **noleap** (*bool*): Whether to exclude leap years.

   - **limit_domain** (*int*): Whether to limit the domain (1 for yes, 0 for no).

   **Returns**

   - *None*

.. function:: TROVA_main(input_file)

   Main function for TROVA.

   This function reads the input configuration file, initializes parameters, and starts the main processing function.

   **Parameters**

   - **input_file** (*str*): Path to the input configuration file.

   **Returns**

   - *None*

TROVA module 2
--------------

.. module:: tensor_operations

Below is the description of the functions within the TROVA complementary module, 
specifically in the *tensor_operations.f90* script.