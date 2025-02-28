What’s New
==========

Releases
--------

v.1.1.1 (2025-02-15)
~~~~~~~~~~~~~~~~~~~~

- Developed in Python, responsible for reading files, configuring TROVA, and generating 
  the outputs of the moisture balance (Evaporation (E)-Precipitation (P)) for the number 
  of days selected in the simulations.

- Developed in Fortran, used in interface with Python to perform computationally 
  demanding calculations in the shortest possible time. It also includes a parallel 
  implementation using the MPI library to reduce TROVA's processing time.

- This new version includes the analysis of moisture sources and sinks by vertical layers.

- This version allows the calculation of the residence time of water vapor in the 
  atmosphere for particles in a target region applying the methodology of Läderach and Sodemann (2016).

- This version has functions that allow the representation of moisture source and sink patterns and 
  the representation in a 2D graph of the residence time values of water vapor in the atmosphere 
  for particles in a target region.

v.1.1 (2023-09-16)
~~~~~~~~~~~~~~~~~~

- Developed in Python, responsible for reading files, configuring TROVA, and generating 
  the outputs of the moisture balance (Evaporation (E)-Precipitation (P)) for the number 
  of days selected in the simulations.

- Developed in Fortran, used in interface with Python to perform computationally 
  demanding calculations in the shortest possible time. It also includes a parallel 
  implementation using the MPI library to reduce TROVA's processing time.

- This new version includes the analysis of moisture sources and sinks by vertical layers.

v.1.0 (2022-12-01)
~~~~~~~~~~~~~~~~~

- Developed in Python, responsible for reading files, configuring TROVA, and generating 
  the outputs of the moisture balance (Evaporation (E)-Precipitation (P)) for the number 
  of days selected in the simulations.

- Developed in Fortran, used in interface with Python to perform computationally 
  demanding calculations in the shortest possible time. It also includes a parallel 
  implementation using the MPI library to reduce TROVA's processing time.