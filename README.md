# Parallel Regional Colocation Mining

Authors: Abigail Kelly, Alec Creasy, Yuan Chen

The serial implementation can be found here: https://github.com/abigamkelly/Map-Based-Colocation-Mining.git

This github includes the code for the parallel regional colocation mining framework.  The following files are included:

* **Parallel Code**: this folder contains the code for our parallel regional colocation mining framework.
   * **border_class.py**: python class that contains all the information for the border region
   * **parallel_c_functions.cpp**: c++ functions for the parallel regional colocation mining process
   * **parallel_colocation.py**: python code that calls the c++ functions defines in **parallel_c_functions.cpp**
   * **parallel_distance.py**: python code that calculates the distance threshold
   * **subregions_class.py**: python class that contains all the information for the subregions
   * **IntermediateData**: folder that contains the star neighborhoods and feature information data that is passed between the python and c++ code

### How to configure:

TODO


### How to run:

TODO


### Dataset
