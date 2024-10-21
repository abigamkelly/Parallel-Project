# Parallel Regional Colocation Mining

Authors: Abigail Kelly, Alec Creasy, Yuan Chen

Original Code: https://github.com/abigamkelly/Map-Based-Colocation-Mining.git

This github includes the code for the parallel regional colocation mining framework.  The following files are included:

* **Parallel Code**: this folder contains the code for our parallel regional colocation mining framework.
   * **border_class.py**: python class that contains all the information for the border region
   * **parallel_c_functions.cpp**: c++ functions for the parallel regional colocation mining process
   * **parallel_colocation.py**: python code that calls the c++ functions defines in **parallel_c_functions.cpp**
   * **parallel_distance.py**: python code that calculates the distance threshold
   * **subregions_class.py**: python class that contains all the information for the subregions
   * **IntermediateData**: folder that contains the star neighborhoods and feature information data that is passed between the python and c++ code
* **Serial Code**: this folder contains the code for our serial regional colocation mining framework.
   * **regional_colocation_compression**: this folder contains the code for our map-based regional colocation mining framework.
        * **c_functions.cpp**: c++ functions used in the regional colocation framework
        * **distance_threshold_calcalculation.ipynb**: python code that estimates the optimal spatial neighborhood relationship constraint
        * **regional_colocation.ipynb**: python code that calls the c++ code in **c_functions.cpp** to perform the colocation mining
        * **required_files**: this folder holds the intermediate data produced by the framework
        * **real_world_data**: this folder contains 3 real-world data sets along with their shapefiles.  It is recommended that you use the **NorthAmerica** data set if you are to run the code due to its shorter run time.
* **colocation_compression**: this folder contains the code for our map-based colocation code (same code as above minus the regional part- this code is only useful for testing synthetic data sets)
        * **c_functions.cpp**: c++ functions called in **compression.ipynb**
        * **compression.ipynb**: python code that calls the c++ code in **c_functions.cpp** to perform the map-based colocation mining
        * **required_files**: intermediate data produced by the code
        * **synthetic_data**: this folder contains 5 synthetic data sets varying in clumpiness.  They are titled TestCase1_#.csv where # represents the clumpiness.  For example, **TestCase1_1.csv** is clumpiness 1.

### How to configure:

TODO


### How to run:

TODO


### Dataset
