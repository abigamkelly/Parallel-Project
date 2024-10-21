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
* **Data**: this folder contains three real-world data sets along with their shapefiles.  It is recommended that you use the **north_america** data set for the fastest run time.

### How to Configure:
1. Open the files in the **IntermediateData** folder and ensure that all the folders and subfolders are empty
2. Open **parallel_distance.py**
3. If running the code using the **north_america** data set, ensure that line 125 is set to: directory = 'data/north_america'.  If using a different data set, change the data set name.  For example, if using the middle_east data set, directory = 'real_world_data/middle_east'.
4. Open **parallel_colocation.py**
5. If running the code using **north_america** data set, ensure that line 134 is set to: shapefile_path = 'data/north_america/shapefile' and line 135 is set to: directory_path = 'data/north_america'
6. Line 136 is a user-defined prevalence threshold (prevalence_threshold).  This variable can be changed to include more or less prevalent patterns.

### How to Compile and Run:
1. Change your current directory to **Parallel Code**
2. Run **parallel_distance.py** by typing the following command in the terminal: "python parallel_distance.py"
3. Type the following command in the terminal to compile the c++ code: "g++ -O3 -fopenmp -shared -o c_functions.so -fPIC parallel_c_functions.cpp"
4. Run **parallel_colocation.py** by typing the following command in the terminal: "python parallel_colocation.py"
