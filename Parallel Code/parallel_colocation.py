from subregions_class import Subregion
from border_class import Border

import pandas as pd
import numpy as np
import math
import csv
import ctypes
import os
import geopandas as gpd
import time
import concurrent.futures
from rtree import index
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor

# global variable
border = Border()

def process_subregion(i, subregion, offset, distance_threshold):
    subregion.calc_feature_info(offset)
    subregion.calc_star_neighbors_parallel(distance_threshold)
    subregion.features = subregion.featureInfo.keys()
    
    feature_info = [(feature, values['count'], values['start'], values['end'])
                    for feature, values in subregion.featureInfo.items()]
    
    star_neighbors = [(feature, ' '.join(map(str, neighbors)))
                      for feature, neighbors in subregion.star_neighbors.items()]

    return i, feature_info, star_neighbors

def parallel_process_subregions(subregions, offsets, distance_threshold):
    results = []
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_subregion, i, subregions[i], offsets[i], distance_threshold)
            for i in range(len(subregions))
        ]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    # Sequential file I/O after parallel processing
    for i, feature_info, star_neighbors in results:
        with open(f'IntermediateData/featureInfoParallel/featureInfo{i}.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['feature', 'count', 'start', 'end'])
            writer.writerows(feature_info)

        sorted_star_neighbors = sorted(star_neighbors, key=lambda x: x[0])
        with open(f'IntermediateData/starNeighborsParallel/starNeighbors{i}.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['feature', 'star neighbors'])
            writer.writerows(sorted_star_neighbors)
            
def read_distance_threshold(directory):
    with open(directory, 'r') as file:   # file to read distance threshold path
        distance_threshold = float(file.read())
    return distance_threshold

def read_data(directory_path):
    subregions = []
    dataframes = []
    offsets = [0]
    number_subregions = 0
    for filename in os.listdir(directory_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(directory_path, filename)
            df = pd.read_csv(file_path)
            # apply the offset to the index
            df.index = range(offsets[number_subregions], offsets[number_subregions] + len(df))
            offsets.append(df.shape[number_subregions])
            subregions.append(Subregion(df))
            dataframes.append(df)
            number_subregions += 1
    offsets.pop()
    return subregions, offsets, number_subregions, dataframes

def shapefile_processing(shapefile_path, distance_threshold):
    shapefile = gpd.read_file(shapefile_path)
    points = np.array(border.combined_df[['latitude', 'longitude']])
    points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(points[:, 1], points[:, 0]))
    points_gdf = points_gdf.set_crs("EPSG:4326")
    shapefile = shapefile.to_crs("EPSG:4326")
    
    featureType = []
    x = []
    y = []
    ID = []
    curr_index = 0
    for point in points_gdf.geometry:
        # Convert radius from kilometers to degrees (assuming a spherical Earth)
        radius_deg = distance_threshold / 111.32  # Approximately 111.32 km per degree of latitude
        # Create a circle geometry
        circle = point.buffer(radius_deg)

        # Find the borders that the point intersects
        intersected_borders = [border for border in shapefile['geometry'] if circle.intersects(border)]

        if len(intersected_borders) >= 2:
            featureType.append(border.combined_df['featureType'].iloc[curr_index])
            x.append(border.combined_df['xCoordinate'].iloc[curr_index])
            y.append(border.combined_df['yCoordinate'].iloc[curr_index])
            ID.append(border.combined_df['ID'].iloc[curr_index])
        curr_index += 1
        
    border.border_df = pd.DataFrame({
        'featureType': featureType,
        'xCoordinate': x,
        'yCoordinate': y,
        'ID': ID})
    
    ids = border.border_df['ID'].to_list()
    return ids
    
    border.calc_feature_info()
    with open('IntermediateData/border_featureInfo/featureInfo.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['feature', 'count', 'start', 'end'])
        for feature, values in border.featureInfo.items():
            writer.writerow([feature, values['count'], values['start'], values['end']])
            
    border.calc_star_neighbors(distance_threshold)
    with open('IntermediateData/border_starNeighbors/starNeighbors.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['feature', 'star neighbors'])
        for feature, neighbors in border.star_neighbors.items():
            values_str = ' '.join(map(str, neighbors))
            writer.writerow([feature, values_str])
            
    return

def main():
    distance_directory = 'IntermediateData/distance_threshold_parameter.txt'
    shapefile_path = '/home/amk7r/Parallel-Project/data/north_america/shapefile'
    directory_path = '/home/amk7r/Parallel-Project/data/north_america'
    prevalence_threshold = 0.55    # set the prevalence threshold
    
    distance_threshold = read_distance_threshold(distance_directory)
    subregions, offsets, number_subregions, dataframes = read_data(directory_path)
    
    s = time.time()
    parallel_process_subregions(subregions, offsets, distance_threshold)
    e = time.time()
    print("TIME TAKEN (SEC) STAR NEIGHBORHOOD:", e - s)
    
    # call the cpp functions for subregion processing
    lib = ctypes.CDLL('./c_functions.so')
    s = time.time()
    lib.subregion_main(ctypes.c_int(number_subregions), ctypes.c_double(prevalence_threshold))
    e = time.time()
    print("TIME TAKEN (SEC) SUBREGIONS:", e - s)
    
    border.combined_df = pd.concat(dataframes, ignore_index=True)
    border.combined_df['ID'] = border.combined_df.index

    # sort the df by featureType
    border.combined_df = border.combined_df.sort_values(by='featureType', ignore_index=True)
    ids = shapefile_processing(shapefile_path, distance_threshold)
    
    number_borders = 1
    s = time.time()
    lib.border_main(ctypes.c_int(number_borders), ctypes.c_double(prevalence_threshold))
    e = time.time()
    print("TIME TAKEN (SEC) BORDER:", e - s)
    
    s = time.time()
    '''arr_len = len(ids)
    arr_type = ctypes.c_int * arr_len
    arr_c = arr_type(*ids)
    border_number = 0
    lib.update_border_info.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.c_int)
    lib.update_border_info(arr_c, arr_len, border_number)
    
    lib.combine_hashmaps.argtypes = (ctypes.c_int, ctypes.c_int)
    lib.combine_hashmaps(number_subregions, number_borders)
    lib.combine_instance_tables.argtypes = (ctypes.c_int, ctypes.c_int)
    lib.combine_instance_tables(number_subregions, number_borders)'''
    
    # Define wrapper functions to run in parallel
    def call_combine_hashmaps():
        lib.combine_hashmaps(number_subregions, number_borders)

    def call_combine_instance_tables():
        lib.combine_instance_tables(number_subregions, number_borders)

    # Execute the functions concurrently
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(call_combine_hashmaps), executor.submit(call_combine_instance_tables)]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()  # Retrieve any result if functions return values
                print("Function completed with result:", result)
            except Exception as e:
                print("An error occurred:", e)
    e = time.time()
    print("TIME TAKEN (SEC) PROCESSING:", e - s)
    
    features = list(border.combined_df['featureType'].unique())
    string_ptrs = (ctypes.c_char_p * len(features))()
    string_ptrs[:] = [s.encode() for s in features]
    lib.region_main.argtypes = (ctypes.c_int, ctypes.c_double, ctypes.POINTER(ctypes.c_char_p), ctypes.c_int)
    lib.region_main(number_subregions, prevalence_threshold, string_ptrs, len(features))
    
if __name__ == "__main__":
    main()

    
'''
North America
    Serial:
        sub_region_main: 0.976 sec
        region_main:

    Parallel:
        sub_region_main: 0.581 sec
        region_main:
        
Middle East
    Serial:
        sub_region_main: 325.825 sec
        region_main:

    Parallel:
        sub_region_main: 145.812 sec
        region_main:


South Asia
    Serial:
        sub_region_main:
        region_main:

    Parallel:
        sub_region_main:
        region_main:

'''