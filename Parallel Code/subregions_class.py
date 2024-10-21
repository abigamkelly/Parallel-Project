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


# This class holds the information and functions pertaining to each subregion
class Subregion:
    def __init__(self, data):
        self.df = pd.DataFrame(data)
        self.df['featureType'] = self.df['featureType'].astype(str)
        self.featureInfo = {}  # Dictionary to store count, start row ID, and end row ID for each feature type
        self.star_neighbors = {}  # Dictionary to store neighbors of different types within distance_threshold for each instance

    # This function will calculate the feature ranges for each type of feature
    def calc_feature_info(self, offset):
        # Initialize variables to track count, start row ID, and end row ID
        count = 0
        start_row_id = 0
        prev_feature = None
        featureInfo = {}

        # Iterate through the DataFrame
        for i, row in self.df.iterrows():
            feature = row['featureType']

            # If feature type changes, update feature_info for the previous feature
            if feature != prev_feature:
                if prev_feature is not None:
                    featureInfo[prev_feature] = {'count': count, 'start': start_row_id, 'end': i - 1}
                count = 1
                start_row_id = i
                prev_feature = feature
            else:
                count += 1

        # Update feature_info for the last feature
        if prev_feature is not None:
            featureInfo[prev_feature] = {'count': count, 'start': start_row_id, 'end': len(self.df) - 1 + offset}

        # Check if the last feature only has one occurrence
        if prev_feature is not None and count == 1:
            featureInfo[prev_feature] = {'count': count, 'start': start_row_id, 'end': start_row_id}
                
        self.featureInfo = featureInfo
        
    def calc_distance(self, point1, point2):
        x1, y1 = point1['xCoordinate'], point1['yCoordinate']
        x2, y2 = point2['xCoordinate'], point2['yCoordinate']
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
           
    def calc_star_neighbors_parallel(self, distance_threshold):
        star_neighbors = {}
        ids = self.df.index.to_numpy()  # Indices of df
        coords = np.array(self.df[['xCoordinate', 'yCoordinate']].values)  # x and y values of the original data
        features = np.array(self.df['featureType'])
        unique_features = np.unique(features)

        # Parallelize over unique features
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(self.process_feature_type, val_to_exclude, features, ids, coords, 
                                       distance_threshold, self.df.copy())  # Copy df for each process
                       for val_to_exclude in unique_features]

            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                star_neighbors.update(result)  # Merge results

        self.star_neighbors = star_neighbors
        
    def process_feature_type(self, val_to_exclude, features, ids, coords, distance_threshold, df):
        star_neighbors = {}
        mask1 = features == val_to_exclude
        mask2 = features != val_to_exclude
        current_feature = coords[mask1]  # Coordinates of current feature type
        current_feature_ids = ids[mask1]  # IDs of current feature type
        exclude_current_feature = coords[mask2]  # Coordinates not including the current feature type
        exclude_current_feature_ids = ids[mask2]  # IDs excluding the current feature type

        idx = index.Index()  # Create a spatial index
        # Insert each point into the spatial index with its unique identifier
        for i, point in enumerate(exclude_current_feature):
            j = exclude_current_feature_ids[i]
            x_coord, y_coord = point[:2]
            idx.insert(j, (x_coord, y_coord, x_coord, y_coord))

        # Iterate over each point of the current feature type
        for i, point in enumerate(current_feature):
            row_id = current_feature_ids[i]
            x_coord, y_coord = point[:2]
            nearby_points = list(idx.intersection((x_coord - distance_threshold, 
                                                   y_coord - distance_threshold, 
                                                   x_coord + distance_threshold, 
                                                   y_coord + distance_threshold)))

            # Filter neighbors based on distance and feature type, and ensure they are greater than the key
            points_to_add = sorted([j for j in nearby_points if j != row_id and 
                                    self.calc_distance(df.loc[j], df.loc[row_id]) <= distance_threshold and 
                                    j > row_id])
            star_neighbors[row_id] = points_to_add

        return star_neighbors

