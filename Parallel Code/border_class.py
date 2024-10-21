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


# This class holds the information pertaining to the border region
class Border:
    def __init__(self):
        self.combined_df = pd.DataFrame()
        self.border_df = pd.DataFrame()
        self.featureInfo = {}
        self.star_neighbors = {}
        
    # This function will calculate the feature ranges for each type of feature
    def calc_feature_info(self):
        # Initialize variables to track count, start row ID, and end row ID
        count = 0
        start_row_id = 0
        prev_feature = None
        featureInfo = {}

        # Iterate through the DataFrame
        for i, row in self.border_df.iterrows():
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
            featureInfo[prev_feature] = {'count': count, 'start': start_row_id, 'end': len(self.border_df) - 1}

        # Check if the last feature only has one occurrence
        if prev_feature is not None and count == 1:
            featureInfo[prev_feature] = {'count': count, 'start': start_row_id, 'end': start_row_id}
                
        self.featureInfo = featureInfo
        
    def calc_distance(self, point1, point2):
        x1, y1 = point1['xCoordinate'], point1['yCoordinate']
        x2, y2 = point2['xCoordinate'], point2['yCoordinate']
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
    # This function will calculate the star neighbors for each instance
    def calc_star_neighbors(self, distance_threshold):
        star_neighbors = {}
        # Create a spatial index
        idx = index.Index()
        # Insert each point into the spatial index with its unique identifier
        for i, row in self.border_df.iterrows():
            x_coord = row['xCoordinate']
            y_coord = row['yCoordinate']
            idx.insert(i, (x_coord, y_coord, x_coord, y_coord))

        # Iterate over each row in the DataFrame
        for i, row in self.border_df.iterrows():
            row_id = i
            feature_type = row['featureType']
            x_coord = row['xCoordinate']
            y_coord = row['yCoordinate']

            # Query the spatial index to find nearby points within the distance threshold
            nearby_points = list(idx.intersection((x_coord - distance_threshold, y_coord - distance_threshold, 
                                                   x_coord + distance_threshold, y_coord + distance_threshold)))

            # Filter neighbors based on distance and feature type, and ensure they are greater than the key
            points_to_add = sorted([j for j in nearby_points if j != row_id and 
                                    self.calc_distance(self.border_df.loc[j], self.border_df.loc[row_id]) <= distance_threshold and 
                                    self.border_df.loc[j, 'featureType'] != feature_type and j > row_id])
            # Store the nearby points in the dictionary
            star_neighbors[row_id] = points_to_add
            
        self.star_neighbors = star_neighbors