try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from kneed import KneeLocator
    import os
    from datetime import datetime
    import time
    import math
    from scipy.spatial.distance import euclidean
    from rtree import index
    from scipy.spatial import distance
    from multiprocessing import Pool
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Please install the missing packages using:")
    print("pip install pandas numpy matplotlib kneed rtree")
    sys.exit(1)


def read_data(directory):
    files = os.listdir(directory)
    csv_files = [file for file in files if file.endswith('.csv')]
    dataframes = []
    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        df = pd.read_csv(file_path)
        dataframes.append(df)
    combined_df = pd.concat(dataframes, ignore_index=True)
    combined_df = combined_df.sort_values(by='featureType', ignore_index=True)
    return combined_df

# Helper function to process each attack type in parallel
def process_attack_type(args):
    X, attack_types, val_to_exclude, max_k = args

    mask1 = attack_types == val_to_exclude
    mask2 = attack_types != val_to_exclude
    current_attack_X = X[mask1]  # X of the current attack type
    dataset = X[mask2]  # Dataset not including the current attack type

    idx = index.Index()  # Create a spatial index
    for i, point in enumerate(dataset):
        x_coord, y_coord = point[:2]
        idx.insert(i, (x_coord, y_coord, x_coord, y_coord))  # Insert each point as a bounding box

    all_distances = np.empty((0, max_k))  # Store the distances for the current attack type

    for point in current_attack_X:
        x_coord, y_coord = point[:2]
        nearest_ids = list(idx.nearest((x_coord, y_coord, x_coord, y_coord), max_k))

        if len(nearest_ids) > 0:
            nearest_points = np.array([dataset[nid][:2] for nid in nearest_ids])
            query_point = np.array([x_coord, y_coord])

            # Vectorized Euclidean distance calculation
            dists = np.linalg.norm(nearest_points - query_point, axis=1)
            sorted_dists = np.sort(dists)[:max_k]
            all_distances = np.concatenate((all_distances, [sorted_dists]), axis=0)
    return all_distances

# Main function to parallelize the outer loop
def rtree_processing(df):
    # x and y values of the original data
    X = np.array(df[['xCoordinate', 'yCoordinate']].values)
    # attack types of the original data
    attack_types = np.array(df['featureType'])

    num_pts = len(X)
    max_k = round(math.sqrt(num_pts)) + 1
    unique_attack_types = np.unique(attack_types)

    # Prepare the arguments for each process
    args_list = [(X, attack_types, val_to_exclude, max_k) for val_to_exclude in unique_attack_types]

    # Use multiprocessing to parallelize the process
    with Pool() as pool:
        results = pool.map(process_attack_type, args_list)

    # Combine the results from each process
    all_distances = np.vstack(results)
    return all_distances, max_k, num_pts

def memoization_table(all_distances, max_k, num_pts):
    dp_table = np.zeros((num_pts, max_k), dtype=float) #dynamic programming table
    # Summing distances for each k from 1 to 3
    dp_table[:, 2] = np.sum(all_distances[:, :3], axis=1)

    # Calculating cumulative sums for each k from 3 to max_k
    for k in range(3, max_k):
        dp_table[:, k] = dp_table[:, k - 1] + all_distances[:, k]

    column_sums = np.sum(dp_table, axis=0)
    averages = []
    for k in range(2, max_k):
        average = column_sums[k] / (num_pts * (k + 1))
        averages.append(average)       
    return averages

def knee_method(averages):
    k_value = 0    
    x = range(3, len(averages) + 3)
    y = averages    
    kn = KneeLocator(x, y, curve='concave', direction='increasing')
    averages = np.array(averages)
    if kn.knee == 3 or kn.knee == len(averages) or kn.knee == None:
        # Calculate differences between neighboring points
        differences = averages[1:] - averages[:-1]

        for i in range(len(differences)):
            if i > 0:
                if temp < differences[i]:
                    k_value = i
                    distance_threshold = averages[i - 1]
                    break
            temp = differences[i]     
    else:
        k_value = kn.knee
        distance_threshold = averages[k_value-3]
        
    return distance_threshold

def main():
    directory = '/home/amk7r/colocation_mining/Parallel-Project/data/north_america/'
    df = read_data(directory)
    s = time.time()
    all_distances, max_k, num_pts = rtree_processing(df)
    averages = memoization_table(all_distances, max_k, num_pts)
    e = time.time()
    distance_threshold = knee_method(averages)

    print('DISTANCE THRESHOLD (KM):', distance_threshold)
    with open('IntermediateData/distance_threshold_parameter.txt', 'w') as file:
        file.write(str(distance_threshold))
    print('TIME TAKEN (SEC):', e - s)
    
if __name__ == "__main__":
    main()
