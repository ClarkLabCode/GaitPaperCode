#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""run_limb_umap.py:
    A script to embed limb timeseries data using UMAP's reference Python implementation.
    Requires that an HDF5 data file has been written to the appropriate location.
"""

import numpy as np
import os
import h5py
import umap
import time

# Entry point of script
if __name__ == "__main__":

    # Start an overall timer
    time_all = time.time()

    # Get the path to the data location
    base_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(base_path, "cache", "umap_data.h5")
    out_path = os.path.join(base_path, "cache", "umap_result.h5")

    # Import the raw data
    start_time = time.time()
    hf = h5py.File(file_path, "r")
    print(f"Importing data from {file_path}")

    # Check if data file contains the appropriately named data array
    if "x_l" in hf:
        x_l = np.array(hf["x_l"]).T
        print(f"Imported limb data array with shape {x_l.shape}")
    else:
        raise ValueError("Failed to import data")

    # Check if the target dimensionality is specified. If not, use 2
    if "limb_components" in hf:
        limb_components = int(hf["limb_components"].value)
    else:
        limb_components = 3

    # Check if the number of neighbors is specified. If not, use 30
    if "n_neighbors" in hf:
        n_neighbors = int(hf["n_neighbors"].value)
    else:
        n_neighbors = 30

    # Close the data file
    hf.close()
    print(f"Completed data import in {time.time() - start_time} seconds.")

    # Compute the embedding
    start_time = time.time()
    limb_reducer = umap.UMAP(n_neighbors=n_neighbors,
                             n_components=limb_components,
                             min_dist=0.0,
                             metric='euclidean'
                             )
    embedding_l = limb_reducer.fit_transform(x_l)
    print(f"Completed limb data embedding in {time.time() - start_time} seconds.")

    # Check for NaN values
    nnan_l = np.count_nonzero(np.isnan(embedding_l))
    print(f"Limb embedding contains {nnan_l} ({nnan_l / embedding_l.size} %) NaN values")

    # Write data to output file
    start_time = time.time()
    with h5py.File(out_path, "w") as hf:
        hf.create_dataset('embedding_l', dtype=np.float64, data=embedding_l)
    print(f"Wrote data to file {out_path} in {time.time() - start_time} seconds.")

    # Print a completion message
    print(f"Completed UMAP embeddings in {time.time() - time_all} seconds.")
