import numpy as np
import h5py
import umap
import time


# Entry point of script
if __name__ == "__main__":
    # Start an overall timer
    time_all = time.time()

    # Load the file paths from the configuration .txt file
    with open('config.txt', 'r') as f:
        cfgstr = f.read()
    [sourcepath, resultpath] = cfgstr.split(',')

    # Import the raw data
    start_time = time.time()
    with h5py.File(sourcepath, "r") as hf:
        if "x_v" in hf:
            x_v = np.array(hf["x_v"]).T
            centroid_components = int(hf["centroid_components"][()])
            print(f"Imported centroid data array with shape {x_v.shape}")
            embed_centroids = True
        else:
            embed_centroids = False

        if "x_l" in hf:
            x_l = np.array(hf["x_l"]).T
            limb_components = int(hf["limb_components"][()])
            print(f"Imported limb data array with shape {x_l.shape}")
            embed_limbs = True
        else:
            embed_limbs = False

    if not embed_centroids and not embed_limbs:
        raise ValueError("No data imported!")
    print(f"Completed data import in {time.time() - start_time} seconds.")

    # Compute embeddings
    n_neighbors = 30
    if embed_centroids:
        start_time = time.time()
        centroid_reducer = umap.UMAP(n_neighbors=n_neighbors,
                                     n_components=centroid_components,
                                     min_dist=0.0,
                                     metric='euclidean'
                                     )
        embedding_v = centroid_reducer.fit_transform(x_v)
        print(f"Completed centroid data embedding in {time.time() - start_time} seconds.")

        # Check for NaN values
        nnan_v = np.count_nonzero(np.isnan(embedding_v))
        print(f"Centroid embedding contains {nnan_v} ({nnan_v / embedding_v.size} %) NaN values")

    if embed_limbs:
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

    # Write data to file
    start_time = time.time()
    with h5py.File(resultpath, "w") as hf:
        if embed_centroids:
            hf.create_dataset('embedding_v', dtype=np.float64, data=embedding_v)
        if embed_limbs:
            hf.create_dataset('embedding_l', dtype=np.float64, data=embedding_l)
    print(f"Wrote data to file in {time.time() - start_time} seconds.")

    # Print a completion message
    print(f"Completed UMAP embeddings in {time.time() - time_all} seconds.")
